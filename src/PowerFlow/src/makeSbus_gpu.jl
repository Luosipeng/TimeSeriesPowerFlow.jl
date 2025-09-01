"""
    makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu; dc=false, Sg=nothing, return_derivative=false)

Build the vector of complex bus power injections using GPU acceleration.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `gen_gpu`: Generator data matrix on GPU with columns representing generator parameters
- `gen`: Generator data matrix on CPU with columns representing generator parameters
- `Vm_gpu`: Vector of bus voltage magnitudes on GPU
- `load_gpu`: Load data matrix on GPU with columns representing load parameters
- `pvarray_gpu`: PV array data matrix on GPU with columns representing PV parameters

# Keyword Arguments
- `dc`: Boolean indicating whether to use DC power flow assumptions (default: false)
- `Sg`: Optional pre-computed generator complex power injections (default: nothing)
- `return_derivative`: Boolean indicating whether to return derivative of Sbus with respect to Vm (default: false)

# Returns
- If `return_derivative=false`: Vector of complex bus power injections (Sbus)
- If `return_derivative=true`: Sparse matrix of partial derivatives of power injections with respect to voltage magnitude (dSbus_dVm)

# Description
This function computes the vector of complex bus power injections (Sbus) for power flow analysis using GPU acceleration.
It accounts for ZIP load models (constant power, constant current, and constant impedance components), generator injections,
and PV array injections.

When `return_derivative=true`, it returns the partial derivatives of the power injections with respect to voltage magnitude,
which is useful for power flow Jacobian calculations.

# Notes
- All power values are converted to per-unit on system MVA base
- The function handles ZIP load models with percentages specified in load_gpu
- Generator status is considered when computing injections
- When dc=true, voltage magnitudes are set to 1.0 p.u.
- PV array injections are calculated separately and added to the total bus injections

# Constants Used (assumed to be defined elsewhere)
- LOAD_CND: Column index for load bus number in load_gpu matrix
- LOADP_PERCENT: Column index for constant power percentage in load_gpu matrix
- LOADI_PERCENT: Column index for constant current percentage in load_gpu matrix
- LOADZ_PERCENT: Column index for constant impedance percentage in load_gpu matrix
- GEN_STATUS: Column index for generator status in gen matrix
- GEN_BUS: Column index for generator bus number in gen matrix
- PG: Column index for real power output in gen_gpu matrix
- QG: Column index for reactive power output in gen_gpu matrix
- PV_IN_SERVICE: Column index for PV array service status in pvarray_gpu matrix
"""
function makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu; dc=false, Sg=nothing, return_derivative=false)

    nb = size(bus_gpu, 1)
    pw_1=CUDA.zeros(size(bus_gpu,1),1)
    pw_2=CUDA.zeros(size(bus_gpu,1),1)
    pw_3=CUDA.zeros(size(bus_gpu,1),1)
    pw_1[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADP_PERCENT]
    pw_2[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADI_PERCENT]
    pw_3[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = PowerFlow.makeSdzip_gpu(baseMVA, bus_gpu,pw_1,pw_2,pw_3)
    dSpv_dVm_gpu = CUDA.CUSPARSE.spzeros(Complex{Float64}, nb)

    # If there are PV arrays and not using DC power flow, calculate their power and derivatives
    if !isnothing(pvarray_gpu)
        Spv, dSpv_dVm = PowerFlow.calculate_pv_power_gpu(pvarray_gpu, bus_gpu, Vm_gpu, baseMVA)
    end

    if return_derivative
        if isempty(Vm_gpu)
            dSbus_dVm = PowerFlow.CUDA.spzeros(nb, nb)
        else
            dSpv_dVm_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(dSpv_dVm)
            diag_elements = Sd.i + 2 .* Vm_gpu .* Sd.z
            dSbus_dVm = -PowerFlow.Diagonal(diag_elements) + dSpv_dVm_gpu
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = CUSPARSE.CuSparseMatrixCSR(sparse(Int64.(gbus), collect(1:ngon), ones(ngon), nb, ngon))

        # element i, j is 1 if gen on(j) at bus i is ON
        if Sg !== nothing
            Sbusg = Cg * Sg[on]
        else
            # Step 1: Create generator complex power vector
            Sg = gen_gpu[on, PG] .+ 1im * gen_gpu[on, QG]

            # Step 2: Create result vector (bus injection power)
            Sbusg = CUDA.zeros(ComplexF64, size(Cg, 1))

            # Step 3: Use CUSPARSE.mv! function to perform matrix-vector multiplication
            # Add extra character parameter 'O' to indicate operation type
            CUDA.CUSPARSE.mv!('N', one(ComplexF64), Cg, Sg, zero(ComplexF64), Sbusg, 'O')

            # Step 4: Divide by base power baseMVA for per-unit normalization
            Sbusg = Sbusg ./ baseMVA
        end

        if dc
            Vm = PowerFlow.CUDA.ones(nb,1)
        end
        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm_gpu .+ Sd.z .* Vm_gpu.^2

        # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads)
        Sbus = Sbusg - Sbusd + Spv
        return Sbus
    end
end

"""
    calculate_pv_power_gpu(pvarray_gpu, bus_gpu, Vm_gpu, baseMVA)

Calculate power injections from PV arrays and their derivatives with respect to voltage magnitude using GPU acceleration.

# Arguments
- `pvarray_gpu`: PV array data matrix on GPU with columns representing PV parameters
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `Vm_gpu`: Vector of bus voltage magnitudes on GPU
- `baseMVA`: Base MVA for the system

# Returns
- `Spv`: Vector of complex power injections from PV arrays
- `dSpv_dVm`: Sparse matrix of partial derivatives of PV power injections with respect to voltage magnitude

# Description
This function computes the power injections from PV arrays using a modified power function model
that relates voltage to current output. It also calculates the derivatives of these injections
with respect to voltage magnitude for use in power flow Jacobian calculations.

The function filters active PV arrays, maps them to their connected buses, and calculates
their power output based on the current voltage conditions.

# Notes
- Power output is calculated using a modified power function model with parameters a, b, and c
- The model accounts for the nonlinear relationship between voltage and current in PV arrays
- Results are converted to per-unit on system MVA base
- Only real power is considered (no reactive power from PV arrays)
- Calculations are partially performed on CPU for complex operations

# Constants Used (assumed to be defined elsewhere)
- PV_IN_SERVICE: Column index for PV array service status
- BUS_I: Column index for bus number in bus_gpu matrix
- PV_BUS: Column index for connected bus number in pvarray_gpu matrix
- BASE_KV: Column index for base voltage in bus_gpu matrix
- PV_VOC: Column index for open circuit voltage in pvarray_gpu matrix
- PV_ISC: Column index for short circuit current in pvarray_gpu matrix
- PV_AREA: Column index for PV array area in pvarray_gpu matrix
- PV_VMPP: Column index for maximum power point voltage in pvarray_gpu matrix
"""
function calculate_pv_power_gpu(pvarray_gpu, bus_gpu, Vm_gpu, baseMVA)
    nb = size(Vm_gpu, 1)
    Spv = CUDA.zeros(Complex{Float64}, nb)
    dSpv_dVm = CUDA.CUSPARSE.spzeros(nb, nb)
    
    # Filter PV arrays that are in service
    active_pv = pvarray_gpu[pvarray_gpu[:, PV_IN_SERVICE] .> 0, :]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # Create mapping from bus number to index
    valid_pv, bus_indices = process_pv_connections(bus_gpu, active_pv, BUS_I, PV_BUS)
    
    # Keep only valid PV arrays
    active_pv = active_pv[valid_pv, :]
    bus_indices = bus_indices[valid_pv]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # Transfer necessary data back to CPU for processing
    bus_cpu = Array(bus_gpu)
    Vm_cpu = Array(Vm_gpu)
    active_pv_cpu = Array(active_pv)
    
    # Get base voltage of buses (kV)
    base_kvs = bus_cpu[bus_indices, BASE_KV]
    
    # Get current bus voltage (per unit)
    V_buses = Vm_cpu[bus_indices]
    
    # Modified power function model parameters
    a = 10.000000
    b = 0.547596
    c = 0.023812
    
    # Extract PV parameters
    Vocs = active_pv_cpu[:, PV_VOC]      # Open circuit voltage (V)
    Iscs = active_pv_cpu[:, PV_ISC]      # Short circuit current (A)
    areas = active_pv_cpu[:, PV_AREA]    # Area or other scaling factor
    Vmpps = active_pv_cpu[:, PV_VMPP]    # Maximum power point voltage (V)
    
    # Convert per unit voltage to actual voltage (V)
    voltage_ratios = 1
    V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # Convert to volts
    
    # Initialize current and derivative arrays
    I_arrays = zeros(size(active_pv_cpu, 1))
    dI_dVs = zeros(size(active_pv_cpu, 1))
    
    # Create masks for vectorized operations
    neg_v_mask = real.(V_arrays) .< 0
    valid_v_mask = (real.(V_arrays) .>= 0) .& (real.(V_arrays) .<= Vocs)
    over_v_mask = real.(V_arrays) .> Vocs
    
    # Calculate ratios
    v_ratios = zeros(size(V_arrays))
    vmpp_ratios = zeros(size(V_arrays))
    v_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vocs[valid_v_mask]
    vmpp_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vmpps[valid_v_mask] .- 1
    
    # Calculate current using modified power function model
    for i in 1:length(V_arrays)
        if valid_v_mask[i]
            # Base power function term
            base_term = (1 - v_ratios[i]^a)^b
            
            # Correction term
            correction_term = (1 - c * vmpp_ratios[i]^2)
            
            # Calculate current
            I_arrays[i] = Iscs[i] * base_term * correction_term
            
            # Ensure current is non-negative
            I_arrays[i] = max(0, I_arrays[i])
        end
    end
    
    # Calculate derivatives (for valid voltage range)
    for i in eachindex(V_arrays)
        if valid_v_mask[i] && real(V_arrays[i]) > 0
            # Derivative of base power function term
            d_base_term_dv = -a * b * (1 - v_ratios[i]^a)^(b-1) * v_ratios[i]^(a-1) / Vocs[i]
            
            # Derivative of correction term
            d_correction_term_dv = -c * 2 * vmpp_ratios[i] / Vmpps[i]
            
            # Calculate total derivative using product rule
            base_term = (1 - v_ratios[i]^a)^b
            correction_term = (1 - c * vmpp_ratios[i]^2)
            
            dI_dVs[i] = Iscs[i] * (
                d_base_term_dv * correction_term + 
                base_term * d_correction_term_dv
            )
            
            # If current is 0, derivative should also be 0
            if I_arrays[i] <= 0
                dI_dVs[i] = 0
            end
        end
    end
    
    # Calculate power (W)
    P_arrays = V_arrays .* I_arrays
    
    # Convert to per unit (considering baseMVA)
    P_pus = P_arrays ./ (baseMVA * 1e6)  # Convert from W to per unit
    
    # Calculate power derivative with respect to voltage
    dP_dVs = I_arrays .+ V_arrays .* dI_dVs
    dP_dV_pus = (dP_dVs .* base_kvs .* 1000 .* voltage_ratios) ./ (baseMVA * 1e6)
    
    # Update power injection and derivative matrix
    Spv_cpu = Array(Spv)
    for i in eachindex(bus_indices)
        Spv_cpu[bus_indices[i]] += P_pus[i] + 0im  # Real power only
    end
    
    # Transfer results back to GPU
    copyto!(Spv, Spv_cpu)
    
    # Process sparse derivative matrix
    I_indices = bus_indices
    J_indices = bus_indices
    V_values = dP_dV_pus
    
    # Create new sparse matrix
    dSpv_dVm = CUDA.CUSPARSE.sparse(I_indices, J_indices, V_values, nb, nb)
    
    return Spv, dSpv_dVm
end

"""
    create_bus_mapping(bus_gpu, BUS_I)

Create a mapping dictionary from bus numbers to their indices in the bus matrix.

# Arguments
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `BUS_I`: Column index for bus number in bus_gpu matrix

# Returns
- `bus_to_idx`: Dictionary mapping bus numbers to their indices in the bus matrix

# Description
This function creates a dictionary that maps bus identification numbers to their 
corresponding indices in the bus matrix. This mapping is useful for quickly finding
the position of a specific bus in the data structures.

# Notes
- The function transfers bus data from GPU to CPU for processing
- Bus numbers may not be sequential or start from 1, hence the need for this mapping
"""
function create_bus_mapping(bus_gpu, BUS_I)
    # Transfer bus data back to CPU to process mapping relationships
    bus_ids = Array(Int.(bus_gpu[:, BUS_I]))
    n = length(bus_ids)
    
    # Create mapping dictionary
    bus_to_idx = Dict{Int, Int}()
    for i in 1:n
        bus_to_idx[bus_ids[i]] = i
    end
    
    return bus_to_idx
end

"""
    find_bus_indices(bus_to_idx, active_pv, PV_BUS)

Find the indices of buses to which PV arrays are connected.

# Arguments
- `bus_to_idx`: Dictionary mapping bus numbers to their indices in the bus matrix
- `active_pv`: PV array data matrix with columns representing PV parameters
- `PV_BUS`: Column index for connected bus number in active_pv matrix

# Returns
- `valid_pv`: Boolean array indicating which PV arrays are connected to valid buses
- `bus_indices`: Array of bus indices corresponding to each valid PV array

# Description
This function identifies which buses in the system have PV arrays connected to them.
It returns a boolean mask indicating which PV arrays are connected to valid buses,
and an array of the corresponding bus indices for those valid connections.

# Notes
- The function processes data entirely on CPU for better dictionary lookup performance
- PV arrays connected to non-existent buses are marked as invalid
"""
function find_bus_indices(bus_to_idx, active_pv, PV_BUS)
    # Transfer PV bus numbers back to CPU
    bus_nums = Array(Int.(active_pv[:, PV_BUS]))
    n = length(bus_nums)
    
    # Create result arrays
    valid_pv = falses(n)
    bus_indices = zeros(Int, n)
    
    # Process each bus number
    for i in 1:n
        bus_id = bus_nums[i]
        if haskey(bus_to_idx, bus_id)
            valid_pv[i] = true
            bus_indices[i] = bus_to_idx[bus_id]
        end
    end
    
    return valid_pv, bus_indices
end

"""
    process_pv_connections(bus_gpu, active_pv, BUS_I, PV_BUS)

Process the connections between PV arrays and buses in the power system.

# Arguments
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `active_pv`: PV array data matrix with columns representing PV parameters
- `BUS_I`: Column index for bus number in bus_gpu matrix
- `PV_BUS`: Column index for connected bus number in active_pv matrix

# Returns
- `valid_pv`: Boolean array indicating which PV arrays are connected to valid buses
- `bus_indices`: Array of bus indices corresponding to each valid PV array

# Description
This function coordinates the process of mapping PV arrays to their connected buses
in the power system. It creates a mapping from bus numbers to indices, then uses this
mapping to identify which PV arrays are connected to valid buses in the system.

# Notes
- This is a wrapper function that calls create_bus_mapping and find_bus_indices
- The function handles the complete process of validating PV-to-bus connections
"""
function process_pv_connections(bus_gpu, active_pv, BUS_I, PV_BUS)
    # Create mapping dictionary
    bus_to_idx = create_bus_mapping(bus_gpu, BUS_I)
    
    # Find bus indices
    valid_pv, bus_indices = find_bus_indices(bus_to_idx, active_pv, PV_BUS)
    
    return valid_pv, bus_indices
end
