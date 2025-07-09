"""
    makeSbus(baseMVA, bus, gen, Vm, load::Matrix{Float64}, pvarray; dc=false, Sg=nothing, return_derivative=false)

Build the vector of complex bus power injections, including PV array contributions.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `gen`: Generator data matrix with columns representing generator parameters
- `Vm`: Vector of bus voltage magnitudes
- `load`: Load data matrix with columns representing load parameters
- `pvarray`: PV array data matrix with columns representing PV parameters

# Keyword Arguments
- `dc`: Boolean indicating whether to use DC power flow assumptions (default: false)
- `Sg`: Optional pre-computed generator complex power injections (default: nothing)
- `return_derivative`: Boolean indicating whether to return derivative of Sbus with respect to Vm (default: false)

# Returns
- If `return_derivative=false`: Vector of complex bus power injections (Sbus)
- If `return_derivative=true`: Sparse matrix of partial derivatives of power injections with respect to voltage magnitude (dSbus_dVm)

# Description
This function computes the vector of complex bus power injections (Sbus) for power flow analysis.
It accounts for ZIP load models, generator injections, and PV array contributions.

When `return_derivative=true`, it returns the partial derivatives of the power injections with respect to voltage magnitude,
which is useful for power flow Jacobian calculations.
"""
function makeSbus(baseMVA, bus, gen, Vm, load::Matrix{Float64}, pvarray; dc=false, Sg=nothing, return_derivative=false)
    nb = size(bus, 1)
    pw_1=zeros(size(bus,1),1)
    pw_2=zeros(size(bus,1),1)
    pw_3=zeros(size(bus,1),1)
    pw_1[Int64.(load[:,LOAD_CND])]=load[:,LOADP_PERCENT]
    pw_2[Int64.(load[:,LOAD_CND])]=load[:,LOADI_PERCENT]
    pw_3[Int64.(load[:,LOAD_CND])]=load[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = makeSdzip(baseMVA, bus,pw_1,pw_2,pw_3)

     # Initialize PV power and derivatives
    Spv = zeros(Complex{Float64}, nb)
    dSpv_dVm = spzeros(nb, nb)
    
    # If there are PV arrays and not using DC power flow, calculate their power and derivatives
    if !isnothing(pvarray)
        Spv, dSpv_dVm = calculate_pv_power(pvarray, bus, Vm, baseMVA)
    end

    if return_derivative
        if isempty(Vm)
            dSbus_dVm = spzeros(nb, nb)
        else
            # Load derivative + PV derivative
            dSbus_dVm = -(spdiagm(0 => Sd.i + 2 .* Vm .* Sd.z)) + dSpv_dVm
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = sparse(gbus, 1:ngon, 1, nb, ngon)  # connection matrix
        # element i, j is 1 if gen on(j) at bus i is ON
        if Sg !== nothing
            Sbusg = Cg * Sg[on]
        else
            Sbusg = Cg * (gen[on, PG] .+ 1im * gen[on, QG]) / baseMVA
        end

        if dc
            Vm = ones(nb,1)
        end
        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm .+ Sd.z .* Vm.^2

       # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads + PV power)
        Sbus = Sbusg - Sbusd + Spv
        return Sbus
    end
end

"""
    makeSbus(baseMVA, bus, gen, Vm, Sg=nothing, return_derivative=false)

Build the vector of complex bus power injections (simplified version without PV arrays).

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `gen`: Generator data matrix with columns representing generator parameters
- `Vm`: Vector of bus voltage magnitudes

# Keyword Arguments
- `Sg`: Optional pre-computed generator complex power injections (default: nothing)
- `return_derivative`: Boolean indicating whether to return derivative of Sbus with respect to Vm (default: false)

# Returns
- If `return_derivative=false`: Vector of complex bus power injections (Sbus)
- If `return_derivative=true`: Sparse matrix of partial derivatives of power injections with respect to voltage magnitude (dSbus_dVm)

# Description
This is a simplified version of the makeSbus function that does not include PV array contributions.
It computes the vector of complex bus power injections (Sbus) for power flow analysis,
accounting for ZIP load models and generator injections.
"""
function makeSbus(baseMVA, bus, gen, Vm, Sg=nothing, return_derivative=false)
    nb = size(bus, 1)

    # Get load parameters
    Sd = makeSdzip(baseMVA, bus)

    if return_derivative
        if isempty(Vm)
            dSbus_dVm = spzeros(nb, nb)
        else
            dSbus_dVm = -(spdiagm(0 => Sd.i + 2 .* Vm .* Sd.z))
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = sparse(gbus, 1:ngon, 1, nb, ngon)  # connection matrix
        # element i, j is 1 if gen on(j) at bus i is ON
        if !isnothing(Sg)
            Sbusg = Cg * Sg[on]
        else
            Sbusg = Cg * (gen[on, PG] .+ 1im * gen[on, QG]) / baseMVA
        end

        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm .+ Sd.z .* Vm.^2

        # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads)
        Sbus = Sbusg - Sbusd
        return Sbus
    end
end

"""
    calculate_pv_power(pvarray, bus, Vm, baseMVA)

Calculate power injection from PV arrays and its derivative with respect to voltage magnitude.

# Arguments
- `pvarray`: PV array data matrix with columns representing PV parameters
- `bus`: Bus data matrix with columns representing bus parameters
- `Vm`: Vector of bus voltage magnitudes
- `baseMVA`: Base MVA for the system

# Returns
- `Spv`: Vector of complex power injections from PV arrays
- `dSpv_dVm`: Sparse matrix of partial derivatives of PV power injections with respect to voltage magnitude

# Description
This function calculates the power injection from PV arrays based on a modified power function model
that accounts for the PV array characteristics including open circuit voltage, short circuit current,
and maximum power point voltage. It also computes the derivatives of these injections with respect
to bus voltage magnitudes, which are needed for power flow Jacobian calculations.

The PV model uses a modified power function with correction terms to better represent the
current-voltage characteristics of PV panels:
I = Isc * (1 - (V/Voc)^a)^b * (1 - c * ((V/Vmpp) - 1)^2)

# Notes
- Only active PV arrays (PV_IN_SERVICE > 0) are considered
- The function maps PV arrays to their respective buses using the bus numbering
- Power output is converted to per-unit on system MVA base
- Only real power (no reactive power) is considered for PV arrays
"""
function calculate_pv_power(pvarray, bus, Vm, baseMVA)
    
    nb = size(Vm, 1)
    Spv = zeros(Complex{Float64}, nb)
    dSpv_dVm = spzeros(nb, nb)
    
    # Filter PV arrays that are in service
    active_pv = pvarray[pvarray[:, PV_IN_SERVICE] .> 0, :]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # Create mapping from bus numbers to indices
    bus_num_to_idx = Dict{Int, Int}()
    for i in eachindex(bus[:,1])
        bus_num_to_idx[Int(bus[i, BUS_I])] = i
    end
    
    # Get bus indices for PV connections
    bus_nums = Int.(active_pv[:, PV_BUS])
    valid_pv = falses(size(active_pv, 1))
    bus_indices = zeros(Int, size(active_pv, 1))
    
    for i in eachindex(bus_nums)
        if haskey(bus_num_to_idx, bus_nums[i])
            valid_pv[i] = true
            bus_indices[i] = bus_num_to_idx[bus_nums[i]]
        end
    end
    
    # Keep only valid PV arrays
    active_pv = active_pv[valid_pv, :]
    bus_indices = bus_indices[valid_pv]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # Get base voltage for buses (kV)
    base_kvs = bus[bus_indices, BASE_KV]
    
    # Get current bus voltages (per unit)
    V_buses = Vm[bus_indices]
    
    # Modified power function model parameters
    a = 10.000000
    b = 0.547596
    c = 0.023812
    
    # Extract PV parameters
    Vocs = active_pv[:, PV_VOC]      # Open circuit voltage (V)
    Iscs = active_pv[:, PV_ISC]      # Short circuit current (A)
    areas = active_pv[:, PV_AREA]    # Area or other scaling factor
    Vmpps = active_pv[:, PV_VMPP]    # Maximum power point voltage (V), ensure this column exists
    
    # Convert per unit voltage to actual voltage (V)
    voltage_ratios = 1
    V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # Convert to volts
    
    # Initialize current and derivative arrays
    I_arrays = zeros(size(active_pv, 1))
    dI_dVs = zeros(size(active_pv, 1))
    
    # Create masks for vectorized operations
    neg_v_mask = V_arrays .< 0
    valid_v_mask = (V_arrays .>= 0) .& (V_arrays .<= Vocs)
    over_v_mask = V_arrays .> Vocs
    
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
            
            # Ensure current is not negative
            I_arrays[i] = max(0, I_arrays[i])
        end
    end
    
    # Calculate derivatives (for valid voltage range)
    for i in eachindex(V_arrays)
        if valid_v_mask[i] && V_arrays[i] > 0
            # Derivative of base power function term
            d_base_term_dv = -a * b * (1 - v_ratios[i]^a)^(b-1) * v_ratios[i]^(a-1) / Vocs[i]
            
            # Derivative of correction term
            d_correction_term_dv = -c * 2 * vmpp_ratios[i] / Vmpps[i]
            
            # Use product rule to calculate total derivative
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
    
    # Update power injection and derivative matrices
    for i in eachindex(bus_indices)
        Spv[bus_indices[i]] += P_pus[i] + 0im  # Real power only
        dSpv_dVm[bus_indices[i], bus_indices[i]] += dP_dV_pus[i]
    end
    
    return Spv, dSpv_dVm
end


# function calculate_pv_power(pvarray, bus, Vm, baseMVA)
    
#     nb = size(Vm, 1)
#     Spv = zeros(Complex{Float64}, nb)
#     dSpv_dVm = spzeros(nb, nb)
    
#     # Filter PV arrays that are in service
#     active_pv = pvarray[pvarray[:, PV_IN_SERVICE] .> 0, :]
    
#     if isempty(active_pv)
#         return Spv, dSpv_dVm
#     end
    
#     # Create mapping from bus numbers to indices
#     bus_num_to_idx = Dict{Int, Int}()
#     for i in 1:size(bus, 1)
#         bus_num_to_idx[Int(bus[i, BUS_I])] = i
#     end
    
#     # Get bus indices for PV connections
#     bus_nums = Int.(active_pv[:, PV_BUS])
#     valid_pv = falses(size(active_pv, 1))
#     bus_indices = zeros(Int, size(active_pv, 1))
    
#     for i in 1:length(bus_nums)
#         if haskey(bus_num_to_idx, bus_nums[i])
#             valid_pv[i] = true
#             bus_indices[i] = bus_num_to_idx[bus_nums[i]]
#         end
#     end
    
#     # Keep only valid PV arrays
#     active_pv = active_pv[valid_pv, :]
#     bus_indices = bus_indices[valid_pv]
    
#     if isempty(active_pv)
#         return Spv, dSpv_dVm
#     end
    
#     # Get base voltage for buses (kV)
#     base_kvs = bus[bus_indices, BASE_KV]
    
#     # Get current bus voltages (per unit)
#     V_buses = Vm[bus_indices]
    
#     # Empirical formula parameters
#     a = 14.5
    
#     # Extract PV parameters
#     Vocs = active_pv[:, PV_VOC]      # Open circuit voltage (V)
#     Iscs = active_pv[:, PV_ISC]      # Short circuit current (A)
#     areas = active_pv[:, PV_AREA]    # Area or other scaling factor
    
#     # Convert per unit voltage to actual voltage (V)
#     voltage_ratios = 1
#     V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # Convert to volts
    
#     # Initialize current and derivative arrays
#     I_arrays = zeros(size(active_pv, 1))
#     dI_dVs = zeros(size(active_pv, 1))
    
#     # Create masks for vectorized operations
#     neg_v_mask = V_arrays .< 0
#     valid_v_mask = (V_arrays .>= 0) .& (V_arrays .<= Vocs)
#     over_v_mask = V_arrays .> Vocs
    
#     # Calculate ratios
#     ratios = zeros(size(V_arrays))
#     ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vocs[valid_v_mask]
    
#     # Calculate current using empirical formula (vectorized operation)
#     I_arrays[valid_v_mask] = Iscs[valid_v_mask] .* areas[valid_v_mask] .* (1 .- ratios[valid_v_mask].^a)
    
#     # Calculate derivatives (vectorized operation)
#     # Avoid derivative issues at V=0
#     nonzero_v = (V_arrays .> 0) .& valid_v_mask
#     dI_dVs[nonzero_v] = -Iscs[nonzero_v] .* areas[nonzero_v] .* a .* ratios[nonzero_v].^(a-1) ./ Vocs[nonzero_v]
    
#     # Calculate power (W)
#     P_arrays = V_arrays .* I_arrays
    
#     # Convert to per unit (considering baseMVA)
#     P_pus = P_arrays ./ (baseMVA * 1e6)  # Convert from W to per unit
    
#     # Calculate power derivative with respect to voltage
#     dP_dVs = I_arrays .+ V_arrays .* dI_dVs
#     dP_dV_pus = (dP_dVs .* base_kvs .* 1000 .* voltage_ratios) ./ (baseMVA * 1e6)
    
#     # Update power injection and derivative matrices
#     for i in 1:length(bus_indices)
#         Spv[bus_indices[i]] += P_pus[i] + 0im  # Real power only
#         dSpv_dVm[bus_indices[i], bus_indices[i]] += dP_dV_pus[i]
#     end
    
#     return Spv, dSpv_dVm
# end
