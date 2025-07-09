"""
    makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu; dc=false, Sg=nothing, return_derivative=false)

Build the vector of complex bus power injections using GPU acceleration.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `gen_gpu`: Generator data matrix on GPU with columns representing generator parameters
- `gen`: Generator data matrix on CPU with columns representing generator parameters
- `Vm_gpu`: Vector of bus voltage magnitudes on GPU
- `load_gpu`: Load data matrix on GPU with columns representing load parameters

# Keyword Arguments
- `dc`: Boolean indicating whether to use DC power flow assumptions (default: false)
- `Sg`: Optional pre-computed generator complex power injections (default: nothing)
- `return_derivative`: Boolean indicating whether to return derivative of Sbus with respect to Vm (default: false)

# Returns
- If `return_derivative=false`: Vector of complex bus power injections (Sbus)
- If `return_derivative=true`: Sparse matrix of partial derivatives of power injections with respect to voltage magnitude (dSbus_dVm)

# Description
This function computes the vector of complex bus power injections (Sbus) for power flow analysis using GPU acceleration.
It accounts for ZIP load models (constant power, constant current, and constant impedance components) and generator injections.

When `return_derivative=true`, it returns the partial derivatives of the power injections with respect to voltage magnitude,
which is useful for power flow Jacobian calculations.

# Notes
- All power values are converted to per-unit on system MVA base
- The function handles ZIP load models with percentages specified in load_gpu
- Generator status is considered when computing injections
- When dc=true, voltage magnitudes are set to 1.0 p.u.

# Constants Used (assumed to be defined elsewhere)
- LOAD_CND: Column index for load bus number in load_gpu matrix
- LOADP_PERCENT: Column index for constant power percentage in load_gpu matrix
- LOADI_PERCENT: Column index for constant current percentage in load_gpu matrix
- LOADZ_PERCENT: Column index for constant impedance percentage in load_gpu matrix
- GEN_STATUS: Column index for generator status in gen matrix
- GEN_BUS: Column index for generator bus number in gen matrix
- PG: Column index for real power output in gen_gpu matrix
- QG: Column index for reactive power output in gen_gpu matrix
"""
function makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu; dc=false, Sg=nothing, return_derivative=false)

    nb = size(bus_gpu, 1)
    pw_1=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_2=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_3=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_1[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADP_PERCENT]
    pw_2[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADI_PERCENT]
    pw_3[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = makeSdzip_gpu(baseMVA, bus_gpu,pw_1,pw_2,pw_3)

    if return_derivative
        if isempty(Vm_gpu)
            dSbus_dVm = PowerFlow.CUDA.spzeros(nb, nb)
        else
            diag_elements = Sd.i + 2 .* Vm_gpu .* Sd.z
            dSbus_dVm = -PowerFlow.Diagonal(diag_elements)
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
        Sbus = Sbusg - Sbusd
        return Sbus
    end
end
