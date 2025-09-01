"""
    newtondcpf_gpu(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")

Solve DC power flow using Newton's method with GPU acceleration.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `gen`: Generator data matrix with columns representing generator parameters
- `load`: Load data matrix with columns representing load parameters
- `pvarray`: PV array data matrix with columns representing PV parameters
- `Ybus`: Admittance matrix of the power system
- `V0`: Initial voltage vector
- `ref`: Index of reference bus
- `p`: Vector of indices of PV and PQ buses (non-reference buses)
- `tol0`: Convergence tolerance for the solution
- `max_it0`: Maximum number of iterations
- `alg`: Linear solver algorithm (default: "")

# Returns
- `V`: Final voltage vector (real part only, as this is DC power flow)
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed

# Description
This function solves the DC power flow problem using Newton's method with GPU acceleration.
It iteratively updates the voltage angles to minimize power mismatches at non-reference buses.
The function transfers all necessary data to the GPU for faster computation.

The algorithm follows these steps:
1. Initialize voltage vector and transfer data to GPU
2. Calculate initial power mismatches
3. For each iteration:
   a. Compute the Jacobian matrix
   b. Solve the linear system to find the voltage update
   c. Update the voltage vector
   d. Check for convergence based on power mismatches

# Notes
- This implementation is specifically for DC power flow, focusing on real power and voltage angles
- All matrices and vectors are transferred to GPU for accelerated computation
- The function handles PV arrays as additional power injections
- The Jacobian accounts for voltage-dependent loads through the makeSbus_gpu function
- The solution returns only the real part of the voltage vector

# Related Functions
- `makeSbus_gpu`: Calculates power injections at buses
- `dSbus_dV`: Calculates derivatives of power injections with respect to voltage
- `julinsolve`: Solves the linear system for the Newton update
"""
function newtondcpf_gpu(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0 .+ im*0.0
    nb = length(V)
    
    # Transfer to GPU
    V_gpu = CUDA.CuVector(V)

    # Transfer Matrix to GPU
    Ybus_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(Ybus)
    bus_gpu = PowerFlow.CuArray(bus)
    gen_gpu = PowerFlow.CuArray(gen)
    load_gpu = PowerFlow.CuArray(load)
    pvarray_gpu = CuArray(pvarray)
    
    # Set up indexing for updating V
    np = length(p)
    j1 = 1; j2 = np; # j1:j2 - V angle of pv buses

    # Create indexing matrix
    Cp_index =sparse(1:np, p, ones(ComplexF64,np), np, nb);
    Cp_index=CuSparseMatrixCSR(Cp_index)
    Cp_index_transpose = sparse( p,1:np, ones(ComplexF64,np), nb,np)
    Cp_index_transpose=CuSparseMatrixCSR(Cp_index_transpose)

    # Evaluate F(x0)
    mis = V_gpu .* conj.(Ybus_gpu * V_gpu) - PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, V_gpu, load_gpu, pvarray_gpu)
    F = real(mis[p])
     # Check tolerance
    normF = norm(F, Inf)
    if normF < tol
        converged = true
    end
    # Do Newton iterations
    while (!converged && i < max_it)

        # Update iteration counter
        i += 1

        # Evaluate Jacobian
        dSbus_dVa, dSbus_dVm = PowerFlow.dSbus_dV(Ybus_gpu, V_gpu)
        neg_dSd_dVm = PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, V_gpu, load_gpu, pvarray_gpu, return_derivative=true)
        dSbus_dVm -= neg_dSd_dVm

        J = real(Cp_index*dSbus_dVm*Cp_index_transpose)

        # Compute update step
        # @time begin
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        
        
        # end
        #precision control
        #dx = round.(dx, digits=6)

        # Update voltage
        if np > 0
            V_gpu[p] .+= dx[j1:j2]
        end
        

        # Evaluate F(x)
        mis = V_gpu .* conj.(Ybus_gpu * V_gpu) - PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, V_gpu, load_gpu, pvarray_gpu)
        F = real(mis[p])

        # Check for convergence
        normF = norm(F, Inf)
        # λ_min, v_min, info = KrylovKit.eigsolve(J, 1, :SR, tol=1e-10)
        # println("Iteration $i: normF = $normF, λ_min = $λ_min")
        if normF < tol
            converged = true
            V = real.(Array(V_gpu))
        end
    end

    return V, converged, i
end
