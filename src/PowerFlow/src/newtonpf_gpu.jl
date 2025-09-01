"""
    newtonpf_gpu(baseMVA, bus, gen, load, Ybus, V0, ref, pv, pq, tol0, max_it0, alg="gpuLU")

Solve AC power flow using Newton's method with GPU acceleration.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `load`: Load data matrix
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference bus index
- `pv`: Vector of PV bus indices
- `pq`: Vector of PQ bus indices
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Algorithm specification for linear solver (default: "gpuLU")

# Returns
- `V`: Final voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed

# Description
This function implements the Newton-Raphson method to solve AC power flow equations
using GPU acceleration. It iteratively updates voltage magnitudes and angles until
the power mismatch falls below the specified tolerance or the maximum number of
iterations is reached.

The algorithm:
1. Initializes voltage values and transfers data to GPU
2. Calculates initial power mismatches
3. Constructs the Jacobian matrix for each iteration
4. Updates voltage values using Newton's method
5. Checks convergence based on power mismatch norm

# Notes
- This implementation uses GPU acceleration for improved performance
- The Jacobian matrix is constructed in block form with partial derivatives
- The function handles both real and reactive power balance constraints
- Data is transferred between CPU and GPU as needed for computation
"""
function newtonpf_gpu(baseMVA,bus,gen,load,pvarray, Ybus, V0, ref, pv, pq, tol0, max_it0, alg="gpuLU")
    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0
    Va = angle.(V)
    Vm = abs.(V)
    nb = length(V)

    # Transfer to GPU
    V_gpu = CuVector(V)
    Va_gpu = CuVector(Va)
    Vm_gpu = CuVector(Vm)
    Ybus_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(Ybus)

    #Transfer matrix to GPU
    bus_gpu = CuArray(bus)
    gen_gpu = CuArray(gen)
    load_gpu = CuArray(load)
    pvarray_gpu = CuArray(pvarray)
    # Set up indexing for updating V
    npv = length(pv)
    npq = length(pq)
    j1 = 1; j2 = npv; # j1:j2 - V angle of pv buses
    j3 = j2 + 1; j4 = j2 + npq; # j3:j4 - V angle of pq buses
    j5 = j4 + 1; j6 = j4 + npq; # j5:j6 - V mag of pq buses

    #Create indexing matrix
    Cpv_pq_index =sparse(1:(npv+npq), vcat(pv, pq), ones(ComplexF64,npv+npq), npv+npq, nb);
    Cpv_pq_index=CuSparseMatrixCSR(Cpv_pq_index)
    Cpv_pq_index_transpose =sparse(vcat(pv, pq), 1:(npv + npq), ones(ComplexF64, npv + npq), nb, npv + npq)
    Cpv_pq_index_transpose=CuSparseMatrixCSR(Cpv_pq_index_transpose)

    Cpq_index =sparse(1:npq, pq, ones(ComplexF64,npq), npq, nb);
    Cpq_index=CuSparseMatrixCSR(Cpq_index)
    Cpq_index_transpose = sparse( pq,1:npq, ones(ComplexF64,npq), nb,npq)
    Cpq_index_transpose=CuSparseMatrixCSR(Cpq_index_transpose)

    Cpv_index = sparse(1:npv, pv, ones(ComplexF64,npv), npv, nb);
    Cpv_index=CuSparseMatrixCSR(Cpv_index)
    # Evaluate F(x0)
    mis = V_gpu .* conj.(Ybus_gpu * V_gpu) - PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu)
    F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]
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
        neg_dSd_dVm = PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu, return_derivative=true)
        dSbus_dVm -= neg_dSd_dVm

        j11 = real(Cpv_pq_index*dSbus_dVa*Cpv_pq_index_transpose);
        j12 = real(Cpv_pq_index*dSbus_dVm*Cpq_index_transpose);
        j21 = imag(Cpq_index*dSbus_dVa*Cpv_pq_index_transpose);
        j22 = imag(Cpq_index*dSbus_dVm*Cpq_index_transpose);

        J = CuSparseMatrixCSR(vcat(hcat(j11, j12), hcat(j21, j22)))

        # Compute update step
        # @time begin
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        
        # end
        #precision control
        #dx = round.(dx, digits=6)

        # Update voltage
        if npv > 0
            Va_gpu[pv] .+= dx[j1:j2]
        end
        if npq > 0
            Va_gpu[pq] .+= dx[j3:j4]
            Vm_gpu[pq] .+= dx[j5:j6]
        end
        V_gpu = Vm_gpu .* exp.(1im * Va_gpu)
        Vm_gpu = abs.(V_gpu) # Update Vm and Va again in case we wrapped around with a negative Vm
        Va_gpu = angle.(V_gpu)

        # Evaluate F(x)
        mis = V_gpu .* conj.(Ybus_gpu * V_gpu) - PowerFlow.makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu)
        F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]

        # Check for convergence
        normF = norm(F, Inf)
        if normF < tol
            converged = true
            V=Array(V_gpu)
        end
    end

    return V, converged, i
end
