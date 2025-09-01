

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
