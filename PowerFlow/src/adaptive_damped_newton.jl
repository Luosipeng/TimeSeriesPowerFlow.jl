"""
    adaptive_damped_newton(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")

Solve power flow using an adaptive damped Newton-Raphson method.

This function implements a Newton-Raphson method with adaptive damping factor
to improve convergence in difficult cases. The damping factor is automatically
adjusted based on the mismatch norm during iterations.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `load`: Load data matrix
- `pvarray`: PV array data
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference bus index
- `p`: Vector of bus indices
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Linear solver algorithm (optional)

# Returns
- `V`: Final voltage solution vector
- `converged`: Boolean indicating convergence status
- `i`: Number of iterations performed
"""
function adaptive_damped_newton(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0
    
    # Set up indexing for updating V
    np = length(p)
    j1 = 1; j2 = np; # j1:j2 - V angle of pv buses
    
    # Evaluate F(x0)
    mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
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
        dSbus_dVa, dSbus_dVm = PowerFlow.dSbus_dV(Ybus, V)
        # neg_dSd_dVm = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray, return_derivative=true)
        # dSbus_dVm .-= neg_dSd_dVm

        J = real(dSbus_dVm[p,p])

        # Compute update step
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        # dx = J \ -F
        
        # Create a copy of voltage update for testing different damping factors
        V_temp = copy(V)
        
        # Initial damping factor
        damping = 1.0
        
        # Apply initial update and calculate new mismatch
        if np > 0
            V_temp[p] .+= damping * dx[j1:j2]
        end
        
        mis_temp = V_temp .* conj.(Ybus * V_temp) - PowerFlow.makeSbus(baseMVA, bus, gen, V_temp, load, pvarray)
        F_temp = real(mis_temp[p])
        new_normF = norm(F_temp, Inf)
        
        # Adaptively adjust damping factor
        # If new mismatch is larger, reduce damping factor
        min_damping = 0.1  # Minimum damping factor limit
        max_attempts = 5   # Maximum number of attempts
        attempt = 0
        
        while new_normF > normF && damping > min_damping && attempt < max_attempts
            # Reduce damping factor
            damping *= 0.5
            attempt += 1
            
            # Recalculate with new damping factor
            V_temp = copy(V)
            if np > 0
                V_temp[p] .+= damping * dx[j1:j2]
            end
            
            mis_temp = V_temp .* conj.(Ybus * V_temp) - PowerFlow.makeSbus(baseMVA, bus, gen, V_temp, load, pvarray)
            F_temp = real(mis_temp[p])
            new_normF = norm(F_temp, Inf)
        end
        
        # Apply final determined damping factor to update voltage
        if np > 0
            V[p] .+= damping * dx[j1:j2]
        end
        
        # Evaluate mismatch at new state
        mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        F = real(mis[p])

        # Check convergence
        normF = norm(F, Inf)
        
        # Optional: Output iteration information, including damping factor
        println("Iteration $i: normF = $normF, damping = $damping")
        λ_min, v_min, info = KrylovKit.eigsolve(J, 1, :SR, tol=1e-10)
        println("Minimum eigenvalue: ", λ_min[1])
        
        if normF < tol
            converged = true
        end
    end

    return V, converged, i
end
