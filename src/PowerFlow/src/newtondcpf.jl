"""
    newtondcpf(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")

Solve DC power flow using Newton's method.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `load`: Load data matrix
- `pvarray`: PV array data
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference bus index
- `p`: Vector of indices for buses to be included in the calculation
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Algorithm specification for linear solver (optional)

# Returns
- `V`: Final voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed

# Description
This function implements the Newton-Raphson method to solve DC power flow equations.
It iteratively updates voltage values until the power mismatch falls below the specified
tolerance or the maximum number of iterations is reached.

The algorithm:
1. Initializes voltage values from the provided starting point
2. Calculates initial power mismatches
3. Constructs the Jacobian matrix for each iteration
4. Updates voltage values using Newton's method
5. Checks convergence based on power mismatch norm

# Notes
- The Jacobian matrix represents the partial derivatives of power with respect to voltage
- This implementation handles both real power balance constraints
"""
function newtondcpf(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
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
        neg_dSd_dVm = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray, return_derivative=true)
        dSbus_dVm -= neg_dSd_dVm

        J = real(dSbus_dVm[p,p])

        # Compute update step
        # @time begin
        # dx, info = PowerFlow.julinsolve(J, -F, alg)
        dx = J \ -F
        
        # end
        #precision control
        #dx = round.(dx, digits=6)

        # Update voltage
        if np > 0
            V[p] .+= dx[j1:j2]
        end
        

        # Evaluate F(x)
        mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        F = real(mis[p])

        # Check for convergence
        normF = norm(F, Inf)
        # λ_min, v_min, info = KrylovKit.eigsolve(J, 1, :SR, tol=1e-10)
        # println("Iteration $i: normF = $normF, λ_min = $λ_min")
        if normF < tol
            converged = true
        end
    end

    return V, converged, i
end
