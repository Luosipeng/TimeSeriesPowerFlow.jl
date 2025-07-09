"""
    newtondcpf_sp(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")

Solve DC power flow using Newton's method with sparse matrices.

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
- `alg`: Algorithm specification (optional)

# Returns
- `V`: Final voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed

# Description
This function solves DC power flow using Newton's method. In DC systems, 
power P = V * I, where I = G * V (G is the conductance matrix). The function 
iteratively updates voltage magnitudes until the power mismatch falls below 
the specified tolerance or the maximum number of iterations is reached.

The algorithm:
1. Initializes voltage values from the provided starting point
2. Calculates initial power mismatches
3. Constructs the Jacobian matrix for each iteration
4. Updates voltage values using Newton's method
5. Checks convergence based on power mismatch norm

# Notes
- For DC systems, only active power is considered
- The Jacobian matrix represents the partial derivatives of power with respect to voltage
"""
function newtondcpf_sp(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    
    # Initialize
    converged = false
    i = 0
    V = V0
    
    # Set indices
    np = length(p)
    
    # Calculate initial power mismatch
    # In DC systems, power P = V * I, where I = G * V (G is the conductance matrix)
    Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
    # DC systems only consider active power
    P_calc = real.(V .* (Ybus * V))
    P_spec = real.(Sbus)
    
    # Calculate power mismatch
    F = P_spec[p] - P_calc[p]
    
    # Check convergence
    normF = norm(F, Inf)
    if normF < tol
        converged = true
    end
    
    # Newton iteration
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1
        
        # Build Jacobian matrix
        # For DC systems, J = ∂P/∂V
        J = zeros(np, np)
        for j in 1:np
            for k in 1:np
                if j == k
                    # Diagonal elements
                    J[j,j] = 2 * real(Ybus[p[j],p[j]]) * V[p[j]] + 
                             sum(real(Ybus[p[j],m]) * V[m] for m in 1:length(V) if m != p[j])
                else
                    # Off-diagonal elements
                    J[j,k] = real(Ybus[p[j],p[k]]) * V[p[j]]
                end
            end
        end
        
        # Calculate update step
        dV = J \ F
        
        # Update voltage
        V[p] += dV
        
        # Recalculate power mismatch
        Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        P_calc = real.(V .* (Ybus * V))
        P_spec = real.(Sbus)
        F = P_spec[p] - P_calc[p]
        
        # Check convergence
        normF = norm(F, Inf)
        println("Iteration $i: normF = $normF")
        if normF < tol
            converged = true
        end
    end
    
    return V, converged, i
end
