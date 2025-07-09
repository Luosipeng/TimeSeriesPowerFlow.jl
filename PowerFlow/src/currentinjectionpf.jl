"""
    currentinjectionpf(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")

Solve power flow using the Current Injection method, particularly suitable for resistive networks.

# Arguments
- `baseMVA`: Base MVA for power system normalization
- `bus`: Matrix containing bus data
- `gen`: Matrix containing generator data
- `load`: Matrix containing load data
- `pvarray`: Array of PV bus information
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference (slack) bus index
- `p`: Array of PQ bus indices
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Algorithm variant (optional)

# Returns
- `V`: Final complex voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed

# Description
This function implements the Current Injection power flow method, which is particularly
effective for networks with high R/X ratios (resistive networks). The method works by:
1. Converting power injections to current injections
2. Calculating current mismatches
3. Solving for voltage updates using the real part of the admittance matrix
4. Iterating until convergence or maximum iterations reached
"""
function currentinjectionpf(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    
    # Initialize
    converged = false
    i = 0
    V = V0
    
    # Set up indexing for updating V
    np = length(p)
    
    # Do Current Injection iterations
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1
        
        # Calculate power injections
        Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        
        # Calculate current injections: I = S*/V*
        Ibus = conj.(Sbus ./ V)
        
        # Calculate current mismatches: ΔI = Ybus*V - Ibus
        Icalc = Ybus * V
        dI = Icalc - Ibus
        
        # Extract real part of current mismatches for PQ buses
        F = real(dI[p])
        
        # Check for convergence
        normF = norm(F, Inf)
        if normF < tol
            converged = true
            break
        end
        
        # For pure resistive networks, Ybus is real
        # Form the Jacobian matrix (for current equations)
        # J = real part of Ybus for PQ buses
        J = real(Ybus[p, p])
        
        # Solve for voltage updates
        dV = J \ -F
        
        # Update voltages
        if np > 0
            V[p] .+= dV
        end
        println("Iteration $i: normF = $normF")
    end
    
    return V, converged, i
end
