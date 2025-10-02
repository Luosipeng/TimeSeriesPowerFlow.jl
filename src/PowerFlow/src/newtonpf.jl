"""
    newtonpf(baseMVA::Float64, bus::Matrix{Float64}, gen::Matrix{Float64}, 
        load::Matrix{Float64}, pvarray, Ybus::SparseArrays.SparseMatrixCSC{ComplexF64}, V0::Vector{ComplexF64}, 
        ref::Vector{Int}, pv::Vector{Int}, pq::Vector{Int}, 
        tol0::Float64, max_it0::Int, alg::String="bicgstab")

Solve AC power flow using Newton-Raphson method with load and PV array models.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `load`: Load data matrix
- `pvarray`: PV array data
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference bus indices
- `pv`: PV bus indices
- `pq`: PQ bus indices
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Algorithm specification for linear solver (default: "bicgstab")

# Returns
- `V`: Final voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed
- `norm_history`: Array of norm values for each iteration

# Description
This function implements the Newton-Raphson method to solve AC power flow equations
with load and PV array models. It iteratively updates voltage magnitudes and angles until
the power mismatch falls below the specified tolerance or the maximum number of
iterations is reached.
"""
function newtonpf(baseMVA::Float64, bus::Matrix{Float64}, gen::Matrix{Float64}, 
    load::Matrix{Float64}, pvarray, Ybus::SparseArrays.SparseMatrixCSC{ComplexF64}, V0::Vector{ComplexF64}, 
    ref::Vector{Int}, pv::Vector{Int}, pq::Vector{Int}, 
    tol0::Float64, max_it0::Int, alg::String="bicgstab")

    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0
    Va = angle.(V)
    Vm = abs.(V)
    
    # Create array to record norm at each iteration
    norm_history = Float64[]
    
    # Set up indexing for updating V
    npv = length(pv)
    npq = length(pq)
    j1 = 1; j2 = npv; # j1:j2 - V angle of pv buses
    j3 = j2 + 1; j4 = j2 + npq; # j3:j4 - V angle of pq buses
    j5 = j4 + 1; j6 = j4 + npq; # j5:j6 - V mag of pq buses
    
    # Evaluate F(x0)
    mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, Vm, load, pvarray)
    F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]
    
    # Check tolerance
    normF = norm(F, Inf)
    push!(norm_history, normF)  # Record initial norm
    
    if normF < tol
        converged = true
    end
    
    # Do Newton iterations
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1

        # Evaluate Jacobian
        dSbus_dVa, dSbus_dVm = PowerFlow.dSbus_dV(Ybus, V)
        neg_dSd_dVm = PowerFlow.makeSbus(baseMVA, bus, gen, Vm, load, pvarray, return_derivative=true)
        dSbus_dVm .-= neg_dSd_dVm

        j11 = real(dSbus_dVa[vcat(pv, pq), vcat(pv, pq)])
        j12 = real(dSbus_dVm[vcat(pv, pq), pq])
        j21 = imag(dSbus_dVa[pq, vcat(pv, pq)])
        j22 = imag(dSbus_dVm[pq, pq])

        J = [j11 j12; j21 j22]

        # Compute update step
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        
        # Update voltage
        if npv > 0
            Va[pv] .+= dx[j1:j2]
        end
        if npq > 0
            Va[pq] .+= dx[j3:j4]
            Vm[pq] .+= dx[j5:j6]
        end
        V = Vm .* exp.(1im * Va)
        Vm = abs.(V) # Update Vm and Va again in case we wrapped around with a negative Vm
        Va = angle.(V)

        # Evaluate F(x)
        mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, Vm, load, pvarray)
        F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]

        # Check for convergence
        normF = norm(F, Inf)
        push!(norm_history, normF)  # Record current iteration norm
        
        if normF < tol
            converged = true
        end
    end

    return V, converged, i, norm_history
end

"""
    newtonpf(baseMVA::Float64, bus::Matrix{Float64}, gen::Matrix{Float64}, 
             Ybus::SparseArrays.SparseMatrixCSC{ComplexF64}, V0::Vector{ComplexF64}, 
             ref::Vector{Int}, pv::Vector{Int}, pq::Vector{Int}, 
             tol0::Float64, max_it0::Int, alg::String="bicgstab")

Solve AC power flow using Newton-Raphson method.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `Ybus`: Bus admittance matrix
- `V0`: Initial voltage vector
- `ref`: Reference bus indices
- `pv`: PV bus indices
- `pq`: PQ bus indices
- `tol0`: Convergence tolerance
- `max_it0`: Maximum number of iterations
- `alg`: Algorithm specification for linear solver (default: "bicgstab")

# Returns
- `V`: Final voltage vector solution
- `converged`: Boolean indicating whether the algorithm converged
- `i`: Number of iterations performed
- `norm_history`: Array of norm values for each iteration

# Description
This function implements the Newton-Raphson method to solve AC power flow equations.
It iteratively updates voltage magnitudes and angles until the power mismatch falls
below the specified tolerance or the maximum number of iterations is reached.
"""
function newtonpf(baseMVA::Float64, bus::Matrix{Float64}, gen::Matrix{Float64}, 
                  Ybus::SparseArrays.SparseMatrixCSC{ComplexF64}, V0::Vector{ComplexF64}, 
                  ref::Vector{Int}, pv::Vector{Int}, pq::Vector{Int}, 
                  tol0::Float64, max_it0::Int, alg::String="bicgstab")
    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0
    Va = angle.(V)
    Vm = abs.(V)
    
    # Create array to record norm at each iteration
    norm_history = Float64[]
    
    # Set up indexing for updating V
    npv = length(pv)
    npq = length(pq)
    j1 = 1; j2 = npv; # j1:j2 - V angle of pv buses
    j3 = j2 + 1; j4 = j2 + npq; # j3:j4 - V angle of pq buses
    j5 = j4 + 1; j6 = j4 + npq; # j5:j6 - V mag of pq buses
    
    # Evaluate F(x0)
    mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, Vm)
    F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]
    
    # Check tolerance
    normF = norm(F, Inf)
    push!(norm_history, normF)  # Record initial norm
    
    if normF < tol
        converged = true
    end
    
    # Do Newton iterations
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1

        # Evaluate Jacobian
        dSbus_dVa, dSbus_dVm = PowerFlow.dSbus_dV(Ybus, V)
        neg_dSd_dVm = PowerFlow.makeSbus(baseMVA, bus, gen, Vm, return_derivative=true)
        dSbus_dVm .-= neg_dSd_dVm

        j11 = real(dSbus_dVa[vcat(pv, pq), vcat(pv, pq)])
        j12 = real(dSbus_dVm[vcat(pv, pq), pq])
        j21 = imag(dSbus_dVa[pq, vcat(pv, pq)])
        j22 = imag(dSbus_dVm[pq, pq])

        J = [j11 j12; j21 j22]

        # Compute update step
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        
        # Update voltage
        if npv > 0
            Va[pv] .+= dx[j1:j2]
        end
        if npq > 0
            Va[pq] .+= dx[j3:j4]
            Vm[pq] .+= dx[j5:j6]
        end
        V = Vm .* exp.(1im * Va)
        Vm = abs.(V) # Update Vm and Va again in case we wrapped around with a negative Vm
        Va = angle.(V)

        # Evaluate F(x)
        mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, Vm)
        F = [real(mis[vcat(pv, pq)]); imag(mis[pq])]

        # Check for convergence
        normF = norm(F, Inf)
        push!(norm_history, normF)  # Record current iteration norm
        
        if normF < tol
            converged = true
        end
    end

    return V, converged, i, norm_history
end
