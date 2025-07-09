# Import necessary libraries
using LinearAlgebra
include("types.jl")
# Define the objective function

function interior_point_method(nonlinear::NonConvexOPT, ipm::IPM)
    # define the initial point, we have four types of initial points: x, z, λ, μ
    x = nonlinear.x0
    # evaluate the cost, constraints, derivaties
    obj = nonlinear.f(x)
    rg = nonlinear.g(x)
    rh = nonlinear.h(x)
    # get the size of the problem
    n = length(x) # number of variables
    niq = length(rh) # number of inequality constraints
    neq = length(rg) # number of equality constraints
    # initialize the multipliers
    z = ones(niq)
    λ = zeros(neq)
    μ = ones(niq)
    z0 = 1.0
    eflag = 0
    γ = 1.0
    σ = 0.1
    # Algorithm constants
    ξ = ipm.ξ
    alpha_min = 1e-8
    max_stepsize = 1e10
    
    # update z and μ with better initialization
    for i in 1:niq
        if rh[i] < -z0
            z[i] = max(-rh[i], 1e-8)  # Ensure positive z
        else
            z[i] = z0
        end
    end
    # Initialize mu more carefully
    for i in 1:niq
        if γ / z[i] > z0
            μ[i] = γ / z[i]
        else
            μ[i] = γ / z0
        end
    end
    
    e = ones(niq)
    obj0 = obj
    
    # Compute initial derivatives
    dg = nonlinear.∇g(x)
    dh = nonlinear.∇h(x)
    df = nonlinear.∇f(x)
    Lx = nonlinear.Lx(x, λ, μ)
    Lxx = nonlinear.Lxx(x, λ, μ)
    
    # initializate
    converged = false
    maxh = niq > 0 ? maximum(rh) : 0.0
    feascond = maximum([norm(rg, Inf), maxh]) / (1 + maximum([norm(x, Inf), norm(z, Inf)]))
    gradcond = norm(Lx, Inf) / (1 + maximum([norm(μ, Inf),norm(λ, Inf)]))
    compcond = niq > 0 ? (z'*μ) / (1 + norm(x, Inf)) : 0.0
    costcond = abs(obj - obj0) / (1 + abs(obj0))
    
    # record the results - initialize with max_iter + 1 to include iteration 0
    hist = History(zeros(n, ipm.max_iter + 1), zeros(neq, ipm.max_iter + 1), zeros(niq, ipm.max_iter + 1), zeros(niq, ipm.max_iter + 1), zeros(ipm.max_iter + 1), zeros(ipm.max_iter + 1), zeros(ipm.max_iter + 1), zeros(ipm.max_iter + 1), zeros(ipm.max_iter + 1))
    
    # Record initial iteration (iteration 0)
    hist.x_record[:, 1] = x
    hist.λ_record[:, 1] = λ
    hist.μ_record[:, 1] = μ
    hist.z_record[:, 1] = z
    hist.obj_record[1] = obj
    hist.feascond_record[1] = feascond
    hist.gradcond_record[1] = gradcond
    hist.compcond_record[1] = compcond
    hist.costcond_record[1] = costcond
    
    # Check initial convergence
    if feascond < ipm.feasible_tol && gradcond < ipm.feasible_tol && compcond < ipm.feasible_tol && costcond < ipm.feasible_tol
        converged = true
        eflag = 1
        hist.x_record = hist.x_record[:, 1:1]
        hist.λ_record = hist.λ_record[:, 1:1]
        hist.μ_record = hist.μ_record[:, 1:1]
        hist.z_record = hist.z_record[:, 1:1]
        hist.obj_record = hist.obj_record[1:1]
        hist.feascond_record = hist.feascond_record[1:1]
        hist.gradcond_record = hist.gradcond_record[1:1]
        hist.compcond_record = hist.compcond_record[1:1]
        hist.costcond_record = hist.costcond_record[1:1]
        return IPM_Solution(x, λ, μ, nonlinear.f(x), eflag, hist)
    end
    
    # start the iteration using the while funciton
    iter = 1
    
    while iter <= ipm.max_iter && !converged
        # calculate the KKT step with improved error handling
        try
            Δx, Δλ, Δμ, Δz = KKT_step(μ, z, e, rg, rh, dg, dh, Lx, Lxx, γ, n, neq)
        catch e
            println("Warning: KKT system solve failed at iteration $iter: $e")
            # Try to continue with a smaller step or exit gracefully
            eflag = -1
            break
        end
        
        # Check for numerical issues with more tolerance
        if any(isnan.(Δx)) || any(isnan.(Δλ)) || any(isnan.(Δμ)) || any(isnan.(Δz))
            println("Warning: NaN detected in step at iteration $iter")
            eflag = -1
            break
        end
        
        # Check for excessively large steps
        step_norm = norm([Δx; Δλ; Δμ; Δz])
        if step_norm > max_stepsize
            println("Warning: Step size too large ($step_norm) at iteration $iter")
            # Scale down the step instead of failing
            scale_factor = max_stepsize / step_norm
            Δx *= scale_factor
            Δλ *= scale_factor  
            Δμ *= scale_factor
            Δz *= scale_factor
        end
        
        # calculate the step length with improved fraction-to-boundary rule
        k = findall(<(0), Δz)
        if length(k) > 0
            α_primal = min(ξ * minimum(-z[k]./Δz[k]), 1.0)
        else
            α_primal = 1.0
        end
        k = findall(<(0), Δμ)
        if length(k) > 0
            α_dual = min(ξ * minimum(-μ[k]./Δμ[k]), 1.0)
        else
            α_dual = 1.0
        end
        
        # Ensure minimum step size with more tolerance
        α_primal = max(α_primal, alpha_min)
        α_dual = max(α_dual, alpha_min)
        
        # update the point
        x = x + α_primal*Δx
        z = z + α_primal*Δz
        λ = λ + α_dual*Δλ
        μ = μ + α_dual*Δμ
        
        # Ensure z and μ remain positive
        z = max.(z, 1e-12)
        μ = max.(μ, 1e-12)

        # Update barrier parameter more carefully
        if niq > 0
            γ = σ * z'*μ / niq
            # Bound gamma to reasonable values
            γ = max(min(γ, 1e6), 1e-12)
        end

        # update the cost, constraints, derivaties   
        obj = nonlinear.f(x) 
        rg = nonlinear.g(x)
        rh = nonlinear.h(x)
        dg = nonlinear.∇g(x)
        dh = nonlinear.∇h(x)
        df = nonlinear.∇f(x)
        Lx = nonlinear.Lx(x, λ, μ)
        Lxx = nonlinear.Lxx(x, λ, μ)

        # Check for NaN in solution
        if any(isnan.(x))
            println("Warning: NaN detected in solution at iteration $iter")
            eflag = -1
            break
        end

        maxh = niq > 0 ? maximum(rh) : 0.0
        feascond = maximum([norm(rg, Inf), maxh]) / (1 + maximum([norm(x, Inf), norm(z, Inf)]))
        gradcond = norm(Lx, Inf) / (1 + maximum([norm(μ, Inf),norm(λ, Inf)]))
        compcond = niq > 0 ? (z'*μ) / (1 + norm(x, Inf)) : 0.0
        costcond = abs(obj - obj0) / (1 + abs(obj0))
        
        # record the results (iter+1 because we start from iteration 0)
        hist.x_record[:, iter + 1] = x
        hist.λ_record[:, iter + 1] = λ
        hist.μ_record[:, iter + 1] = μ
        hist.z_record[:, iter + 1] = z
        hist.obj_record[iter + 1] = obj
        hist.feascond_record[iter + 1] = feascond
        hist.gradcond_record[iter + 1] = gradcond
        hist.compcond_record[iter + 1] = compcond
        hist.costcond_record[iter + 1] = costcond
        obj0 = obj

        # Enhanced convergence check
        if feascond < ipm.feasible_tol && gradcond < ipm.feasible_tol && compcond < ipm.feasible_tol && costcond < ipm.feasible_tol
            eflag = 1
            converged = true
            break
        end
        iter += 1
    end
    
    # Final status check - set eflag to 0 if not converged and no error occurred
    if !converged && eflag == 0
        eflag = 0  # Did not converge within max iterations
    end
    
    # reduce the size of the record (iter+1 because we include iteration 0)
    actual_iters = min(iter + 1, size(hist.x_record, 2))
    hist.x_record = hist.x_record[:, 1:actual_iters]
    hist.λ_record = hist.λ_record[:, 1:actual_iters]
    hist.μ_record = hist.μ_record[:, 1:actual_iters]
    hist.z_record = hist.z_record[:, 1:actual_iters]
    hist.obj_record = hist.obj_record[1:actual_iters]
    hist.feascond_record = hist.feascond_record[1:actual_iters]
    hist.gradcond_record = hist.gradcond_record[1:actual_iters]
    hist.compcond_record = hist.compcond_record[1:actual_iters]
    hist.costcond_record = hist.costcond_record[1:actual_iters]
    
    # Return with iterations count
    return IPM_Solution(x, λ, μ, nonlinear.f(x), eflag, iter, hist)

end


# define the KKT iteration for the 
function KKT_step(μ, z, e, rg, rh, dg, dh, Lx, Lxx, γ, n, neq)
    # Ensure z is positive to avoid division by zero
    z_safe = max.(z, 1e-12)
    zinvdiag = Diagonal(1.0 ./ z_safe)
    mudiag = Diagonal(μ)
    dh_zinv = dh * zinvdiag
    
    # 3.39
    M = Lxx + dh_zinv * mudiag * dh'
    # 3.42 
    N = Lx + dh_zinv * (mudiag*rh + γ*e)
    # 3.43 - Build KKT system more carefully
    KKT = [M dg; dg' zeros(neq, neq)]
    rhs = [-N; -rg]
    
    # Initialize Δ to handle scoping
    Δ = nothing
    
    # Check for conditioning issues and solve
    try
        # Try direct solve first
        Δ = KKT \ rhs
    catch
        # If direct solve fails, try with regularization
        reg = 1e-8 * max(norm(M, Inf), 1.0) * I(size(KKT, 1))
        try
            Δ = (KKT + reg) \ rhs
        catch
            # If that fails too, use a larger regularization
            reg = 1e-6 * max(norm(M, Inf), 1.0) * I(size(KKT, 1))
            Δ = (KKT + reg) \ rhs
        end
    end
    
    # Check if solve was successful
    if Δ === nothing
        error("KKT system solve failed completely")
    end
    
    Δx = Δ[1:n]
    Δλ = Δ[n+1:n+neq]
    # 3.36
    Δz = -rh - z - dh' * Δx
    # 3.35
    Δμ = -μ + zinvdiag*(γ.*e - mudiag * Δz)
    return Δx, Δλ, Δμ, Δz
end