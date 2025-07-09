# Import necessary libraries
using LinearAlgebra
include("types.jl")

function sequential_quadratic_programming(nonlinear::NonConvexOPT, sqp::SQP)
    # Initialize variables
    x = copy(nonlinear.x0)
    n = length(x)
    
    # Evaluate initial functions
    obj = nonlinear.f(x)
    rg = nonlinear.g(x)
    rh = nonlinear.h(x)
    
    # Get problem dimensions
    neq = length(rg)
    niq = length(rh)
    
    # Initialize multipliers more carefully
    λ = zeros(neq)
    μ = max.(zeros(niq), 0.1)  # Start with small positive values
    
    # Initialize Hessian approximation (BFGS) - start with scaled identity
    obj_scale = max(abs(obj), 1.0)
    B = Matrix(obj_scale * I, n, n)
    
    # Evaluate gradients
    df = nonlinear.∇f(x)
    dg = neq > 0 ? nonlinear.∇g(x) : zeros(n, 0)
    dh = niq > 0 ? nonlinear.∇h(x) : zeros(n, 0)
    
    # Initialize convergence flags
    converged = false
    eflag = 0
    iter = 0
    
    # Calculate initial convergence conditions
    maxh = niq > 0 ? maximum(max.(rh, 0.0)) : 0.0
    feascond = maximum([norm(rg, Inf), maxh]) / (1 + norm(x, Inf))
    gradL = df + (neq > 0 ? dg * λ : zeros(n)) + (niq > 0 ? dh * μ : zeros(n))
    gradcond = norm(gradL, Inf) / (1 + norm(x, Inf))
    compcond = niq > 0 ? norm(min.(μ, -rh), Inf) / (1 + norm(x, Inf)) : 0.0
    costcond = 0.0
    
    # Initialize history
    hist = History(zeros(n, sqp.max_iter + 1), zeros(neq, sqp.max_iter + 1), 
                  zeros(niq, sqp.max_iter + 1), zeros(niq, sqp.max_iter + 1),
                  zeros(sqp.max_iter + 1), zeros(sqp.max_iter + 1), 
                  zeros(sqp.max_iter + 1), zeros(sqp.max_iter + 1), zeros(sqp.max_iter + 1))
    
    # Record initial iteration
    hist.x_record[:, 1] = x
    hist.λ_record[:, 1] = λ
    hist.μ_record[:, 1] = μ
    hist.z_record[:, 1] = zeros(niq)  # Not used in SQP but kept for consistency
    hist.obj_record[1] = obj
    hist.feascond_record[1] = feascond
    hist.gradcond_record[1] = gradcond
    hist.compcond_record[1] = compcond
    hist.costcond_record[1] = costcond
    
    # Check initial convergence
    if feascond < sqp.feasible_tol && gradcond < sqp.feasible_tol && compcond < sqp.feasible_tol
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
        return SQP_Solution(x, λ, μ, obj, eflag, 0, hist)
    end
    
    obj0 = obj
    
    # Main SQP iteration
    while iter < sqp.max_iter && !converged
        iter += 1
        
        # Solve QP subproblem
        p = nothing
        λ_new = nothing
        μ_new = nothing
        
        try
            p, λ_new, μ_new = solve_qp_subproblem(B, df, dg, dh, rg, rh, n, neq, niq)
        catch e
            println("Warning: QP subproblem failed at iteration $iter: $e")
            eflag = -1
            break
        end
        
        # Check for numerical issues
        if p === nothing || any(isnan.(p))
            println("Warning: Invalid search direction at iteration $iter")
            eflag = -1
            break
        end
        
        # Line search with merit function
        α = line_search_merit(nonlinear, x, p, λ, μ, rg, rh, sqp.penalty_param)
        
        if α < sqp.min_step_size
            println("Warning: Step size too small at iteration $iter")
            eflag = -2
            break
        end
        
        # Update variables
        x_new = x + α * p
        obj_new = nonlinear.f(x_new)
        rg_new = nonlinear.g(x_new)
        rh_new = nonlinear.h(x_new)
        
        # Evaluate new gradients
        df_new = nonlinear.∇f(x_new)
        dg_new = neq > 0 ? nonlinear.∇g(x_new) : zeros(n, 0)
        dh_new = niq > 0 ? nonlinear.∇h(x_new) : zeros(n, 0)
        
        # Update Hessian approximation using BFGS
        s = x_new - x
        gradL_new = df_new + (neq > 0 ? dg_new * λ_new : zeros(n)) + (niq > 0 ? dh_new * μ_new : zeros(n))
        y = gradL_new - gradL
        
        B = bfgs_update(B, s, y)
        
        # Update variables for next iteration
        x = x_new
        obj = obj_new
        rg = rg_new
        rh = rh_new
        df = df_new
        dg = dg_new
        dh = dh_new
        λ = λ_new
        μ = μ_new
        gradL = gradL_new
        
        # Check convergence
        maxh = niq > 0 ? maximum(max.(rh, 0.0)) : 0.0
        feascond = maximum([norm(rg, Inf), maxh]) / (1 + norm(x, Inf))
        gradcond = norm(gradL, Inf) / (1 + norm(x, Inf))
        compcond = niq > 0 ? norm(min.(μ, -rh), Inf) / (1 + norm(x, Inf)) : 0.0
        costcond = abs(obj - obj0) / (1 + abs(obj0))
        
        # Record iteration results
        hist.x_record[:, iter + 1] = x
        hist.λ_record[:, iter + 1] = λ
        hist.μ_record[:, iter + 1] = μ
        hist.z_record[:, iter + 1] = zeros(niq)
        hist.obj_record[iter + 1] = obj
        hist.feascond_record[iter + 1] = feascond
        hist.gradcond_record[iter + 1] = gradcond
        hist.compcond_record[iter + 1] = compcond
        hist.costcond_record[iter + 1] = costcond
        
        obj0 = obj
        
        # Check convergence
        if feascond < sqp.feasible_tol && gradcond < sqp.feasible_tol && compcond < sqp.feasible_tol
            converged = true
            eflag = 1
            break
        end
    end
    
    # Set final status
    if !converged && eflag == 0
        eflag = 0  # Maximum iterations reached
    end
    
    # Trim history arrays
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
    
    return SQP_Solution(x, λ, μ, obj, eflag, iter, hist)
end

function solve_qp_subproblem(B, df, dg, dh, rg, rh, n, neq, niq)
    # Solve: min 0.5*p'*B*p + df'*p
    # s.t.: dg'*p + rg = 0  (equality constraints)
    #       dh'*p + rh ≤ 0  (inequality constraints)
    
    # Use active set method for QP
    active_ineq = findall(x -> x > -1e-6, rh)  # Near-active inequalities
    
    # Form KKT system with active constraints
    A = [dg'; dh'[active_ineq, :]]
    b = [-rg; -rh[active_ineq]]
    
    nactive = length(active_ineq)
    ncon = neq + nactive
    
    if ncon > 0
        # Solve KKT system: [B A'; A 0] [p; λμ] = [-df; b]
        KKT = [B A'; A zeros(ncon, ncon)]
        rhs = [-df; b]
        
        try
            sol = KKT \ rhs
            p = sol[1:n]
            λμ = sol[n+1:end]
            
            λ_new = neq > 0 ? λμ[1:neq] : Float64[]
            μ_active = nactive > 0 ? λμ[neq+1:end] : Float64[]
            
            # Reconstruct full μ vector
            μ_new = zeros(niq)
            if nactive > 0
                μ_new[active_ineq] = μ_active
            end
            
        catch
            # Fallback: regularized solve
            reg = 1e-8 * I(size(KKT, 1))
            sol = (KKT + reg) \ rhs
            p = sol[1:n]
            λμ = sol[n+1:end]
            
            λ_new = neq > 0 ? λμ[1:neq] : Float64[]
            μ_active = nactive > 0 ? λμ[neq+1:end] : Float64[]
            μ_new = zeros(niq)
            if nactive > 0
                μ_new[active_ineq] = μ_active
            end
        end
    else
        # Unconstrained QP
        p = -B \ df
        λ_new = Float64[]
        μ_new = zeros(niq)
    end
    
    return p, λ_new, μ_new
end

function line_search_merit(nonlinear, x, p, λ, μ, rg, rh, penalty_param)
    # Merit function: f(x) + penalty_param * (||g(x)||₁ + ||max(0,h(x))||₁)
    
    α = 1.0
    c1 = 1e-4  # Armijo parameter
    ρ = 0.5    # Backtracking parameter
    
    # Current merit function value
    merit_current = nonlinear.f(x) + penalty_param * (norm(rg, 1) + norm(max.(rh, 0.0), 1))
    
    # Directional derivative (more robust calculation)
    df = nonlinear.∇f(x)
    dg = length(rg) > 0 ? nonlinear.∇g(x) : zeros(length(x), 0)
    dh = length(rh) > 0 ? nonlinear.∇h(x) : zeros(length(x), 0)
    
    # More accurate directional derivative
    dmerit = df' * p
    if length(rg) > 0
        dmerit += penalty_param * sum(sign.(rg) .* (dg' * p))
    end
    if length(rh) > 0
        active_h = rh .> 0
        if any(active_h)
            dmerit += penalty_param * sum((dh' * p)[active_h])
        end
    end
    
    max_backtracks = 30  # Increased from 20
    backtrack_count = 0
    min_alpha = 1e-10    # More aggressive minimum
    
    while α > min_alpha && backtrack_count < max_backtracks
        x_trial = x + α * p
        
        try
            obj_trial = nonlinear.f(x_trial)
            rg_trial = nonlinear.g(x_trial)
            rh_trial = nonlinear.h(x_trial)
            
            merit_trial = obj_trial + penalty_param * (norm(rg_trial, 1) + norm(max.(rh_trial, 0.0), 1))
            
            # Armijo condition with better numerical handling
            armijo_rhs = merit_current + c1 * α * dmerit
            if merit_trial <= armijo_rhs || α < 1e-8
                return max(α, 1e-12)  # Ensure minimum step size
            end
        catch
            # Function evaluation failed, reduce step size more aggressively
            α *= 0.1
            backtrack_count += 1
            continue
        end
        
        α *= ρ
        backtrack_count += 1
    end
    
    return max(α, 1e-12)  # Return minimum acceptable step size
end

function bfgs_update(B, s, y)
    # BFGS update: B_new = B - (B*s*s'*B)/(s'*B*s) + (y*y')/(y'*s)
    
    if abs(s' * y) < 1e-12
        return B  # Skip update if curvature condition fails
    end
    
    Bs = B * s
    sBs = s' * Bs
    sy = s' * y
    
    if sBs < 1e-12 || sy < 1e-12
        return B  # Skip update for numerical stability
    end
    
    B_new = B - (Bs * Bs') / sBs + (y * y') / sy
    
    # Ensure positive definiteness (simple check)
    if any(diag(B_new) .< 1e-12)
        return B  # Revert to previous Hessian if update fails
    end
    
    return B_new
end
