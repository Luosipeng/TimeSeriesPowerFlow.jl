using LinearAlgebra
using SparseArrays
using Printf

"""
    mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)

MATPOWER Interior Point Solver (MIPS) - Julia implementation.

Primal-dual interior point method for NLP (nonlinear programming).
Minimize a function f(x) beginning from a starting point x0, subject to 
optional linear and nonlinear constraints and variable bounds.

# Arguments
- `f_fcn`: Function handle that evaluates objective function f(x), its gradients and Hessian
- `x0`: Starting value of optimization vector x
- `A`, `l`, `u`: Matrix A and vectors l, u to define linear constraints
- `xmin`, `xmax`: Lower and upper bounds on optimization variables
- `gh_fcn`: Function handle for nonlinear constraints g(x) and h(x)
- `hess_fcn`: Function handle for Hessian of Lagrangian
- `opt`: Options structure

# Returns
- `x`: Solution vector
- `f`: Final objective function value
- `exitflag`: Exit flag (1=converged, 0=max iter, -1=failed)
- `output`: Output structure with iteration info
- `lambda`: Structure with Lagrange multipliers
"""
function mips(f_fcn, x0::Vector{Float64}; 
              A::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int}}=spzeros(0,length(x0)),
              l::Vector{Float64}=Float64[],
              u::Vector{Float64}=Float64[], 
              xmin::Vector{Float64}=fill(-Inf, length(x0)),
              xmax::Vector{Float64}=fill(Inf, length(x0)),
              gh_fcn=nothing,
              hess_fcn=nothing,
              opt=Dict{Symbol,Any}())

    # Input argument handling
    nx = length(x0)
    
    # Set default argument values if missing
    if isempty(A) || (isempty(l) || all(l .== -Inf)) && (isempty(u) || all(u .== Inf))
        A = spzeros(0, nx)
    end
    nA = size(A, 1)
    
    if isempty(u)
        u = fill(Inf, nA)
    end
    if isempty(l)
        l = fill(-Inf, nA)
    end
    
    nonlinear = !isnothing(gh_fcn)
    
    # Default options
    verbose = get(opt, :verbose, 0)
    linsolver = get(opt, :linsolver, "")
    feastol = get(opt, :feastol, 1e-6)
    gradtol = get(opt, :gradtol, 1e-6)
    comptol = get(opt, :comptol, 1e-6)
    costtol = get(opt, :costtol, 1e-6)
    max_it = get(opt, :max_it, 150)
    step_control = get(opt, :step_control, false)
    red_it = get(opt, :red_it, 20)
    cost_mult = get(opt, :cost_mult, 1.0)
    
    # Algorithm constants
    xi = get(opt, :xi, 0.99995)
    sigma = get(opt, :sigma, 0.1)
    z0 = get(opt, :z0, 1.0)
    alpha_min = get(opt, :alpha_min, 1e-8)
    rho_min = get(opt, :rho_min, 0.95)
    rho_max = get(opt, :rho_max, 1.05)
    mu_threshold = get(opt, :mu_threshold, 1e-5)
    max_stepsize = get(opt, :max_stepsize, 1e10)
    full_hist = get(opt, :full_hist, false)
    
    # Validation
    if xi >= 1 || xi < 0.5
        error("mips: opt.xi ($xi) must be a number slightly less than 1")
    end
    if sigma > 1 || sigma <= 0
        error("mips: opt.sigma ($sigma) must be a number between 0 and 1")
    end
    
    # Initialize
    i = 0
    converged = false
    eflag = 0
    
    # Initialize history
    hist = Vector{Dict{Symbol,Any}}()
    
    # Add var limits to linear constraints
    AA = [I(nx); A]
    ll = [xmin; l]
    uu = [xmax; u]
    
    # Split up linear constraints
    ieq = findall(abs.(uu - ll) .<= eps())
    igt = findall((uu .>= 1e10) .& (ll .> -1e10))
    ilt = findall((ll .<= -1e10) .& (uu .< 1e10))
    ibx = findall((abs.(uu - ll) .> eps()) .& (uu .< 1e10) .& (ll .> -1e10))
    
    Ae = AA[ieq, :]
    be = uu[ieq]
    Ai = vcat(AA[ilt, :], -AA[igt, :], AA[ibx, :], -AA[ibx, :])
    bi = vcat(uu[ilt], -ll[igt], uu[ibx], -ll[ibx])
    
    # Evaluate cost f(x0) and constraints g(x0), h(x0)
    x = copy(x0)
    f, df = f_fcn(x)
    f *= cost_mult
    df *= cost_mult
    
    if nonlinear
        hn, gn, dhn, dgn = gh_fcn(x)
        h = vcat(hn, Ai * x - bi)
        g = vcat(gn, Ae * x - be)
        dh = hcat(dhn, Ai')
        dg = hcat(dgn, Ae')
    else
        h = Ai * x - bi
        g = Ae * x - be
        dh = Ai
        dg = Ae
    end
    
    # Grab dimensions
    neq = length(g)
    niq = length(h)
    neqnln = nonlinear ? length(gn) : 0
    niqnln = nonlinear ? length(hn) : 0
    nlt = length(ilt)
    ngt = length(igt)
    nbx = length(ibx)
    
    # Initialize gamma, lam, mu, z, e
    gamma = 1.0
    lam = zeros(neq)
    z = fill(z0, niq)
    mu = copy(z)
    e = ones(niq)
    
    k = findall(h .< -z0)
    z[k] = -h[k]
    k = findall(gamma ./ z .> z0)
    if !isempty(k)
        mu[k] = gamma ./ z[k]
    end
    
    # Check tolerance
    f0 = f
    if step_control
        L = f + dot(lam, g) + dot(mu, h + z) - gamma * sum(log.(z))
    end
    Lx = df + dg * lam + dh * mu
    
    maxh = isempty(h) ? 0.0 : maximum(h)
    feascond = max(norm(g, Inf), maxh) / (1 + max(norm(x, Inf), norm(z, Inf)))
    gradcond = norm(Lx, Inf) / (1 + max(norm(lam, Inf), norm(mu, Inf)))
    compcond = dot(z, mu) / (1 + norm(x, Inf))
    costcond = abs(f - f0) / (1 + abs(f0))
    
    # Save initial history
    hist_entry = Dict(:feascond => feascond, :gradcond => gradcond,
                     :compcond => compcond, :costcond => costcond,
                     :gamma => gamma, :stepsize => 0.0, :obj => f/cost_mult,
                     :alphap => 0.0, :alphad => 0.0)
    if full_hist
        hist_entry[:x] = copy(x)
        hist_entry[:z] = copy(z)
        hist_entry[:g] = copy(g)
        hist_entry[:h] = copy(h)
        hist_entry[:lam] = copy(lam)
        hist_entry[:mu] = copy(mu)
    end
    push!(hist, hist_entry)
    
    if verbose > 0
        s = step_control ? "-sc" : ""
        println("MATPOWER Interior Point Solver -- MIPS$s (Julia)")
        if verbose > 1
            println(" it    objective   step size   feascond     gradcond     compcond     costcond")
            println("----  ------------ --------- ------------ ------------ ------------ ------------")
            @printf("%3d  %12.8g %10s %12g %12g %12g %12g\n",
                   i, f/cost_mult, "", feascond, gradcond, compcond, costcond)
        end
    end
    
    if feascond < feastol && gradcond < gradtol && compcond < comptol && costcond < costtol
        converged = true
        if verbose > 0
            println("Converged!")
        end
    end
    
    # Main Newton iterations
    while !converged && i < max_it
        i += 1
        
        # Compute update step
        lambda_struct = Dict(:eqnonlin => lam[1:neqnln], :ineqnonlin => mu[1:niqnln])
        
        if nonlinear
            if isnothing(hess_fcn)
                error("mips: Hessian evaluation via finite differences not yet implemented.\n       Please provide your own hessian evaluation function.")
            end
            Lxx = hess_fcn(x, lambda_struct, cost_mult)
        else
            _, _, d2f = f_fcn(x)
            Lxx = d2f * cost_mult
        end
        
        zinvdiag = spdiagm(1.0 ./ z)
        mudiag = spdiagm(mu)
        dh_zinv = dh * zinvdiag
        M = Lxx + dh_zinv * mudiag * dh'
        N = Lx + dh_zinv * (mudiag * h + gamma * e)
        
        KKT = [M dg; dg' spzeros(neq, neq)]
        rhs = [-N; -g]
        dxdlam = KKT \ rhs
        
        if any(isnan.(dxdlam)) || norm(dxdlam) > max_stepsize
            if verbose > 0
                println("Numerically Failed")
            end
            eflag = -1
            break
        end
        
        dx = dxdlam[1:nx]
        dlam = dxdlam[nx+1:nx+neq]
        dz = -h - z - dh' * dx
        dmu = -mu + zinvdiag * (gamma * e - mudiag * dz)
        
        # Optional step-size control
        sc = false
        if step_control
            x1 = x + dx
            
            # Evaluate cost, constraints, derivatives at x1
            f1, df1 = f_fcn(x1)
            f1 *= cost_mult
            df1 *= cost_mult
            
            if nonlinear
                hn1, gn1, dhn1, dgn1 = gh_fcn(x1)
                h1 = vcat(hn1, Ai * x1 - bi)
                g1 = vcat(gn1, Ae * x1 - be)
                dh1 = hcat(dhn1, Ai')
                dg1 = hcat(dgn1, Ae')
            else
                h1 = Ai * x1 - bi
                g1 = Ae * x1 - be
                dh1 = dh
                dg1 = dg
            end
            
            # Check tolerance
            Lx1 = df1 + dg1 * lam + dh1 * mu
            maxh1 = isempty(h1) ? 0.0 : maximum(h1)
            feascond1 = max(norm(g1, Inf), maxh1) / (1 + max(norm(x1, Inf), norm(z, Inf)))
            gradcond1 = norm(Lx1, Inf) / (1 + max(norm(lam, Inf), norm(mu, Inf)))
            
            if feascond1 > feascond && gradcond1 > gradcond
                sc = true
            end
        end
        
        if sc
            alpha = 1.0
            for j in 1:red_it
                dx1 = alpha * dx
                x1 = x + dx1
                f1, _ = f_fcn(x1)
                f1 *= cost_mult
                
                if nonlinear
                    hn1, gn1 = gh_fcn(x1)
                    h1 = vcat(hn1, Ai * x1 - bi)
                    g1 = vcat(gn1, Ae * x1 - be)
                else
                    h1 = Ai * x1 - bi
                    g1 = Ae * x1 - be
                end
                
                L1 = f1 + dot(lam, g1) + dot(mu, h1 + z) - gamma * sum(log.(z))
                
                if verbose > 2
                    @printf("\n   %3d            %10g", -j, norm(dx1))
                end
                
                rho = (L1 - L) / (dot(Lx, dx1) + 0.5 * dot(dx1, Lxx * dx1))
                
                if rho > rho_min && rho < rho_max
                    break
                else
                    alpha = alpha / 2
                end
            end
            dx = alpha * dx
            dz = alpha * dz
            dlam = alpha * dlam
            dmu = alpha * dmu
        end
        
        # Update step sizes
        k = findall(dz .< 0)
        alphap = isempty(k) ? 1.0 : min(xi * minimum(z[k] ./ (-dz[k])), 1.0)
        
        k = findall(dmu .< 0)
        alphad = isempty(k) ? 1.0 : min(xi * minimum(mu[k] ./ (-dmu[k])), 1.0)
        
        # Do the update
        x += alphap * dx
        z += alphap * dz
        lam += alphad * dlam
        mu += alphad * dmu
        
        if niq > 0
            gamma = sigma * dot(z, mu) / niq
        end
        
        # Evaluate cost, constraints, derivatives
        f, df = f_fcn(x)
        f *= cost_mult
        df *= cost_mult
        
        if nonlinear
            hn, gn, dhn, dgn = gh_fcn(x)
            h = vcat(hn, Ai * x - bi)
            g = vcat(gn, Ae * x - be)
            dh = hcat(dhn, Ai')
            dg = hcat(dgn, Ae')
        else
            h = Ai * x - bi
            g = Ae * x - be
        end
        
        # Check tolerance
        Lx = df + dg * lam + dh * mu
        maxh = isempty(h) ? 0.0 : maximum(h)
        feascond = max(norm(g, Inf), maxh) / (1 + max(norm(x, Inf), norm(z, Inf)))
        gradcond = norm(Lx, Inf) / (1 + max(norm(lam, Inf), norm(mu, Inf)))
        compcond = dot(z, mu) / (1 + norm(x, Inf))
        costcond = abs(f - f0) / (1 + abs(f0))
        
        # Save history
        hist_entry = Dict(:feascond => feascond, :gradcond => gradcond,
                         :compcond => compcond, :costcond => costcond,
                         :gamma => gamma, :stepsize => norm(dx), :obj => f/cost_mult,
                         :alphap => alphap, :alphad => alphad)
        if full_hist
            hist_entry[:x] = copy(x)
            hist_entry[:z] = copy(z)
            hist_entry[:g] = copy(g)
            hist_entry[:h] = copy(h)
            hist_entry[:lam] = copy(lam)
            hist_entry[:mu] = copy(mu)
        end
        push!(hist, hist_entry)
        
        if verbose > 1
            @printf("%3d  %12.8g %10.5g %12g %12g %12g %12g\n",
                   i, f/cost_mult, norm(dx), feascond, gradcond, compcond, costcond)
        end
        
        if feascond < feastol && gradcond < gradtol && compcond < comptol && costcond < costtol
            converged = true
            if verbose > 0
                println("Converged!")
            end
        else
            if any(isnan.(x)) || alphap < alpha_min || alphad < alpha_min || gamma < eps() || gamma > 1/eps()
                if verbose > 0
                    println("Numerically Failed")
                end
                eflag = -1
                break
            end
            f0 = f
            if step_control
                L = f + dot(lam, g) + dot(mu, h + z) - gamma * sum(log.(z))
            end
        end
    end
    
    if verbose > 0 && !converged
        println("Did not converge in $i iterations.")
    end
    
    # Package up results
    if eflag != -1
        eflag = converged ? 1 : 0
    end
    
    message = if eflag == 0
        "Did not converge"
    elseif eflag == 1
        "Converged"
    elseif eflag == -1
        "Numerically failed"
    else
        "Please hang up and dial again"
    end
    
    output = Dict(:iterations => i, :hist => hist, :message => message)
    
    # Zero out multipliers on non-binding constraints
    mu[(h .< -feastol) .& (mu .< mu_threshold)] .= 0
    
    # Un-scale cost and prices
    f /= cost_mult
    lam /= cost_mult
    mu /= cost_mult
    
    # Re-package multipliers into struct
    lam_lin = lam[(neqnln+1):neq]
    mu_lin = mu[(niqnln+1):niq]
    kl = findall(lam_lin .< 0)
    ku = findall(lam_lin .> 0)
    
    mu_l = zeros(nx + nA)
    mu_u = zeros(nx + nA)
    
    if !isempty(ieq) && !isempty(kl)
        mu_l[ieq[kl]] = -lam_lin[kl]
    end
    if !isempty(ieq) && !isempty(ku)
        mu_u[ieq[ku]] = lam_lin[ku]
    end
    
    if !isempty(igt)
        mu_l[igt] = mu_lin[nlt+1:nlt+ngt]
    end
    if !isempty(ilt)
        mu_u[ilt] = mu_lin[1:nlt]
    end
    if !isempty(ibx)
        mu_l[ibx] = mu_lin[nlt+ngt+nbx+1:nlt+ngt+2*nbx]
        mu_u[ibx] = mu_lin[nlt+ngt+1:nlt+ngt+nbx]
    end
    
    lambda = Dict(
        :mu_l => mu_l[(nx+1):end],
        :mu_u => mu_u[(nx+1):end],
        :lower => mu_l[1:nx],
        :upper => mu_u[1:nx]
    )
    
    if niqnln > 0
        lambda[:ineqnonlin] = mu[1:niqnln]
    end
    if neqnln > 0
        lambda[:eqnonlin] = lam[1:neqnln]
    end
    
    return x, f, eflag, output, lambda
end

# Convenience function for problem struct input
function mips(nonlinear::NonConvexOPT, ipm::IPM=IPM(1e-6, 100, 1e-6, 1e-6, true, 0.99995))
    # Extract problem dimensions
    x0 = nonlinear.x0
    nx = length(x0)
    
    # Create objective function wrapper that returns (f, df) tuple
    function f_fcn(x)
        f = nonlinear.f(x)
        df = nonlinear.∇f(x)
        return f, df
    end
    
    # Create constraint function wrapper
    function gh_fcn(x)
        # Get equality and inequality constraints
        g = nonlinear.g(x)  # equality constraints
        h = nonlinear.h(x)  # inequality constraints (h(x) ≤ 0)
        
        # Get constraint gradients
        dg = nonlinear.∇g(x)  # gradient of equality constraints
        dh = nonlinear.∇h(x)  # gradient of inequality constraints
        
        return h, g, dh, dg  # Note: MIPS expects transposed gradients
    end
    
    # Create Hessian function wrapper
    function hess_fcn(x, lambda, cost_mult)
        # Extract multipliers
        λ = haskey(lambda, :eqnonlin) ? lambda[:eqnonlin] : Float64[]
        μ = haskey(lambda, :ineqnonlin) ? lambda[:ineqnonlin] : Float64[]
        
        # Ensure correct dimensions
        neq = length(nonlinear.g(x))
        niq = length(nonlinear.h(x))
        
        if length(λ) != neq
            λ = zeros(neq)
        end
        if length(μ) != niq
            μ = zeros(niq)
        end
        
        # Compute Hessian of Lagrangian
        Lxx = nonlinear.Lxx(x, λ, μ)
        return Lxx * cost_mult
    end
    
    # Convert IPM struct to MIPS options dictionary
    opt = Dict{Symbol,Any}(
        :verbose => 0,
        :feastol => ipm.feasible_tol,
        :gradtol => ipm.feasible_tol,
        :comptol => ipm.feasible_tol,
        :costtol => ipm.feasible_tol,
        :max_it => ipm.max_iter,
        :xi => ipm.ξ
    )
    
    # Set up linear constraints (none in this case, but keep structure)
    A = spzeros(0, nx)
    l = Float64[]
    u = Float64[]
    
    # Variable bounds (assume unbounded if not specified)
    xmin = fill(-Inf, nx)
    xmax = fill(Inf, nx)
    
    # Call main MIPS function
    x_opt, f_opt, eflag, output, lambda = mips(f_fcn, x0;
                                               A=A, l=l, u=u,
                                               xmin=xmin, xmax=xmax,
                                               gh_fcn=gh_fcn,
                                               hess_fcn=hess_fcn,
                                               opt=opt)
    
    # Extract multipliers for return
    λ_opt = haskey(lambda, :eqnonlin) ? lambda[:eqnonlin] : Float64[]
    μ_opt = haskey(lambda, :ineqnonlin) ? lambda[:ineqnonlin] : Float64[]
    
    # Create compatible return structure (similar to IPM_Solution)
    return MIPS_Solution(x_opt, λ_opt, μ_opt, f_opt, eflag, output[:iterations], output)
end
