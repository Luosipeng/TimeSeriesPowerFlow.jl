# Import necessary modules
using LinearAlgebra, SparseArrays
using ..PowerFlow: ext2int, int2ext, bustypes
using ..PowerFlow: PQ, PV, REF, NONE
using ..PowerFlow: BUS_I, BUS_TYPE, PD, QD, GS, BS
using ..PowerFlow: BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN
using ..PowerFlow: F_BUS, T_BUS, BR_R, BR_X, BR_B
using ..PowerFlow: RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS
using ..PowerFlow: ANGMIN, ANGMAX, PF, QF, PT, QT
using ..PowerFlow: GEN_BUS, PG, QG, QMAX, QMIN
using ..PowerFlow: VG, MBASE, GEN_STATUS, PMAX, PMIN, NCOST, COST

include("bounds_validation.jl")

"""
    rundcopf(jpc, opt)

Run DC Optimal Power Flow using Interior Point Method.

Solves the DC optimal power flow problem:
    min sum of generator costs
    subject to:
        power balance constraints
        generator limits
        line flow limits
        voltage angle limits

# Arguments
- `jpc`: Julia Power Case structure (Dict format)
- `opt`: Options dictionary

# Returns
- Updated `jpc` with optimal solution
"""
function rundcopf(jpc::Dict{String,Any}, opt::Dict{String,Any})
    # Load solver path
    include(joinpath(dirname(@__DIR__), "solvers", "types.jl"))
    include(joinpath(dirname(@__DIR__), "solvers", "interior_point_method.jl"))
    # Define helper function for sparse identity matrix  
    speye(n) = sparse(I, n, n)
    
    # Get base MVA
    baseMVA = jpc["baseMVA"]
    
    # Get matrices
    bus = copy(jpc["busAC"])
    gen = copy(jpc["genAC"])
    branch = copy(jpc["branchAC"])
    load = haskey(jpc, "loadAC") ? copy(jpc["loadAC"]) : zeros(0, 8)
    
    # ROBUST BOUNDS VALIDATION
    try
        bus = validate_voltage_bounds!(bus, opt)
        gen = validate_generator_bounds!(gen, baseMVA, opt)
        branch = validate_line_limits!(branch, baseMVA, opt)
        check_problem_feasibility(bus, gen, branch, baseMVA, opt)
    catch e
        if get(opt, "VERBOSE", 1) > 0
            println("Problem validation failed: $(e)")
        end
        jpc["success"] = false
        jpc["f"] = Inf
        return jpc
    end
    
    # Convert to internal numbering using matrix-based signature
    bus, gen, branch, load, _, i2e = ext2int(bus, gen, branch, load, zeros(0, 8))
    
    # Get dimensions
    nb = size(bus, 1)        # number of buses
    ng = size(gen, 1)        # number of generators
    nl = size(branch, 1)     # number of lines
    
    # Get bus types
    ref, pv, pq = bustypes(bus, gen)
    # find the reference bus indexes
    ref = findall(bus[:, BUS_TYPE] .== REF)
    
    # Get generator and load data
    on = findall(gen[:, GEN_STATUS] .> 0)
    gbus = gen[on, GEN_BUS]
    
    # Build generator connection matrix
    Cg = sparse(Int.(gbus), 1:length(on), ones(length(on)), nb, length(on))  # Fix: convert to Int
    # Build transmission from and to matrices
    f = Int.(branch[:, F_BUS])                      # list of "from" buses
    t = Int.(branch[:, T_BUS])                      # list of "to" buses
    i = [1:nl; 1:nl]                               # double set of row indices
    Cft = sparse(i, [f; t], [ones(nl); -ones(nl)], nl, nb)  # connection matrix for the electrical network

    # Set up optimization variables: [va; pij; pg]
    ng = length(on)  # number of generators in service
    va = 0;
    pij = va + nb;
    pg = pij + nl;
    nx = pg + ng; 
    # Variable bounds - Fix: handle angle limits properly
    # voltage angle limits
    angle_limit_deg = get(opt, "ANGLE_LIMIT_DEG", 60.0)  # More conservative default
    Va_min = fill(-angle_limit_deg, nb)
    Va_max = fill(angle_limit_deg, nb)
    Va_min[ref] .= 0.0  # Reference bus angle is fixed at 0
    Va_max[ref] .= 0.0  # Reference bus angle is fixed at 0
    # transmission line flow limits
    Pij_max = zeros(nl)
    if size(branch, 2) >= RATE_A
        Pij_max = branch[:, RATE_A] / baseMVA  # line flow limits in p.u.
        # More sophisticated default limit calculation
        zero_limits = findall(Pij_max .<= 0)
        if !isempty(zero_limits)
            # Calculate default based on system characteristics
            total_gen_capacity = sum(gen[on, PMAX]) / baseMVA
            avg_line_capacity = total_gen_capacity / nl * get(opt, "AVG_LINE_FACTOR", 2.0)
            Pij_max[zero_limits] .= avg_line_capacity
            
            if get(opt, "VERBOSE", 1) > 0
                println("Set $(length(zero_limits)) missing line limits to $(avg_line_capacity) p.u.")
            end
        end
    else
        # Intelligent default based on impedance
        for i in 1:nl
            x = branch[i, BR_X]
            if x > 0
                # Capacity inversely related to reactance
                Pij_max[i] = min(5.0 / x, 10.0)  # Cap at 10 p.u.
            else
                Pij_max[i] = 2.0  # Conservative default
            end
        end
    end

    Pij_min = -Pij_max  # assuming symmetric limits
    # generator power limits
    Pg_min = max.(gen[on, PMIN] / baseMVA, 0.0)  # Ensure non-negative minimum
    Pg_max = gen[on, PMAX] / baseMVA
    
    # Ensure feasible generation bounds
    if any(Pg_min .>= Pg_max)
        problematic = findall(Pg_min .>= Pg_max)
        for i in problematic
            Pg_max[i] = Pg_min[i] + 0.1  # Add small margin
        end
    end
    
    # Combine bounds
    xmin = [Va_min; Pij_min; Pg_min]
    xmax = [Va_max; Pij_max; Pg_max]
    
    # IMPROVED INITIALIZATION
    Va0 = zeros(nb)
    Pij0 = zeros(nl)
    
    # Better initial generation distribution
    total_load = sum(bus[:, PD]) / baseMVA
    gen_capacity = Pg_max - Pg_min
    if sum(gen_capacity) > 0
        # Distribute load based on capacity
        Pg0 = Pg_min + gen_capacity * (total_load / sum(gen_capacity))
    else
        Pg0 = (Pg_min + Pg_max) / 2
    end
    
    x0 = [Va0; Pij0; Pg0]
    
    # VALIDATE AND PROJECT INITIAL POINT
    x0 = validate_feasible_initialization(x0, xmin, xmax, opt)
    
    # Build objective function coefficients - Fix: handle cost structure properly
    c = zeros(nx)
    H = spzeros(nx, nx)
    
    # Generator cost coefficients
    npg = length(on)  # number of active generators
    for i in 1:npg
        gen_idx = on[i]
        if size(gen, 2) >= NCOST && gen[gen_idx, NCOST] > 0
            ncost = Int(gen[gen_idx, NCOST])
            
            if ncost >= 3 && size(gen, 2) >= COST+1
                # Quadratic cost: c2*Pg^2 + c1*Pg + c0
                c2 = gen[gen_idx, COST] * baseMVA^2      # quadratic coefficient
                c1 = gen[gen_idx, COST + 1] * baseMVA      # linear coefficient
                # Set up quadratic term in H matrix
                H[pg + i, pg + i] = 2 * c2
                # Set up linear term in c vector
                c[pg + i] = c1
            elseif ncost == 2 && size(gen, 2) >= COST
                # Linear cost: c1*Pg + c0
                c1 = gen[gen_idx, COST] * baseMVA
                c[pg + i] = c1
            else
                # Default linear cost
                c[pg + i] = 10.0  # Default cost per MW
            end
        else
            # Default linear cost
            c[pg + i] = 10.0  # Default cost per MW
        end
    end

    # Load injection at each bus
    Pd = bus[:, PD] / baseMVA
    Gs = bus[:, GS] / baseMVA
    Pinj = Pd + Gs
    
    # Equality constraint: 
    # Cft * Pij = Cg * Pg - Pinj
    Aeq = [spzeros(nb, nb) -Cft' Cg]
    beq = Pinj

    # Pij = diag(1./(branch(:,BR_X)))*Cft
    # Build branch susceptance matrix
    Aeq_power_flow = spzeros(nl, nx)
    beq_power_flow = zeros(nl)
    Aeq_power_flow[:, va + 1:va + nb] = -Diagonal(1 ./ branch[:, BR_X]) * Cft
    Aeq_power_flow[:, pij + 1:pij + nl] = I(nl)
    Aeq = [Aeq; Aeq_power_flow]
    beq = [beq; beq_power_flow]

    # Transform inequality constraints h(x) <= 0
    # Variable bounds: xmin <= x <= xmax becomes
    # x - xmax <= 0 and -x + xmin <= 0
    h_ineq = function(x)
        return [x - xmax; -x + xmin]  # h(x) <= 0
    end
    
    # Gradient of inequality constraints
    ∇h_ineq = function(x)
        return [I(nx); -I(nx)]'  # ∇h returns matrix where each column is ∇h_i
    end
    
    # Transform equality constraints g(x) = 0
    g_eq = function(x)
        return Aeq * x - beq  # g(x) = 0
    end
    
    # Gradient of equality constraints
    ∇g_eq = function(x)
        return Aeq'  # ∇g returns matrix where each column is ∇g_i
    end
    
    # Objective function f(x) = 0.5 * x' * H * x + c' * x
    f_obj = function(x)
        if !isempty(H.nzval)
            return (0.5 * x' * H * x + c' * x)[1]  # Ensure scalar return
        else
            return (c' * x)[1]  # Ensure scalar return  
        end
    end
    
    # Gradient of objective function
    ∇f_obj = function(x)
        if !isempty(H.nzval)
            return H * x + c
        else
            return c
        end
    end
    
    # Hessian of objective function
    ∇2f_obj = function(x)
        if !isempty(H.nzval)
            return H
        else
            return spzeros(nx, nx)
        end
    end
    
    # Lagrangian: L(x,λ,μ) = f(x) + λ'g(x) + μ'h(x)
    Lx = function(x, λ, μ)
        return ∇f_obj(x) + ∇g_eq(x) * λ + ∇h_ineq(x) * μ
    end
    
    # Hessian of Lagrangian (for DC OPF, only objective contributes)
    Lxx = function(x, λ, μ)
        return ∇2f_obj(x)  # g and h are linear, so their Hessians are zero
    end
    
    # Set up the NonConvexOPT problem structure
    nonlinear = NonConvexOPT(
        f_obj,      # objective function
        ∇f_obj,     # gradient of objective
        ∇2f_obj,    # Hessian of objective
        g_eq,       # equality constraints g(x) = 0
        ∇g_eq,      # gradient of equality constraints
        h_ineq,     # inequality constraints h(x) <= 0
        ∇h_ineq,    # gradient of inequality constraints
        Lx,         # gradient of Lagrangian
        Lxx,        # Hessian of Lagrangian
        x0          # initial point
    )
    
    # Set up IPM solver options
    solver_tol = get(opt, "OPF_VIOLATION", 1e-5)  # Slightly relaxed default
    max_iterations = get(opt, "OPF_MAX_IT", 150)   # More iterations
    
    ipm = IPM(
        solver_tol,        # tol
        max_iterations,    # max_iter
        1e-8,             # positive_tol
        solver_tol,       # feasible_tol
        false,            # initial_point_projection
        0.995             # ξ (more conservative step size)
    )
    
    # Solve using Interior Point Method
    result = interior_point_method(nonlinear, ipm)
    
    # Extract solution
    success = result.eflag == 1
    x_sol = result.x
    f_val = result.obj
    
    if success
        # Extract voltage angles and generator outputs  
        Va_sol = x_sol[1:nb]
        Pij_sol = x_sol[nb+1:nb+nl]
        Pg_sol = x_sol[nb+nl+1:end]
        
        # Update bus voltage angles
        bus[:, VA] = Va_sol * 180/pi
        
        # Update generator outputs
        gen[on, PG] = Pg_sol * baseMVA
        
        # Ensure branch matrix has enough columns for power flow results
        if size(branch, 2) < PT
            # Expand branch matrix to include power flow columns
            branch_new = zeros(size(branch, 1), max(PT, QT))
            branch_new[:, 1:size(branch, 2)] = branch
            branch = branch_new
        end
        
        # Update line flows
        branch[:, PF] = Pij_sol * baseMVA
        branch[:, PT] = -Pij_sol * baseMVA
        if size(branch, 2) >= QT
            branch[:, QF] .= 0.0
            branch[:, QT] .= 0.0
        end
        
        if get(opt, "VERBOSE", 1) > 0
            println("DC OPF converged using Interior Point Method")
            println("Objective function value: $(f_val)")
            println("Iterations: $(result.iterations)")
        end
    else
        if get(opt, "VERBOSE", 1) > 0
            println("DC OPF failed to converge: exitflag = $(result.eflag)")
            println("Iterations: $(result.iterations)")
            # Additional diagnostic information
            if haskey(result, "hist") && !isempty(result.hist.x_record)
                final_violation = norm(g_eq(result.x), Inf)
                println("Final constraint violation: $(final_violation)")
            end
        end
    end
    
    # Convert back to external numbering
    bus_ext, gen_ext, branch_ext, load_ext = int2ext(i2e, bus, gen, branch, load, zeros(0, 8))
    
    # Update jpc
    jpc["busAC"] = bus_ext
    jpc["genAC"] = gen_ext
    jpc["branchAC"] = branch_ext
    if !isempty(load_ext)
        jpc["loadAC"] = load_ext
    end
    jpc["success"] = success
    jpc["iterations"] = result.iterations  # Use the iterations field directly
    jpc["f"] = success ? f_val : Inf
    
    return jpc
end