# Import necessary modules
using LinearAlgebra, SparseArrays, Statistics
using ..PowerFlow: ext2int, int2ext, makeYbus
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
    runacopf(jpc, opt)

Run AC Optimal Power Flow using Interior Point Method.

Solves the AC optimal power flow problem:
    min sum of generator costs
    subject to:
        AC power balance constraints
        generator limits (P, Q, V)
        line flow limits
        voltage magnitude limits

# Arguments
- `jpc`: Julia Power Case structure (Dict format)
- `opt`: Options dictionary

# Returns
- Updated `jpc` with optimal solution
"""
function runacopf(jpc::Dict{String,Any}, opt::Dict{String,Any})
    # Load solver path
    include(joinpath(dirname(@__DIR__), "solvers", "types.jl"))
    include(joinpath(dirname(@__DIR__), "solvers", "interior_point_method.jl"))
    # Get base MVA
    baseMVA = get(jpc, "baseMVA", 100.0)
    
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
            println("AC OPF validation failed: $(e)")
        end
        jpc["success"] = false
        jpc["f"] = Inf
        return jpc
    end
    
    # Convert to internal numbering
    bus, gen, branch, load, _, i2e = ext2int(bus, gen, branch, load, zeros(0, 8))
    
    # Get dimensions
    nb = size(bus, 1)        # number of buses
    ng = size(gen, 1)        # number of generators
    nl = size(branch, 1)     # number of lines
    
    # Get bus types and indices
    ref = findall(bus[:, BUS_TYPE] .== REF)
    if isempty(ref)
        error("No reference bus found")
    end

    # Get generator data
    on = findall(gen[:, GEN_STATUS] .> 0)
    gbus = Int.(gen[on, GEN_BUS])
    
    # Build admittance matrix
    Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)
    
    # Build generator connection matrix
    Cg = sparse(gbus, 1:length(on), ones(length(on)), nb, length(on))
    # Build transmission line and node connection matrices
    Cf = sparse(branch[:, F_BUS], 1:nl, ones(nl), nb, nl)'
    # Load data (per unit)
    Pd = bus[:, PD] / baseMVA
    Qd = bus[:, QD] / baseMVA
    
    # Set up optimization variables: [va_non_ref; vm; pg; qg]
    ng = length(on)  # number of generators in service

    # Variable indices (FIXED: exclude reference bus angle)
    va_idx = 1:nb                    # non-reference bus angles only
    vm_idx = nb+1:nb+nb             # all voltage magnitudes
    pg_idx = nb+nb+1:nb+nb+ng       # generator active power
    qg_idx = nb+nb+ng+1:nb+nb+2*ng  # generator reactive power
    nx = nb + nb + 2*ng              # total number of variables
    
    # Variable bounds
    # Voltage angles in RADIANS (FIXED)
    angle_limit_deg = get(opt, "ANGLE_LIMIT_DEG", 60.0)
    angle_limit_rad = deg2rad(angle_limit_deg)
    Va_min = fill(-angle_limit_rad, nb)
    Va_max = fill(angle_limit_rad, nb)
    Va_min[ref] .= 0.0  # Reference bus angle is fixed at 0
    Va_max[ref] .= 0.0  # Reference bus angle is fixed at 0
    # Voltage magnitude limits from bus data if available
    Vm_min = bus[:, VMIN]
    Vm_max = bus[:, VMAX]
    
    # If limits not specified, use defaults
    if all(Vm_min .== 0)
        Vm_min = 0.9 * ones(nb)
    end
    if all(Vm_max .== 0)
        Vm_max = 1.1 * ones(nb)
    end
    
    # Generator power limits from actual data (FIXED)
    if all(gen[on, PMIN] .== 0)
        gen[on, PMIN] .= 0.0
    end
    if all(gen[on, PMAX] .== 0)
        gen[on, PMAX] .= 1.0 * baseMVA  # Set default max if not specified
    end
    if all(gen[on, QMIN] .== 0)
        gen[on, QMIN] .= -1.0 * baseMVA  # Set default min if not specified
    end
    if all(gen[on, QMAX] .== 0)
        gen[on, QMAX] .= 1.0 * baseMVA   # Set default max if not specified
    end
    
    # Combine bounds
    Pg_min = gen[on, PMIN] / baseMVA
    Pg_max = gen[on, PMAX] / baseMVA
    Qg_min = gen[on, QMIN] / baseMVA
    Qg_max = gen[on, QMAX] / baseMVA
    
    # Validate bounds consistency
    if any(Pg_min .>= Pg_max)
        problematic = findall(Pg_min .>= Pg_max)
        Pg_max[problematic] = Pg_min[problematic] .+ 0.1
    end
    
    if any(Qg_min .>= Qg_max)
        problematic = findall(Qg_min .>= Qg_max)
        Qg_max[problematic] = Qg_min[problematic] .+ 0.1
    end
    total_Pd = sum(Pd)
    # Transmission line flow limits
    # If branch data is available, validate flow limits
    if size(branch, 2) >= RATE_A
        # Use RATE_A as flow limit if available
        branch_rate = branch[:, RATE_A] / baseMVA
        if all(branch_rate .== 0)
            # If no limits specified, use a large default value
            branch_rate = total_Pd * ones(nl)
        end
    else
        # Default to a large value if no limits specified
        branch_rate = total_Pd * ones(nl)
    end
    branch_rate = total_Pd * ones(nl) # To do: with practical limits
    
    # Set variable bounds
    xmin = [Va_min; Vm_min; Pg_min; Qg_min]
    xmax = [Va_max; Vm_max; Pg_max; Qg_max]
    
    # INTELLIGENT INITIALIZATION STRATEGY
    # Start with flat voltage profile but respect bounds
    Va0 = zeros(nb)
    Vm0 = clamp.(ones(nb), Vm_min .+ 1e-3, Vm_max .- 1e-3)
    
    # Smart generation initialization based on load distribution
    total_Pd = sum(Pd)
    total_Qd = sum(Qd)
    
    # Distribute P generation based on capacity and proximity to load
    if total_Pd > 0
        Pg_capacity = Pg_max - Pg_min
        capacity_share = Pg_capacity / sum(Pg_capacity)
        Pg0 = Pg_min + capacity_share * total_Pd * 1.05  # 5% margin for losses
    else
        Pg0 = (Pg_min + Pg_max) / 2
    end
    
    # Distribute Q generation
    if total_Qd > 0 && sum(Qg_max - Qg_min) > 0
        Qg_capacity = Qg_max - Qg_min
        capacity_share = Qg_capacity / sum(Qg_capacity)
        Qg0 = Qg_min + capacity_share * total_Qd * 1.1  # 10% margin
    else
        Qg0 = (Qg_min + Qg_max) / 2
    end
    
    # Ensure bounds are satisfied
    Pg0 = clamp.(Pg0, Pg_min, Pg_max)
    Qg0 = clamp.(Qg0, Qg_min, Qg_max)
    
    x0 = [Va0; Vm0; Pg0; Qg0]
    
    # VALIDATE INITIAL POINT
    x0 = validate_feasible_initialization(x0, xmin, xmax, opt)
    
    # Build objective function coefficients
    c = zeros(nx)
    H = spzeros(nx, nx)
    
    # Generator cost coefficients
    for i in 1:ng
        gen_idx = on[i]
        if size(gen, 2) >= NCOST && gen[gen_idx, NCOST] > 0
            ncost = Int(gen[gen_idx, NCOST])
            
            if ncost >= 3 && size(gen, 2) >= COST+1
                # Quadratic cost: c2*Pg^2 + c1*Pg + c0
                c2 = gen[gen_idx, COST] * baseMVA^2
                c1 = gen[gen_idx, COST + 1] * baseMVA
                # Set up quadratic term in H matrix
                H[pg_idx[i], pg_idx[i]] = 2 * c2
                # Set up linear term in c vector
                c[pg_idx[i]] = c1
            elseif ncost == 2 && size(gen, 2) >= COST
                # Linear cost: c1*Pg + c0
                c1 = gen[gen_idx, COST] * baseMVA
                c[pg_idx[i]] = c1
            else
                # Default linear cost
                c[pg_idx[i]] = 10.0
            end
        else
            # Default linear cost if no cost data
            c[pg_idx[i]] = 10.0
        end
    end
    
    # Power balance equations - FIXED to handle reference bus correctly
    function power_balance_constraints(x)
        va = x[va_idx]  # Non-reference angles
        vm = x[vm_idx]
        pg = x[pg_idx]
        qg = x[qg_idx]

        # Calculate complex voltages and power injections
        V = vm .* exp.(1im * va)
        S_calc = V .* conj.(Ybus * V)
        P_calc = real(S_calc)
        Q_calc = imag(S_calc)
        
        # Generator injections
        P_gen = Cg * pg
        Q_gen = Cg * qg
        
        # Standard power balance: P_calc - P_gen + P_load = 0
        P_balance = P_calc - P_gen + Pd
        Q_balance = Q_calc - Q_gen + Qd
        
        return [P_balance; Q_balance]
    end
    
    # Jacobian of power balance constraints - FIXED to handle reference bus
    function power_balance_jacobian(x)
        va = x[va_idx]  
        vm = x[vm_idx]
        V = vm .* exp.(1im * va)  # Complex voltages
        V_inverse = 1 ./ V  # Inverse voltages for current calculations
        Ibus = Ybus * V  # Complex current injections
        dSbus_dVa = 1im * Diagonal(V) * ( Diagonal(conj(Ibus)) - conj(Ybus)*Diagonal(conj(V)) )
        dSbus_dVm = Diagonal(V) * (Diagonal(conj(Ibus)) + conj(Ybus)*Diagonal(conj(V))) * Diagonal(V_inverse)
        # ∂P/∂θ block (nb × nb)
        dP_dVa = real(dSbus_dVa)
        # ∂P/∂V block (nb × nb)
        dP_dVm = real(dSbus_dVm)
        # ∂Q/∂θ block (nb × nb)
        dQ_dVa = imag(dSbus_dVa)
        # ∂Q/∂V block (nb × nb)
        dQ_dVm = imag(dSbus_dVm)
        # Generator derivatives (negative because P_gen is subtracted)
        dP_dpg = -Cg
        dP_dqg = spzeros(nb, ng)
        dQ_dpg = spzeros(nb, ng)
        dQ_dqg = -Cg
        # Combine all Jacobian blocks
        J = [dP_dVa dP_dVm dP_dpg dP_dqg;
             dQ_dVa dQ_dVm dQ_dpg dQ_dqg]
        
        return J'
    end
    # 
    function power_balance_second_order_derivative(x, lambda)
        # This function computes the second derivative of the power flow equations
        # with respect to the optimization variables.
        # It is used in the Hessian of the Lagrangian.
        va = x[va_idx]
        vm = x[vm_idx]
        V = vm .* exp.(1im * va)  # Complex voltages
        V_inverse = 1 ./ V  # Inverse voltages for current calculations
        Ibus = Ybus * V  # Complex current injections
        A = Diagonal(lambda) * Diagonal(V)
        B = Ybus * Diagonal(V)
        C = A * conj(B)
        D = conj(Ybus) * Diagonal(V)
        E = Diagonal(conj(V))*(D*Diagonal(lambda) - Diagonal(D*lambda))
        F = C - A*Diagonal(conj(Ibus))
        G = Diagonal(V_inverse)
        dG2_dVa2 = E + F
        dG2_dVmdVa = 1im * G * (E-F)
        dG2_dVadVm = dG2_dVmdVa'
        dG2_dVm2 = G * (C + C')*G
        dG2 = [
            dG2_dVa2 dG2_dVmdVa spzeros(nb, ng) spzeros(nb, ng);
            dG2_dVadVm dG2_dVm2 spzeros(nb, ng) spzeros(nb, ng);
            spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng);
            spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng)
        ]
        return dG2
    end
    # We need to compute the Hessian of the power balance constraints
    function power_balance_hessian(x, lambda, mu)
        lambda_P = lambda[1:nb]
        lambda_Q = lambda[nb+1:2*nb]
        H_eq_P = real(power_balance_second_order_derivative(x, lambda_P))
        H_eq_Q = imag(power_balance_second_order_derivative(x, lambda_Q))
        # Combine Hessian blocks into a single matrix
        H_eq = H_eq_P + H_eq_Q
        return H_eq
    end

    function power_flow_constraints(x)
        # This function computes the power balance constraints
        # It is used in the equality constraints of the optimization problem.
        va = x[va_idx]
        vm = x[vm_idx]
        V = vm .* exp.(1im * va)  # Complex voltages
        Vf = Cf * V  # Voltage at the sending end of lines
        Sf = Vf .* conj.(Yf * V)  # Complex power at sending end
        return abs2.(Sf) - branch_rate.^2  # Power flow magnitude constraints
    end
    # Jacobian of power flow constraints
    function power_flow_jacobian(x)
        va = x[va_idx]
        vm = x[vm_idx]
        V = vm .* exp.(1im * va)  # Complex voltages
        Vf = Cf * V  # Voltage at the sending end of lines
        If = Yf * V  # Current at the sending end of lines
        Sf = Vf .* conj(If)  # Complex power at sending end
        dSf_dVa = 1im *( Diagonal(conj(If)) * Cf * Diagonal(V) - Diagonal( Cf * V) * conj(Yf) * Diagonal(conj(V)) )
        dSf_dVm = Diagonal(conj(If)) * Cf * Diagonal(va) + Diagonal(Cf * V) * conj(Yf) * Diagonal(conj(va))
        dSf_dPg = spzeros(nl, ng)
        dSf_dQg = spzeros(nl, ng)
        # Combine all Jacobian blocks
        dSf_dx = [dSf_dVa dSf_dVm dSf_dPg dSf_dQg;]
        dSf2_dx = 2*(Diagonal(real(Sf))*real(dSf_dx) + Diagonal(imag(Sf))*imag(dSf_dx))
        return dSf2_dx'
    end
    # Hessian matrix of power flow constraints
    function power_flow_second_order_derivative(x, mu)       
        va = x[va_idx]
        vm = x[vm_idx]
        V = vm .* exp.(1im * va)  # Complex voltages
        Vf = Cf * V  # Voltage at the sending end of lines
        If = Yf * V  # Current at the sending end of lines
        Sf = Vf .* conj(If)  # Complex power at sending end
        # Compute the second derivatives
        Af = conj(Yf)' * Diagonal(Diagonal(conj(Sf)) * mu) * Cf
        Bf = Diagonal(conj(V)) * Af * Diagonal(V)
        Df = Diagonal(Af*V)*Diagonal(conj(V))
        Ef = Diagonal(Af'*conj(V))*Diagonal(V)
        Ff = Bf + Bf'
        Gf = Diagonal(1.0 ./vm)
        dSf_dVa2 = Ff - Df - Ef
        dSf_dVmdVa = 1im * Gf * (Bf - Bf' - Df + Ef)
        dSf_dVadVm = dSf_dVmdVa'
        dSf_dVm2 = Gf * Ff * Gf

        dSf_dVa = 1im *( Diagonal(conj(If)) * Cf * Diagonal(V) - Diagonal( Cf * V) * conj(Yf) * Diagonal(conj(V)) )
        dSf_dVm = Diagonal(conj(If)) * Cf * Diagonal(va) + Diagonal(Cf * V) * conj(Yf) * Diagonal(conj(va))

        HSf_dVa2 = 2 * real(dSf_dVa2 + dSf_dVa'*Diagonal(mu)*conj(dSf_dVa))
        HSf_dVadVm = 2 * real(dSf_dVmdVa + dSf_dVm'*Diagonal(mu)*conj(dSf_dVa))
        HSf_dVmdVa = 2 * real(dSf_dVadVm + dSf_dVa'*Diagonal(mu)*conj(dSf_dVm))
        HSf_dVm2 = 2 * real(dSf_dVm2 + dSf_dVm'*Diagonal(mu)*conj(dSf_dVm))

        Hsf = [
            HSf_dVa2 HSf_dVmdVa spzeros(nb, ng) spzeros(nb, ng);
            HSf_dVadVm HSf_dVm2 spzeros(nb, ng) spzeros(nb, ng);
            spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng);
            spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng)
        ]
        return Hsf
    end
    
    # Combined equality constraints g(x) = 0
    g_eq = power_balance_constraints
    
    # Gradient of equality constraints
    ∇g_eq = power_balance_jacobian
    
    # Inequality constraints h(x) ≤ 0
    h_ineq = function(x)
        # Variable bounds
        h_bounds = [x - xmax; -x + xmin]
        # Add transmission line flow limits
        h_bounds = [h_bounds; power_flow_constraints(x)]
        return h_bounds
    end
    
    # Gradient of inequality constraints
    ∇h_ineq = function(x)
        grad_bounds = [I(nx); -I(nx)]'
        # Add transmission line flow limits gradient
        # return grad_bounds
        return [grad_bounds power_flow_jacobian(x)]
    end
    
    # Objective function f(x) = 0.5 * x' * H * x + c' * x
    f_obj = function(x)
        if !isempty(H.nzval)
            return (0.5 * x' * H * x + c' * x)[1]
        else
            return (c' * x)[1]
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
    
    # Lagrangian gradient
    Lx = function(x, λ, μ)
        grad_f = ∇f_obj(x)
        grad_g = ∇g_eq(x)
        grad_h = ∇h_ineq(x)
        
        return grad_f + grad_g * λ + grad_h * μ
    end
    
    # Hessian of Lagrangian (simplified - only objective contributes)
    Lxx = function(x, λ, μ)
        return ∇2f_obj(x) + power_balance_hessian(x, λ, μ) + power_flow_second_order_derivative(x, μ[2 * nx + 1:end])
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
    
    # ENHANCED SOLVER SETTINGS
    solver_tol = get(opt, "OPF_VIOLATION", 5e-5)  # Reasonable tolerance for AC OPF
    max_iterations = get(opt, "OPF_MAX_IT", 200)   # More iterations for nonlinear problem
    
    ipm = IPM(
        solver_tol,        # tol
        max_iterations,    # max_iter
        1e-4,             # positive_tol (relaxed for AC OPF)
        solver_tol,       # feasible_tol
        false,            # initial_point_projection
        0.99              # ξ (more conservative for nonlinear)
    )
    
    # Initialize result variables
    success = false
    x_sol = x0
    f_val = Inf
    result = nothing  # Initialize result variable
    
    # ROBUST SOLUTION ATTEMPT
    try
        result = interior_point_method(nonlinear, ipm)
        success = result.eflag == 1
        x_sol = result.x
        f_val = result.obj
        
        if success
            # SOLUTION VALIDATION
            constraint_violation = maximum(abs.(g_eq(x_sol)))
            if constraint_violation > solver_tol * 10
                @warn "Large constraint violation in solution: $(constraint_violation)"
                success = false
            end
        end
        
    catch e
        if get(opt, "VERBOSE", 1) > 0
            println("AC OPF solver error: $(e)")
        end
        success = false
        result = nothing
    end
    
    if success
        # Extract solution variables
        Va_sol = x_sol[va_idx]
        Vm_sol = x_sol[vm_idx]
        Pg_sol = x_sol[pg_idx]
        Qg_sol = x_sol[qg_idx]
        
        # Update bus voltage magnitudes and angles
        bus[:, VM] = Vm_sol
        bus[:, VA] = rad2deg.(Va_sol)
        
        # Update generator outputs
        gen[on, PG] = Pg_sol * baseMVA
        gen[on, QG] = Qg_sol * baseMVA
                
        # Calculate line flows
        V = Vm_sol .* exp.(1im * Va_sol)
        If = Yf * V
        It = Yt * V
        Sf = V[Int.(branch[:, F_BUS])] .* conj.(If) * baseMVA
        St = V[Int.(branch[:, T_BUS])] .* conj.(It) * baseMVA
                
        # Ensure branch matrix has enough columns
        if size(branch, 2) < QT
            branch_new = zeros(size(branch, 1), QT)
            branch_new[:, 1:size(branch, 2)] = branch
            branch = branch_new
        end
                
        # Update line flows
        branch[:, PF] = real(Sf)
        branch[:, QF] = imag(Sf)
        branch[:, PT] = real(St)
        branch[:, QT] = imag(St)
        
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
        jpc["hist"] = result.hist
        jpc["f"] = f_val 
        jpc["mis_match"] = g_eq(x_sol)
        
        if get(opt, "VERBOSE", 1) > 0
            println("AC OPF converged using Interior Point Method")
            println("Objective function value: $(f_val)")
            println("Iterations: $(result.iterations)")
            println("Final constraint violation: $(maximum(abs.(g_eq(x_sol))))")
        end
    else
        # Update jpc for failure case
        jpc["success"] = false
        jpc["f"] = Inf
        jpc["iterations"] = result !== nothing ? result.iterations : 0
        
        if get(opt, "VERBOSE", 1) > 0
            println("AC OPF failed to converge")
            if result !== nothing
                println("Iterations: $(result.iterations)")
                println("Exit flag: $(result.eflag)")
            else
                println("Solver failed to initialize")
            end
        end
    end
    
    return jpc
end
