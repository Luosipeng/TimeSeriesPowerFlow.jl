"""
AC Optimal Power Flow formulation implementation.
Full nonlinear AC power flow equations with voltage magnitudes and reactive power.
Following run_acopf.jl exactly, just adapted for unified interface.
"""

# Import PowerFlow constants
using ..PowerFlow: BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN
using ..PowerFlow: F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS
using ..PowerFlow: ANGMIN, ANGMAX, PF, QF, PT, QT
using ..PowerFlow: GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, NCOST, COST
using ..PowerFlow: PQ, PV, REF, NONE

# Include bounds validation functions
include("../bounds_validation.jl")

# ============================================================================
# Problem Creation
# ============================================================================

"""
Create AC OPF problem from case data.
"""
function create_problem(::Type{ACFormulation}, jpc::Dict{String, Any}, options::OPFOptions)
    # Convert options to legacy format for bounds validation
    legacy_opt = Dict{String, Any}(
        "DEFAULT_VMIN" => options.voltage_min_default,
        "DEFAULT_VMAX" => options.voltage_max_default,
        "DEFAULT_PMIN_FRACTION" => 0.0,
        "DEFAULT_PMAX_FRACTION" => 1.0,
        "DEFAULT_QMIN_FRACTION" => -0.5,
        "DEFAULT_QMAX_FRACTION" => 0.5,
        "MIN_POWER_GAP" => options.min_power_gap,
        "MIN_VOLTAGE_GAP" => 0.05,
        "BOUNDS_MARGIN" => options.bounds_margin,
        "DEFAULT_RATE_MULTIPLIER" => 5.0,
        "INIT_PROJECTION_WARN" => 0.1,
        "VERBOSE" => options.verbose
    )
    
    # Validate and prepare data using bounds validation
    validated_jpc = deepcopy(jpc)
    
    # Apply bounds validation following run_acopf.jl pattern
    try
        validated_jpc["busAC"] = validate_voltage_bounds!(validated_jpc["busAC"], legacy_opt)
        validated_jpc["genAC"] = validate_generator_bounds!(validated_jpc["genAC"], validated_jpc["baseMVA"], legacy_opt)
        validated_jpc["branchAC"] = validate_line_limits!(validated_jpc["branchAC"], validated_jpc["baseMVA"], legacy_opt)
        check_problem_feasibility(validated_jpc["busAC"], validated_jpc["genAC"], validated_jpc["branchAC"], validated_jpc["baseMVA"], legacy_opt)
    catch e
        if options.verbose > 0
            println("AC OPF validation failed: $(e)")
        end
        error("AC OPF problem validation failed: $(e)")
    end
    
    # Extract and convert data to internal numbering
    bus, gen, branch, _, _, i2e = PowerFlow.ext2int(
        validated_jpc["busAC"], 
        validated_jpc["genAC"],
        validated_jpc["branchAC"],
        zeros(0, 8),  # load
        zeros(0, 8)   # areas
    )
    
    baseMVA = validated_jpc["baseMVA"]
    nb, ng_total, nl = size(bus, 1), size(gen, 1), size(branch, 1)
    
    # Get bus types and active generators - following run_acopf.jl
    ref_buses = findall(bus[:, BUS_TYPE] .== REF)
    pv_buses = findall(bus[:, BUS_TYPE] .== PV)
    pq_buses = findall(bus[:, BUS_TYPE] .== PQ)
    on = findall(gen[:, GEN_STATUS] .> 0) # active generators
    ng = length(on)
    gbus = Int.(gen[on, GEN_BUS])
    
    if isempty(ref_buses)
        error("No reference bus found")
    end
    
    # Build network matrices - following run_acopf.jl
    Ybus, Yf, Yt = PowerFlow.makeYbus(baseMVA, bus, branch)
    
    # Build generator connection matrix
    Cg = sparse(gbus, 1:ng, ones(ng), nb, ng)
    
    # Variable bounds using validated data
    var_bounds = create_ac_bounds(bus, gen, branch, baseMVA, on, options)
    
    # Cost coefficients
    cost_coeffs = create_ac_costs(gen, on, baseMVA)
    
    # Create problem instance
    return OPFProblem{ACFormulation, Float64}(
        ACFormulation(),
        nb, ng, nl,
        baseMVA,
        Ybus, Yf, Yt,
        Cg,
        Float64.(bus), Float64.(gen), Float64.(branch),
        ref_buses, pv_buses, pq_buses, on,
        var_bounds, cost_coeffs, i2e
    )
end

"""
Create variable bounds for AC formulation using validated data - following run_acopf.jl.
"""
function create_ac_bounds(bus, gen, branch, baseMVA, on, options)
    nb, ng, nl = size(bus, 1), length(on), size(branch, 1)
    
    # Voltage angle bounds (in radians) - following run_acopf.jl
    angle_limit_deg = options.angle_limit_deg
    angle_limit_rad = deg2rad(angle_limit_deg)
    angle_limit_rad = 180.0
    va_min = fill(-angle_limit_rad, nb)
    va_max = fill(angle_limit_rad, nb)
    
    # Reference bus angle is fixed at 0
    ref_buses = findall(bus[:, BUS_TYPE] .== REF)
    va_min[ref_buses] .= 0.0
    va_max[ref_buses] .= 0.0
    
    # Voltage magnitude bounds - more relaxed defaults
    vm_min = bus[:, VMIN]
    vm_max = bus[:, VMAX]
    
    # If bounds are too tight or missing, use more relaxed defaults
    invalid_min = findall((vm_min .<= 0) .| (vm_min .> 0.98))
    invalid_max = findall((vm_max .<= 0) .| (vm_max .< 1.02))
    
    if !isempty(invalid_min)
        vm_min[invalid_min] .= 0.90  # More relaxed than 0.95
    end
    if !isempty(invalid_max)
        vm_max[invalid_max] .= 1.10  # More relaxed than 1.05
    end
    
    # Ensure minimum voltage gap
    tight_v = findall((vm_max - vm_min) .< 0.1)
    if !isempty(tight_v)
        vm_max[tight_v] = vm_min[tight_v] .+ 0.1
    end
    
    # Line flow bounds - set ALL to total load as in run_acopf.jl
    total_load = sum(bus[:, PD]) / baseMVA
    pij_max = fill(total_load, nl)  # Set ALL line limits to total load
    

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
    # Generator bounds - ensure they are feasible
    pg_min = gen[on, PMIN] / baseMVA
    pg_max = gen[on, PMAX] / baseMVA
    qg_min = gen[on, QMIN] / baseMVA
    qg_max = gen[on, QMAX] / baseMVA
    
    # Fix any infeasible generator bounds
    problematic_p = findall(pg_min .>= pg_max)
    if !isempty(problematic_p)
        pg_max[problematic_p] = pg_min[problematic_p] .+ 0.5  # Larger margin
    end
    
    problematic_q = findall(qg_min .>= qg_max)
    if !isempty(problematic_q)
        qg_max[problematic_q] = qg_min[problematic_q] .+ 0.5  # Larger margin
    end
    
    # Ensure adequate generation capacity with generous margin
    total_pd = sum(bus[:, PD]) / baseMVA
    total_qd = sum(bus[:, QD]) / baseMVA
    total_pg_capacity = sum(pg_max)
    total_qg_capacity_pos = sum(qg_max)
    total_qg_capacity_neg = sum(abs.(qg_min))
    
    # Check and warn about capacity
    if total_pd > total_pg_capacity * 0.90  # 10% margin
        @warn "Tight P generation capacity: load=$total_pd, capacity=$total_pg_capacity"
        # Increase generation capacity if needed
        scale_factor = total_pd / (total_pg_capacity * 0.85)
        pg_max .*= scale_factor
    end
    
    if options.verbose > 0
        println("  AC bounds summary:")
        println("    Voltage magnitude: [$(minimum(vm_min)), $(maximum(vm_max))] p.u.")
        println("    Line flow limits: ALL set to $(total_load) p.u.")
        println("    P generation: [$(minimum(pg_min)), $(maximum(pg_max))] p.u.")
        println("    Q generation: [$(minimum(qg_min)), $(maximum(qg_max))] p.u.")
    end
    
    return (
        va_min = va_min, va_max = va_max,
        vm_min = vm_min, vm_max = vm_max,
        pij_max = pij_max,
        pg_min = pg_min, pg_max = pg_max,
        qg_min = qg_min, qg_max = qg_max
    )
end

"""
Create cost coefficients for AC formulation.
"""
function create_ac_costs(gen, on, baseMVA)
    ng = length(on)
    
    # Initialize cost vectors
    c_linear = zeros(ng)
    c_quadratic = zeros(ng)
    
    for (i, gen_idx) in enumerate(on)
        if size(gen, 2) >= NCOST && gen[gen_idx, NCOST] > 0
            ncost = Int(gen[gen_idx, NCOST])
            
            if ncost >= 3 && size(gen, 2) >= COST+1
                # Quadratic cost: c2*Pg^2 + c1*Pg + c0
                c_quadratic[i] = gen[gen_idx, COST] * baseMVA^2
                c_linear[i] = gen[gen_idx, COST + 1] * baseMVA
            elseif ncost == 2 && size(gen, 2) >= COST
                # Linear cost: c1*Pg + c0
                c_linear[i] = gen[gen_idx, COST] * baseMVA
            end
        else
            # Default linear cost
            c_linear[i] = 10.0
        end
    end
    
    return (linear = c_linear, quadratic = c_quadratic)
end

# ============================================================================
# Constraint Functions (EXACTLY following run_acopf.jl)
# ============================================================================

"""
Evaluate AC power balance constraints (EXACTLY from run_acopf.jl).
"""
function evaluate_ac_power_balance(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    # Extract variables following run_acopf.jl variable ordering: [va; vm; pg; qg]
    nb, ng = problem.nb, problem.ng
    va_idx = 1:nb
    vm_idx = nb+1:nb+nb
    pg_idx = nb+nb+1:nb+nb+ng
    qg_idx = nb+nb+ng+1:nb+nb+2*ng
    
    va = x[va_idx]
    vm = x[vm_idx]
    pg = x[pg_idx]
    qg = x[qg_idx]

    # Calculate complex voltages and power injections (EXACTLY from run_acopf.jl)
    V = vm .* exp.(1im * va)
    S_calc = V .* conj.(problem.Ybus * V)
    P_calc = real(S_calc)
    Q_calc = imag(S_calc)
    
    # Load data (per unit) - EXACTLY from run_acopf.jl
    Pd = problem.bus_data[:, PD] / problem.baseMVA
    Qd = problem.bus_data[:, QD] / problem.baseMVA
    
    # Generator injections - EXACTLY from run_acopf.jl
    P_gen = problem.Cg * pg
    Q_gen = problem.Cg * qg
    
    # Power balance equations: P_calc - P_gen + P_load = 0 (EXACTLY from run_acopf.jl)
    P_balance = P_calc - P_gen + Pd
    Q_balance = Q_calc - Q_gen + Qd
    
    return [P_balance; Q_balance]
end

"""
Evaluate AC line flow constraints (EXACTLY from run_acopf.jl).
"""
function evaluate_ac_line_flow(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    # Extract variables following run_acopf.jl
    nb, ng = problem.nb, problem.ng
    va_idx = 1:nb
    vm_idx = nb+1:nb+nb
    
    va = x[va_idx]
    vm = x[vm_idx]
    
    # EXACTLY from run_acopf.jl power_flow_constraints function
    V = vm .* exp.(1im * va)
    f_bus = Int.(problem.branch_data[:, F_BUS])
    Cf = sparse(f_bus, 1:problem.nl, ones(problem.nl), nb, problem.nl)'
    Vf = Cf * V
    Sf = Vf .* conj.(problem.Yf * V)
    
    # Use ALL line limits set to total load (EXACTLY as in run_acopf.jl)
    branch_rate = problem.var_bounds.pij_max  # All set to total load
    
    return abs2.(Sf) - branch_rate.^2
end

"""
Evaluate AC inequality constraints (EXACTLY from run_acopf.jl h_ineq function).
"""
function evaluate_bounds(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    bounds = problem.var_bounds
    
    # Variable bounds (EXACTLY from run_acopf.jl)
    xmin = [bounds.va_min; bounds.vm_min; bounds.pg_min; bounds.qg_min]
    xmax = [bounds.va_max; bounds.vm_max; bounds.pg_max; bounds.qg_max]
    
    # EXACTLY from run_acopf.jl h_ineq function
    h_bounds = [x - xmax; -x + xmin]
    # Add transmission line flow limits (EXACTLY from run_acopf.jl)
    h_flow = evaluate_ac_line_flow(problem, x)
    
    return [h_bounds; h_flow]
end

# ============================================================================
# Jacobian Functions (EXACTLY following run_acopf.jl)
# ============================================================================

"""
Compute AC power balance Jacobian (EXACTLY from run_acopf.jl).
"""
function ac_power_balance_jacobian(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    # Extract variables following run_acopf.jl
    nb, ng = problem.nb, problem.ng
    va_idx = 1:nb
    vm_idx = nb+1:nb+nb
    
    va = x[va_idx]
    vm = x[vm_idx]
    
    # EXACTLY from run_acopf.jl power_balance_jacobian function
    V = vm .* exp.(1im * va)
    V_inverse = 1 ./ V
    Ibus = problem.Ybus * V
    dSbus_dVa = 1im * Diagonal(V) * (Diagonal(conj(Ibus)) - conj(problem.Ybus) * Diagonal(conj(V)))
    dSbus_dVm = Diagonal(V) * (Diagonal(conj(Ibus)) + conj(problem.Ybus) * Diagonal(conj(V))) * Diagonal(V_inverse)
    
    # Real and imaginary parts (EXACTLY from run_acopf.jl)
    dP_dVa = real(dSbus_dVa)
    dP_dVm = real(dSbus_dVm)
    dQ_dVa = imag(dSbus_dVa)
    dQ_dVm = imag(dSbus_dVm)
    
    # Generator derivatives (EXACTLY from run_acopf.jl)
    dP_dpg = -problem.Cg
    dP_dqg = spzeros(nb, ng)
    dQ_dpg = spzeros(nb, ng)
    dQ_dqg = -problem.Cg
    
    # Combine Jacobian blocks (EXACTLY from run_acopf.jl)
    J = [dP_dVa dP_dVm dP_dpg dP_dqg;
         dQ_dVa dQ_dVm dQ_dpg dQ_dqg]
    
    return J'
end

"""
Compute AC line flow Jacobian (EXACTLY from run_acopf.jl).
"""
function ac_line_flow_jacobian(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    # Extract variables following run_acopf.jl
    nb, ng = problem.nb, problem.ng
    va_idx = 1:nb
    vm_idx = nb+1:nb+nb
    
    va = x[va_idx]
    vm = x[vm_idx]
    
    # EXACTLY from run_acopf.jl power_flow_jacobian function
    V = vm .* exp.(1im * va)
    f_bus = Int.(problem.branch_data[:, F_BUS])
    Cf = sparse(f_bus, 1:problem.nl, ones(problem.nl), nb, problem.nl)'
    Vf = Cf * V
    If = problem.Yf * V
    Sf = Vf .* conj.(If)
    
    dSf_dVa = 1im * (Diagonal(conj(If)) * Cf * Diagonal(V) - Diagonal(Cf * V) * conj(problem.Yf) * Diagonal(conj(V)))
    dSf_dVm = Diagonal(conj(If)) * Cf * Diagonal(va) + Diagonal(Cf * V) * conj(problem.Yf) * Diagonal(conj(va))
    dSf_dPg = spzeros(problem.nl, ng)
    dSf_dQg = spzeros(problem.nl, ng)
    
    # Combine all Jacobian blocks (EXACTLY from run_acopf.jl)
    dSf_dx = [dSf_dVa dSf_dVm dSf_dPg dSf_dQg]
    dSf2_dx = 2 * (Diagonal(real(Sf)) * real(dSf_dx) + Diagonal(imag(Sf)) * imag(dSf_dx))
    
    return dSf2_dx'
end

"""
Compute bounds Jacobian (EXACTLY from run_acopf.jl ∇h_ineq function).
"""
function bounds_jacobian(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    nx = length(x)
    # EXACTLY from run_acopf.jl ∇h_ineq function
    grad_bounds = [I(nx); -I(nx)]'
    # Add transmission line flow limits gradient (EXACTLY from run_acopf.jl)
    grad_flow = ac_line_flow_jacobian(problem, x)
    
    return [grad_bounds grad_flow]
end

# ============================================================================
# Hessian Functions
# ============================================================================

"""
Compute AC power balance Hessian for Lagrangian.
"""
function ac_power_balance_hessian(problem::OPFProblem{ACFormulation}, x::Vector{Float64}, lambda::Vector{Float64})
    indices = get_variable_indices(ACFormulation, problem.nb, problem.ng, problem.nl)
    
    va = x[indices.va]
    vm = x[indices.vm]
    
    # Split Lagrange multipliers
    lambda_P = lambda[1:problem.nb]
    lambda_Q = lambda[problem.nb+1:2*problem.nb]
    
    # Complex voltage
    V = vm .* exp.(1im * va)
    V_inv = 1 ./ V
    
    # Second derivatives for power balance constraints
    H_P = compute_power_hessian_real(problem, V, V_inv, lambda_P)
    H_Q = compute_power_hessian_imag(problem, V, V_inv, lambda_Q)
    
    # Total Hessian
    H_total = H_P + H_Q
    
    return H_total
end

"""
Compute real part of power flow Hessian.
"""
function compute_power_hessian_real(problem::OPFProblem{ACFormulation}, V::Vector{ComplexF64}, V_inv::Vector{ComplexF64}, lambda::Vector{Float64})
    nb = problem.nb
    ng = problem.ng
    
    # Current injections
    Ibus = problem.Ybus * V
    
    # Hessian components for real power
    A = Diagonal(lambda) * Diagonal(V)
    B = problem.Ybus * Diagonal(V)
    C = A * conj(B)
    D = conj(problem.Ybus) * Diagonal(V)
    E = Diagonal(conj(V)) * (D * Diagonal(lambda) - Diagonal(D * lambda))
    F = C - A * Diagonal(conj(Ibus))
    G = Diagonal(V_inv)
    
    # Second derivatives
    d2P_dVa2 = real(E + F)
    d2P_dVmdVa = real(1im * G * (E - F))
    d2P_dVadVm = d2P_dVmdVa'
    d2P_dVm2 = real(G * (C + C') * G)
    
    # Assemble full Hessian
    H = [d2P_dVa2 d2P_dVmdVa spzeros(nb, ng) spzeros(nb, ng);
         d2P_dVadVm d2P_dVm2 spzeros(nb, ng) spzeros(nb, ng);
         spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng);
         spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng)]
    
    return H
end

"""
Compute imaginary part of power flow Hessian.
"""
function compute_power_hessian_imag(problem::OPFProblem{ACFormulation}, V::Vector{ComplexF64}, V_inv::Vector{ComplexF64}, lambda::Vector{Float64})
    nb = problem.nb
    ng = problem.ng
    
    # Current injections
    Ibus = problem.Ybus * V
    
    # Hessian components for reactive power
    A = Diagonal(lambda) * Diagonal(V)
    B = problem.Ybus * Diagonal(V)
    C = A * conj(B)
    D = conj(problem.Ybus) * Diagonal(V)
    E = Diagonal(conj(V)) * (D * Diagonal(lambda) - Diagonal(D * lambda))
    F = C - A * Diagonal(conj(Ibus))
    G = Diagonal(V_inv)
    
    # Second derivatives
    d2Q_dVa2 = imag(E + F)
    d2Q_dVmdVa = imag(1im * G * (E - F))
    d2Q_dVadVm = d2Q_dVmdVa'
    d2Q_dVm2 = imag(G * (C + C') * G)
    
    # Assemble full Hessian
    H = [d2Q_dVa2 d2Q_dVmdVa spzeros(nb, ng) spzeros(nb, ng);
         d2Q_dVadVm d2Q_dVm2 spzeros(nb, ng) spzeros(nb, ng);
         spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng);
         spzeros(ng, nb) spzeros(ng, nb) spzeros(ng, ng) spzeros(ng, ng)]
    
    return H
end

"""
Compute AC line flow limits Hessian for Lagrangian.
"""
function ac_line_limits_hessian(problem::OPFProblem{ACFormulation}, x::Vector{Float64}, mu::Vector{Float64})
    indices = get_variable_indices(ACFormulation, problem.nb, problem.ng, problem.nl)
    
    va = x[indices.va]
    vm = x[indices.vm]
    
    # Complex voltage
    V = vm .* exp.(1im * va)
    V_inv = 1 ./ V
    
    # Line quantities
    f_bus = Int.(problem.branch_data[:, F_BUS])
    Cf = sparse(1:problem.nl, f_bus, ones(problem.nl), problem.nl, problem.nb)
    
    Vf = Cf * V
    If = problem.Yf * V
    Sf = Vf .* conj(If)
    
    # Compute second derivatives of line flows
    Af = conj(problem.Yf)' * Diagonal(Diagonal(conj(Sf)) * mu) * Cf
    Bf = Diagonal(conj(V)) * Af * Diagonal(V)
    Df = Diagonal(Af * V) * Diagonal(conj(V))
    Ef = Diagonal(Af' * conj(V)) * Diagonal(V)
    Ff = Bf + Bf'
    Gf = Diagonal(V_inv)
    
    d2Sf_dVa2 = Ff - Df - Ef
    d2Sf_dVmdVa = 1im * Gf * (Bf - Bf' - Df + Ef)
    d2Sf_dVadVm = d2Sf_dVmdVa'
    d2Sf_dVm2 = Gf * Ff * Gf
    
    # Derivatives of line flows
    dSf_dVa = 1im * (Diagonal(conj(If)) * Cf * Diagonal(V) - 
                     Diagonal(Vf) * conj(problem.Yf) * Diagonal(conj(V)))
    dSf_dVm = (Diagonal(conj(If)) * Cf * Diagonal(V_inv) + 
               Diagonal(Vf) * conj(problem.Yf) * Diagonal(conj(V)) * Diagonal(V_inv))
    
    # Hessian of flow magnitude squared
    HSf_dVa2 = 2 * real(d2Sf_dVa2 + dSf_dVa' * Diagonal(mu) * conj(dSf_dVa))
    HSf_dVadVm = 2 * real(d2Sf_dVmdVa + dSf_dVm' * Diagonal(mu) * conj(dSf_dVa))
    HSf_dVmdVa = 2 * real(d2Sf_dVadVm + dSf_dVa' * Diagonal(mu) * conj(dSf_dVm))
    HSf_dVm2 = 2 * real(d2Sf_dVm2 + dSf_dVm' * Diagonal(mu) * conj(dSf_dVm))
    
    # Assemble full Hessian
    H = [HSf_dVa2 HSf_dVmdVa spzeros(problem.nb, problem.ng) spzeros(problem.nb, problem.ng);
         HSf_dVadVm HSf_dVm2 spzeros(problem.nb, problem.ng) spzeros(problem.nb, problem.ng);
         spzeros(problem.ng, problem.nb) spzeros(problem.ng, problem.nb) spzeros(problem.ng, problem.ng) spzeros(problem.ng, problem.ng);
         spzeros(problem.ng, problem.nb) spzeros(problem.ng, problem.nb) spzeros(problem.ng, problem.ng) spzeros(problem.ng, problem.ng)]
    
    return H
end

# ============================================================================
# Solution Update Functions (following run_acopf.jl exactly)
# ============================================================================

"""
Update case data for AC solution (following run_acopf.jl).
"""
function update_ac_case_data!(case_data, problem, solution)
    baseMVA = case_data["baseMVA"]
    
    # Convert back to external numbering
    bus_int = copy(problem.bus_data)
    gen_int = copy(problem.gen_data)
    branch_int = copy(problem.branch_data)
    
    # Update bus data
    bus_int[:, VM] = solution.variables.vm
    bus_int[:, VA] = rad2deg.(solution.variables.va)
    
    # Update generator outputs
    gen_int[problem.active_gens, PG] = solution.variables.pg * baseMVA
    gen_int[problem.active_gens, QG] = solution.variables.qg * baseMVA
    
    # Calculate and update line flows
    V = solution.variables.vm .* exp.(1im * solution.variables.va)
    If = problem.Yf * V
    It = problem.Yt * V
    
    f_bus = Int.(branch_int[:, F_BUS])
    t_bus = Int.(branch_int[:, T_BUS])
    
    Sf = V[f_bus] .* conj(If) * baseMVA
    St = V[t_bus] .* conj(It) * baseMVA
    
    # Ensure branch matrix has enough columns
    if size(branch_int, 2) < QT
        branch_new = zeros(size(branch_int, 1), QT)
        branch_new[:, 1:size(branch_int, 2)] = branch_int
        branch_int = branch_new
    end
    
    # Update line flows
    branch_int[:, PF] = real(Sf)
    branch_int[:, QF] = imag(Sf)
    branch_int[:, PT] = real(St)
    branch_int[:, QT] = imag(St)
    
    # Convert to external numbering
    bus_ext, gen_ext, branch_ext, _ = PowerFlow.int2ext(
        problem.int_to_ext, bus_int, gen_int, branch_int, zeros(0,8), zeros(0,8)
    )
    
    case_data["busAC"] = bus_ext
    case_data["genAC"] = gen_ext
    case_data["branchAC"] = branch_ext
end

"""
Create intelligent initial point for AC optimization using validated bounds.
"""
function create_ac_initial_point(problem::OPFProblem{ACFormulation}, options::OPFOptions)
    indices = get_variable_indices(ACFormulation, problem.nb, problem.ng, problem.nl)
    nx = indices.total
    
    # AC initialization following run_acopf.jl EXACTLY
    # Start with flat voltage profile
    va0 = zeros(problem.nb)
    vm0 = ones(problem.nb)  # Flat 1.0 p.u. voltage magnitudes
    
    # Ensure voltage magnitudes are within bounds
    vm0 = clamp.(vm0, 
                 problem.var_bounds.vm_min .+ 1e-4, 
                 problem.var_bounds.vm_max .- 1e-4)
    
    # Generation initialization - EXACTLY as in run_acopf.jl
    total_pd = sum(problem.bus_data[:, PD]) / problem.baseMVA
    total_qd = sum(problem.bus_data[:, QD]) / problem.baseMVA
    
    # P generation - distribute load equally among generators
    pg0 = fill(total_pd / problem.ng, problem.ng)
    pg0 = clamp.(pg0, problem.var_bounds.pg_min, problem.var_bounds.pg_max)
    
    # Q generation - start at zero or small values
    qg0 = zeros(problem.ng)
    qg0 = clamp.(qg0, problem.var_bounds.qg_min, problem.var_bounds.qg_max)
    
    x0 = zeros(nx)
    x0[indices.va] = va0
    x0[indices.vm] = vm0
    x0[indices.pg] = pg0
    x0[indices.qg] = qg0
    
    return x0
end