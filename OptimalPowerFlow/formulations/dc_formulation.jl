"""
DC Optimal Power Flow formulation implementation.
Linear approximation with voltage angles and active power only.
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
Create DC OPF problem from case data.
"""
function create_problem(::Type{DCFormulation}, jpc::Dict{String, Any}, options::OPFOptions)
    # Convert options to legacy format for bounds validation
    legacy_opt = Dict{String, Any}(
        "DEFAULT_VMIN" => options.voltage_min_default,
        "DEFAULT_VMAX" => options.voltage_max_default,
        "DEFAULT_PMIN_FRACTION" => 0.0,
        "DEFAULT_PMAX_FRACTION" => 1.0,
        "DEFAULT_QMIN_FRACTION" => -0.5,
        "DEFAULT_QMAX_FRACTION" => 0.5,
        "MIN_POWER_GAP" => options.min_power_gap,
        "BOUNDS_MARGIN" => options.bounds_margin,
        "DEFAULT_RATE_MULTIPLIER" => 5.0,
        "INIT_PROJECTION_WARN" => 0.1,
        "VERBOSE" => options.verbose
    )
    
    # Validate and prepare data using bounds validation
    validated_jpc = deepcopy(jpc)
    
    # Apply bounds validation
    try
        validated_jpc["busAC"] = validate_voltage_bounds!(validated_jpc["busAC"], legacy_opt)
        validated_jpc["genAC"] = validate_generator_bounds!(validated_jpc["genAC"], validated_jpc["baseMVA"], legacy_opt)
        validated_jpc["branchAC"] = validate_line_limits!(validated_jpc["branchAC"], validated_jpc["baseMVA"], legacy_opt)
        check_problem_feasibility(validated_jpc["busAC"], validated_jpc["genAC"], validated_jpc["branchAC"], validated_jpc["baseMVA"], legacy_opt)
    catch e
        if options.verbose > 0
            println("Problem validation failed: $(e)")
        end
        error("DC OPF problem validation failed: $(e)")
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
    
    # Get active generators only (following run_dcopf.jl)
    on = findall(gen[:, GEN_STATUS] .> 0)
    ng = length(on)
    gbus = Int.(gen[on, GEN_BUS])
    
    # Get bus types and indices
    ref_buses = findall(bus[:, BUS_TYPE] .== REF)
    pv_buses = findall(bus[:, BUS_TYPE] .== PV)  
    pq_buses = findall(bus[:, BUS_TYPE] .== PQ)
    
    # Build network matrices (simplified for DC - only need connectivity)
    Ybus, _, _ = PowerFlow.makeYbus(baseMVA, bus, branch)
    
    # Build generator connection matrix (following run_dcopf.jl pattern)
    Cg = sparse(gbus, 1:ng, ones(ng), nb, ng)
    
    # Variable bounds using validated data
    var_bounds = create_dc_bounds(bus, gen, branch, baseMVA, on, options)
    
    # Cost coefficients  
    cost_coeffs = create_dc_costs(gen, on, baseMVA)
    
    # Create problem instance with correct dimensions
    return OPFProblem{DCFormulation, Float64}(
        DCFormulation(),
        nb, ng, nl,
        baseMVA,
        Ybus, spzeros(ComplexF64, 0, 0), spzeros(ComplexF64, 0, 0),
        Cg,
        Float64.(bus), Float64.(gen), Float64.(branch),
        ref_buses, pv_buses, pq_buses, on,  # Use 'on' for active_gens
        var_bounds, cost_coeffs, i2e
    )
end

"""
Create variable bounds for DC formulation using validated data.
"""
function create_dc_bounds(bus, gen, branch, baseMVA, on, options)
    nb, ng, nl = size(bus, 1), length(on), size(branch, 1)
    
    # Voltage angle bounds (in radians) - conservative defaults
    angle_limit_rad = deg2rad(options.angle_limit_deg)
    va_min = fill(-angle_limit_rad, nb)
    va_max = fill(angle_limit_rad, nb)
    
    # Reference bus angle is fixed at 0
    ref_buses = findall(bus[:, BUS_TYPE] .== REF)
    va_min[ref_buses] .= 0.0
    va_max[ref_buses] .= 0.0
    
    # Line flow bounds - use validated RATE_A data
    pij_max = branch[:, RATE_A] / baseMVA  # Already validated in bounds_validation
    pij_min = -pij_max  # symmetric limits
    
    # Generator bounds - use validated data
    pg_min = gen[on, PMIN] / baseMVA  # Already validated
    pg_max = gen[on, PMAX] / baseMVA  # Already validated
    
    # Final feasibility check
    total_load = sum(bus[:, PD]) / baseMVA
    total_gen_capacity = sum(pg_max)
    
    if total_load > total_gen_capacity * 0.98  # Strict check
        @warn "Tight generation capacity: load=$total_load, capacity=$total_gen_capacity"
    end
    
    return (
        va_min = va_min, va_max = va_max,
        pij_min = pij_min, pij_max = pij_max,
        pg_min = pg_min, pg_max = pg_max
    )
end

"""
Create cost coefficients for DC formulation (following run_dcopf.jl).
"""
function create_dc_costs(gen, on, baseMVA)
    ng = length(on)
    
    # Initialize cost vectors
    c_linear = zeros(ng)
    c_quadratic = zeros(ng)
    
    # Following the exact pattern from run_dcopf.jl
    for i in 1:ng
        gen_idx = on[i]
        if size(gen, 2) >= NCOST && gen[gen_idx, NCOST] > 0
            ncost = Int(gen[gen_idx, NCOST])
            
            if ncost >= 3 && size(gen, 2) >= COST+1
                # Quadratic cost: c2*Pg^2 + c1*Pg + c0
                c_quadratic[i] = gen[gen_idx, COST] * baseMVA^2
                c_linear[i] = gen[gen_idx, COST + 1] * baseMVA
            elseif ncost == 2 && size(gen, 2) >= COST
                # Linear cost: c1*Pg + c0
                c_linear[i] = gen[gen_idx, COST] * baseMVA
            else
                # Default linear cost
                c_linear[i] = 10.0
            end
        else
            # Default linear cost
            c_linear[i] = 10.0
        end
    end
    
    return (linear = c_linear, quadratic = c_quadratic)
end

# ============================================================================
# Constraint Functions
# ============================================================================

"""
Evaluate DC power balance constraints (following run_dcopf.jl).
"""
function evaluate_power_balance(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    indices = get_variable_indices(DCFormulation, problem.nb, problem.ng, problem.nl)
    
    va = x[indices.va]
    pij = x[indices.pij] 
    pg = x[indices.pg]
    
    # Load injection at each bus (following run_dcopf.jl)
    Pd = problem.bus_data[:, PD] / problem.baseMVA
    Gs = problem.bus_data[:, GS] / problem.baseMVA
    Pinj = Pd + Gs
    
    # Build Cft matrix exactly as in run_dcopf.jl
    f = Int.(problem.branch_data[:, F_BUS])
    t = Int.(problem.branch_data[:, T_BUS])
    i = [1:problem.nl; 1:problem.nl]
    Cft = sparse(i, [f; t], [ones(problem.nl); -ones(problem.nl)], problem.nl, problem.nb)
    
    # Power balance: Cft' * pij = Cg * pg - Pinj
    balance = -Cft' * pij + problem.Cg * pg - Pinj
    
    return balance
end

"""
Evaluate DC line flow definition constraints (following run_dcopf.jl).
"""
function evaluate_line_flow(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    indices = get_variable_indices(DCFormulation, problem.nb, problem.ng, problem.nl)
    
    va = x[indices.va]
    pij = x[indices.pij]
    
    # Line flows: pij = (va_from - va_to) / x_branch
    f = Int.(problem.branch_data[:, F_BUS])
    t = Int.(problem.branch_data[:, T_BUS])
    x_branch = problem.branch_data[:, BR_X]
    
    i = [1:problem.nl; 1:problem.nl]
    Cft = sparse(i, [f; t], [ones(problem.nl); -ones(problem.nl)], problem.nl, problem.nb)
    
    # Flow definition: pij = diag(1./branch(:,BR_X))*Cft*va
    flow_def = pij - Diagonal(1 ./ x_branch) * Cft * va
    
    return flow_def
end

"""
Evaluate DC inequality constraints (bounds).
"""
function evaluate_bounds(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    bounds = problem.var_bounds
    
    # Combined bounds: [x - xmax; -x + xmin] <= 0
    xmin = [bounds.va_min; bounds.pij_min; bounds.pg_min]
    xmax = [bounds.va_max; bounds.pij_max; bounds.pg_max]
    
    return [x - xmax; -x + xmin]
end

# ============================================================================
# Jacobian Functions (following run_dcopf.jl pattern)
# ============================================================================

"""
Compute Jacobian of DC power balance constraints.
"""
function power_balance_jacobian(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    # Build Cft matrix
    f = Int.(problem.branch_data[:, F_BUS])
    t = Int.(problem.branch_data[:, T_BUS])
    i = [1:problem.nl; 1:problem.nl]
    Cft = sparse(i, [f; t], [ones(problem.nl); -ones(problem.nl)], problem.nl, problem.nb)
    
    # Jacobian: ∂(power_balance)/∂x = [0 | -Cft' | Cg]
    J = [spzeros(problem.nb, problem.nb) -Cft' problem.Cg]
    
    return J'
end

"""
Compute Jacobian of DC line flow constraints.
"""
function line_flow_jacobian(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    # Build matrices
    f = Int.(problem.branch_data[:, F_BUS])
    t = Int.(problem.branch_data[:, T_BUS])
    x_branch = problem.branch_data[:, BR_X]
    
    i = [1:problem.nl; 1:problem.nl]
    Cft = sparse(i, [f; t], [ones(problem.nl); -ones(problem.nl)], problem.nl, problem.nb)
    
    # Jacobian: ∂(line_flow)/∂x = [-diag(1/x)*Cft | I | 0]
    J = [-Diagonal(1 ./ x_branch) * Cft I(problem.nl) spzeros(problem.nl, problem.ng)]
    
    return J'
end

"""
Compute bounds Jacobian (identity matrices).
"""
function bounds_jacobian(problem::OPFProblem{DCFormulation}, x::Vector{Float64})
    nx = length(x)
    return [I(nx); -I(nx)]'
end

# ============================================================================
# Solution Update Functions
# ============================================================================

"""
Update case data for DC solution (following run_dcopf.jl).
"""
function update_dc_case_data!(case_data, problem, solution)
    baseMVA = case_data["baseMVA"]
    
    # Convert back to external numbering
    bus_int = copy(problem.bus_data)
    gen_int = copy(problem.gen_data)  
    branch_int = copy(problem.branch_data)
    
    # Update voltage angles (convert from radians to degrees)
    bus_int[:, VA] = rad2deg.(solution.variables.va)
    
    # Update generator outputs
    gen_int[problem.active_gens, PG] = solution.variables.pg * baseMVA
    
    # Ensure branch matrix has enough columns for power flow results
    if size(branch_int, 2) < PT
        branch_new = zeros(size(branch_int, 1), max(PT, QT))
        branch_new[:, 1:size(branch_int, 2)] = branch_int
        branch_int = branch_new
    end
    
    # Update line flows
    branch_int[:, PF] = solution.variables.pij * baseMVA
    branch_int[:, PT] = -solution.variables.pij * baseMVA
    if size(branch_int, 2) >= QT
        branch_int[:, QF] .= 0.0
        branch_int[:, QT] .= 0.0
    end
    
    # Convert to external numbering
    bus_ext, gen_ext, branch_ext, _ = PowerFlow.int2ext(
        problem.int_to_ext, bus_int, gen_int, branch_int, zeros(0,8), zeros(0,8)
    )
    
    case_data["busAC"] = bus_ext
    case_data["genAC"] = gen_ext  
    case_data["branchAC"] = branch_ext
end

"""
Create intelligent initial point for DC optimization (following validated bounds).
"""
function create_dc_initial_point(problem::OPFProblem{DCFormulation}, options::OPFOptions)
    indices = get_variable_indices(DCFormulation, problem.nb, problem.ng, problem.nl)
    nx = indices.total
    
    # DC initialization following run_dcopf.jl pattern
    va0 = zeros(problem.nb)
    pij0 = zeros(problem.nl)
    
    # Smart generation initialization using validated bounds
    total_load = sum(problem.bus_data[:, PD]) / problem.baseMVA
    pg_capacity = problem.var_bounds.pg_max - problem.var_bounds.pg_min
    
    if sum(pg_capacity) > 0
        # Distribute load based on validated capacity
        pg0 = problem.var_bounds.pg_min + pg_capacity * (total_load / sum(pg_capacity))
    else
        pg0 = (problem.var_bounds.pg_min + problem.var_bounds.pg_max) / 2
    end
    
    # Ensure bounds satisfaction
    pg0 = clamp.(pg0, problem.var_bounds.pg_min, problem.var_bounds.pg_max)
    
    x0 = zeros(nx)
    x0[indices.va] = va0
    x0[indices.pij] = pij0  
    x0[indices.pg] = pg0
    
    # Use bounds validation for initial point
    bounds = problem.var_bounds
    xmin = [bounds.va_min; bounds.pij_min; bounds.pg_min]
    xmax = [bounds.va_max; bounds.pij_max; bounds.pg_max]
    
    # Create legacy options for validate_feasible_initialization
    legacy_opt = Dict{String, Any}(
        "BOUNDS_MARGIN" => options.bounds_margin,
        "INIT_PROJECTION_WARN" => 0.1
    )
    
    x0_validated = validate_feasible_initialization(x0, xmin, xmax, legacy_opt)
    
    return x0_validated
end
