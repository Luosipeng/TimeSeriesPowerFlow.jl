"""
Unified solver interface for OptimalPowerFlow module.
Provides consistent API across all formulations with automatic solver selection.
"""

# Load existing solver
include(joinpath(dirname(@__DIR__), "..", "solvers", "interior_point_method.jl"))

# ============================================================================
# Main Solve Function (Public API)
# ============================================================================

"""
    solve_opf(case_data, formulation, options=OPFOptions())

Solve optimal power flow problem with unified interface.

# Arguments
- `case_data::Dict`: Power system case data
- `formulation::Type{<:AbstractOPFFormulation}`: OPF formulation type
- `options::OPFOptions`: Solver options (optional)

# Returns
- `OPFResult`: Comprehensive result structure

# Examples
```julia
# DC OPF
result = solve_opf(case30(), DCFormulation)

# AC OPF with custom options
opts = OPFOptions(max_iterations=150, tolerance=1e-6)
result = solve_opf(case30(), ACFormulation, opts)
```
"""
function solve_opf(case_data::Dict{String, Any}, 
                   formulation::Type{F}, 
                   options::OPFOptions = OPFOptions()) where {F<:AbstractOPFFormulation}
    
    # Start timing
    solve_start_time = time()
    
    try
        # Validate options for formulation
        validated_options = validate_options(options, formulation)
        
        # Create problem
        problem = create_problem(formulation, case_data, validated_options)
        
        # Solve problem
        solution = solve_problem(problem, validated_options)
        
        # Extract and update case data
        updated_case = extract_solution(problem, solution, case_data)
        
        # Calculate performance metrics
        solve_time = time() - solve_start_time
        memory_usage = @allocated solve_problem(problem, validated_options)
        
        # Create result
        return OPFResult{F}(
            formulation(),
            solution.success,
            solution.objective,
            solution.iterations,
            solve_time,
            updated_case,
            solution,
            memory_usage,
            solution.max_violation
        )
        
    catch e
        solve_time = time() - solve_start_time
        
        if options.verbose > 0
            println("OPF solve failed: $(e)")
        end
        
        # Return failure result
        return OPFResult{F}(
            formulation(),
            false,
            Inf,
            0,
            solve_time,
            case_data,
            OPFSolution{F, Float64}(
                formulation(), false, Inf, 0, solve_time,
                NamedTuple(), Inf, NamedTuple(), Dict{String, Any}()
            ),
            0,
            Inf
        )
    end
end

# ============================================================================
# Problem Solving
# ============================================================================

"""
Solve OPF problem using appropriate solver.
"""
function solve_problem(problem::OPFProblem{F, T}, options::OPFOptions) where {F, T}
    # Build optimization functions
    obj_func, grad_func, hess_func = build_objective(problem)
    eq_func, eq_jac_func = build_equality_constraints(problem)
    ineq_func, ineq_jac_func = build_inequality_constraints(problem)
    
    # Create Lagrangian functions
    lagrangian_grad, lagrangian_hess = build_lagrangian(
        grad_func, hess_func, eq_jac_func, ineq_jac_func, problem
    )
    
    # Initial point
    x0 = create_initial_point(problem, options)
    
    # Create nonlinear optimization structure
    nonlinear = NonConvexOPT(
        obj_func,        # objective function
        grad_func,       # gradient of objective
        hess_func,       # Hessian of objective
        eq_func,         # equality constraints g(x) = 0
        eq_jac_func,     # gradient of equality constraints
        ineq_func,       # inequality constraints h(x) <= 0
        ineq_jac_func,   # gradient of inequality constraints
        lagrangian_grad, # gradient of Lagrangian
        lagrangian_hess, # Hessian of Lagrangian
        x0               # initial point
    )
    
    # Configure solver
    solver_opts = get_solver_options(options, F)
    ipm = IPM(
        solver_opts["tolerance"],
        solver_opts["max_iterations"],
        get(solver_opts, "positive_tol", 1e-8),
        get(solver_opts, "feasibility_tolerance", solver_opts["tolerance"]),
        false,  # initial_point_projection
        get(solver_opts, "barrier_parameter", 0.995)
    )
    
    # Solve
    result = interior_point_method(nonlinear, ipm)
    
    # Process solution
    return process_solution(result, problem, options)
end

# ============================================================================
# Objective Function Building
# ============================================================================

"""
Build objective function for given formulation.
"""
function build_objective(problem::OPFProblem{F, T}) where {F, T}
    costs = problem.cost_coeffs
    indices = get_variable_indices(F, problem.nb, problem.ng, problem.nl)
    
    # Objective function
    obj_func = function(x)
        if F == DCFormulation
            pg = x[indices.pg]
            return sum(costs.quadratic .* pg.^2) + sum(costs.linear .* pg)
        elseif F == ACFormulation
            pg = x[indices.pg]
            return sum(costs.quadratic .* pg.^2) + sum(costs.linear .* pg)
        end
    end
    
    # Gradient
    grad_func = function(x)
        grad = zeros(length(x))
        if F == DCFormulation
            pg = x[indices.pg]
            grad[indices.pg] = 2 * costs.quadratic .* pg + costs.linear
        elseif F == ACFormulation
            pg = x[indices.pg]
            grad[indices.pg] = 2 * costs.quadratic .* pg + costs.linear
        end
        return grad
    end
    
    # Hessian
    hess_func = function(x)
        nx = length(x)
        H = spzeros(nx, nx)
        if !isempty(costs.quadratic)
            H[indices.pg, indices.pg] = Diagonal(2 * costs.quadratic)
        end
        return H
    end
    
    return obj_func, grad_func, hess_func
end

# ============================================================================
# Constraint Building
# ============================================================================

"""
Build equality constraints for given formulation.
"""
function build_equality_constraints(problem::OPFProblem{F, T}) where {F, T}
    if F == DCFormulation
        # DC: power balance + line flow definition
        eq_func = function(x)
            balance = evaluate_power_balance(problem, x)
            flows = evaluate_line_flow(problem, x)
            return [balance; flows]
        end
        
        eq_jac_func = function(x)
            jac_balance = power_balance_jacobian(problem, x)
            jac_flows = line_flow_jacobian(problem, x)
            return [jac_balance; jac_flows]  # Fixed: should be vertical concatenation
        end
        
    elseif F == ACFormulation
        # AC: ONLY power balance constraints (EXACTLY as run_acopf.jl)
        eq_func = function(x)
            return evaluate_ac_power_balance(problem, x)
        end
        
        eq_jac_func = function(x)
            return ac_power_balance_jacobian(problem, x)
        end
    end
    
    return eq_func, eq_jac_func
end

"""
Build inequality constraints (EXACTLY following run_acopf.jl h_ineq function).
"""
function build_inequality_constraints(problem::OPFProblem{F, T}) where {F, T}
    ineq_func = function(x)
        return evaluate_bounds(problem, x)  # This includes both bounds AND line flow for AC
    end
    
    ineq_jac_func = function(x)
        return bounds_jacobian(problem, x)  # This includes both bounds AND line flow Jacobian for AC
    end
    
    return ineq_func, ineq_jac_func
end

# ============================================================================
# Solution Processing
# ============================================================================

"""
Process solver result into OPFSolution.
"""
function process_solution(result, problem::OPFProblem{F, T}, options::OPFOptions) where {F, T}
    success = result.eflag == 1
    x_sol = result.x
    
    if success
        # Extract variables
        indices = get_variable_indices(F, problem.nb, problem.ng, problem.nl)
        variables = extract_variables(F, x_sol, indices)
        
        # Compute constraint violations
        eq_func, _ = build_equality_constraints(problem)
        eq_viol = maximum(abs.(eq_func(x_sol)))
        
        ineq_func, _ = build_inequality_constraints(problem)
        ineq_viol = maximum(max.(ineq_func(x_sol), 0.0))
        
        max_violation = max(eq_viol, ineq_viol)
        
        return OPFSolution{F, T}(
            problem.formulation,
            success,
            T(result.obj),
            result.iterations,
            0.0,  # Will be filled by caller
            variables,
            T(max_violation),
            (equality = T(eq_viol), inequality = T(ineq_viol)),
            Dict("solver_result" => result)
        )
    else
        return OPFSolution{F, T}(
            problem.formulation,
            false,
            T(Inf),
            result.iterations,
            0.0,
            NamedTuple(),
            T(Inf),
            NamedTuple(),
            Dict("solver_result" => result)
        )
    end
end

"""
Extract solution variables based on formulation.
"""
function extract_variables(::Type{DCFormulation}, x::Vector, indices)
    return (
        va = x[indices.va],
        pij = x[indices.pij],
        pg = x[indices.pg]
    )
end

function extract_variables(::Type{ACFormulation}, x::Vector, indices)
    return (
        va = x[indices.va],
        vm = x[indices.vm],
        pg = x[indices.pg],
        qg = x[indices.qg]
    )
end

# ============================================================================
# Solution Extraction
# ============================================================================

"""
Extract solution and update case data.
"""
function extract_solution(problem::OPFProblem{F, T}, solution::OPFSolution{F, T}, 
                         original_case::Dict{String, Any}) where {F, T}
    
    if !solution.success
        updated_case = deepcopy(original_case)
        updated_case["success"] = false
        updated_case["f"] = Inf
        return updated_case
    end
    
    # Update case data with solution
    updated_case = deepcopy(original_case)
    
    if F == DCFormulation
        update_dc_case_data!(updated_case, problem, solution)
    elseif F == ACFormulation
        update_ac_case_data!(updated_case, problem, solution)
    end
    
    updated_case["success"] = true
    updated_case["f"] = solution.objective
    updated_case["iterations"] = solution.iterations
    
    return updated_case
end

"""
Update case data for DC solution.
"""
function update_dc_case_data!(case_data, problem, solution)
    baseMVA = case_data["baseMVA"]
    
    # Convert back to external numbering
    bus_int = copy(problem.bus_data)
    gen_int = copy(problem.gen_data)
    branch_int = copy(problem.branch_data)
    
    # Update voltage angles
    bus_int[:, 9] = rad2deg.(solution.variables.va)  # VA column
    
    # Update generator outputs
    gen_int[problem.active_gens, 2] = solution.variables.pg * baseMVA  # PG
    
    # Update line flows
    if size(branch_int, 2) >= 14
        branch_int[:, 14] = solution.variables.pij * baseMVA  # PF
        branch_int[:, 16] = -solution.variables.pij * baseMVA  # PT
    end
    
    # Convert to external numbering
    bus_ext, gen_ext, branch_ext, _ = PowerFlow.int2ext(
        problem.int_to_ext, bus_int, gen_int, branch_int, zeros(0,8), zeros(0,8)
    )
    
    case_data["busAC"] = bus_ext
    case_data["genAC"] = gen_ext
    case_data["branchAC"] = branch_ext
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
Create intelligent initial point for optimization.
"""
function create_initial_point(problem::OPFProblem{F, T}, options::OPFOptions) where {F, T}
    if F == DCFormulation
        return create_dc_initial_point(problem, options)
    elseif F == ACFormulation
        return create_ac_initial_point(problem, options)
    else
        error("Unknown formulation type: $F")
    end
end

"""
Build Lagrangian gradient and Hessian functions.
"""
function build_lagrangian(grad_func, hess_func, eq_jac_func, ineq_jac_func, problem)
    lagrangian_grad = function(x, λ, μ)
        grad_f = grad_func(x)
        grad_g = eq_jac_func(x)
        grad_h = ineq_jac_func(x)
        
        return grad_f + grad_g * λ + grad_h * μ
    end
    
    lagrangian_hess = function(x, λ, μ)
        H_f = hess_func(x)
        # For most OPF problems, constraint Hessians are zero or simple
        # This can be extended for more complex constraint Hessians
        return H_f
    end
    
    return lagrangian_grad, lagrangian_hess
end

# ============================================================================
# AC-Specific Functions (Placeholder for full AC implementation)
# ============================================================================

"""
Evaluate AC power balance constraints.
"""
function evaluate_ac_power_balance(problem::OPFProblem{ACFormulation}, x::Vector{Float64})
    indices = get_variable_indices(ACFormulation, problem.nb, problem.ng, problem.nl)
    
    va = x[indices.va]
    vm = x[indices.vm]
    pg = x[indices.pg]
    qg = x[indices.qg]
    
    # Complex voltage
    V = vm .* exp.(1im * va)
    
    # Power injections
    S_calc = V .* conj.(problem.Ybus * V)
    P_calc = real(S_calc)
    Q_calc = imag(S_calc)
    
    # Load injections
    Pd = problem.bus_data[:, 3] / problem.baseMVA  # PD
    Qd = problem.bus_data[:, 4] / problem.baseMVA  # QD
    
    # Generator injections
    P_gen = problem.Cg * pg
    Q_gen = problem.Cg * qg
    
    # Power balance equations
    P_balance = P_calc - P_gen + Pd
    Q_balance = Q_calc - Q_gen + Qd
    
    return [P_balance; Q_balance]
end

# ============================================================================
# Legacy Compatibility Functions
# ============================================================================

"""
Legacy DC OPF interface for backward compatibility.
"""
function rundcopf(jpc::Dict{String, Any}, opt::Dict{String, Any})
    # Convert legacy options
    options = from_legacy_options(opt)
    
    # Solve using unified interface
    result = solve_opf(jpc, DCFormulation, options)
    
    # Return in legacy format
    legacy_result = result.case_data
    legacy_result["success"] = result.success
    legacy_result["f"] = result.objective
    legacy_result["iterations"] = result.iterations
    
    return legacy_result
end

"""
Legacy AC OPF interface for backward compatibility.
"""
function runacopf(jpc::Dict{String, Any}, opt::Dict{String, Any})
    # Convert legacy options
    options = from_legacy_options(opt)
    
    # Solve using unified interface
    result = solve_opf(jpc, ACFormulation, options)
    
    # Return in legacy format
    legacy_result = result.case_data
    legacy_result["success"] = result.success
    legacy_result["f"] = result.objective
    legacy_result["iterations"] = result.iterations
    
    return legacy_result
end