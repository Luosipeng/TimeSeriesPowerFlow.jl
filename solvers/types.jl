# define a struct for the problem
struct LP
    c::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
end

struct NonConvexOPT
    f::Function
    ∇f::Function
    ∇2f::Function
    g::Function # g(x) = 0
    ∇g::Function
    h::Function # h(x) \leq 0
    ∇h::Function
    Lx::Function
    Lxx::Function
    x0::Vector{Float64} # initial point
end

# define a struct for the algorithm parameters
struct IPM
    tol::Float64
    max_iter::Int
    positive_tol::Float64
    feasible_tol::Float64
    initial_point_projection::Bool
    ξ::Float64
end

# define a struct for the solution
struct Solution
    x::Vector{Float64}
    λ::Vector{Float64}
    s::Vector{Float64}
    obj::Float64
end

# define the struct of the history
mutable struct History
    x_record::Matrix{Float64}
    λ_record::Matrix{Float64}
    μ_record::Matrix{Float64}
    z_record::Matrix{Float64}
    obj_record::Vector{Float64}
    feascond_record::Vector{Float64}
    gradcond_record::Vector{Float64}
    compcond_record::Vector{Float64}
    costcond_record::Vector{Float64}
end

# define the structure for the solution of the interior point method, x, λ, μ, obj, eflag, x_record, λ_record, μ_record, z_record, obj_record, feascond_record, gradcond_record, compcond_record
struct IPM_Solution
    x::Vector{Float64}
    λ::Vector{Float64}
    μ::Vector{Float64}
    obj::Float64
    eflag::Int
    iterations::Int  # Add iterations field
    hist::History
end

# define a struct for mixed integer linear programming problems
# Standard form: optimize c^T x subject to A*x sense b, lb <= x <= ub, x[i] in {continuous, integer, binary}
struct MILP
    c::Vector{Float64}                    # objective coefficient vector
    A::AbstractMatrix{Float64}            # constraint matrix (can be sparse or dense)
    b::Vector{Float64}                    # right-hand side vector
    sense::Vector{Char}                   # constraint sense ('=', '<', '>')
    lb::Vector{Float64}                   # lower bounds on variables
    ub::Vector{Float64}                   # upper bounds on variables
    vtype::Vector{Char}                   # variable types ('C'=continuous, 'I'=integer, 'B'=binary)
    model_sense::String                   # optimization direction ("min" or "max")
end

# Constructor for MILP with default values
function MILP(c::Vector{Float64}, A::AbstractMatrix{Float64}, b::Vector{Float64};
              sense::Vector{Char} = fill('<', length(b)),
              lb::Vector{Float64} = zeros(length(c)),
              ub::Vector{Float64} = fill(Inf, length(c)),
              vtype::Vector{Char} = fill('C', length(c)),
              model_sense::String = "min")
    return MILP(c, A, b, sense, lb, ub, vtype, model_sense)
end

# define a struct for MILP solution results
struct MILP_Solution
    x::Vector{Float64}                    # optimal solution vector
    objval::Float64                       # optimal objective value
    status::Int                           # solution status (Gurobi status codes)
    runtime::Float64                      # solution time in seconds
    mip_gap::Float64                      # MIP gap (for integer problems)
    iterations::Int                       # number of iterations (if available)
end

# Constructor for MILP_Solution with default values
function MILP_Solution(x::Vector{Float64}, objval::Float64, status::Int;
                      runtime::Float64 = 0.0,
                      mip_gap::Float64 = 0.0,
                      iterations::Int = 0)
    return MILP_Solution(x, objval, status, runtime, mip_gap, iterations)
end

# Utility function to convert MILP struct to dictionary format for compatibility
function milp_to_dict(milp::MILP)
    return Dict(
        :obj => milp.c,
        :A => milp.A,
        :rhs => milp.b,
        :sense => milp.sense,
        :lb => milp.lb,
        :ub => milp.ub,
        :vtype => milp.vtype,
        :modelsense => milp.model_sense
    )
end

# Utility function to convert dictionary result back to MILP_Solution struct
function dict_to_milp_solution(result_dict::Dict)
    return MILP_Solution(
        result_dict[:x],
        result_dict[:objval],
        result_dict[:status],
        runtime = get(result_dict, :runtime, 0.0),
        mip_gap = get(result_dict, :mip_gap, 0.0),
        iterations = get(result_dict, :iterations, 0)
    )
end

# SQP solver parameters
struct SQP
    feasible_tol::Float64    # Feasibility tolerance
    max_iter::Int            # Maximum iterations  
    min_step_size::Float64   # Minimum step size
    penalty_param::Float64   # Penalty parameter for merit function
end

# MIPS solution structure
struct MIPS_Solution
    x::Vector{Float64}       # Optimal solution
    λ::Vector{Float64}       # Equality constraint multipliers
    μ::Vector{Float64}       # Inequality constraint multipliers
    obj::Float64             # Objective value
    eflag::Int              # Exit flag
    iterations::Int         # Number of iterations
    output::Dict            # Additional output information
end

# SQP solution structure  
struct SQP_Solution
    x::Vector{Float64}       # Optimal solution
    λ::Vector{Float64}       # Equality constraint multipliers
    μ::Vector{Float64}       # Inequality constraint multipliers
    obj::Float64             # Objective value
    eflag::Int              # Exit flag
    iterations::Int         # Number of iterations
    hist::History           # Iteration history
end