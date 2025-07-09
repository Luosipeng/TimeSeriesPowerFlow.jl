# A standard form problem is defined as:
# max c^T x
# s.t. Ax sense b (where sense can be '=', '<', or '>')
#      lb <= x <= ub
#     x \in Z^p \times R^q
# Two types of variables are supported: continuous and integer

using Gurobi
using SparseArrays

# Main implementation with full parameter set
function mixed_integer_linear_programming(cobj::Vector, A::AbstractMatrix, b::Vector, sense::Vector, lb::Vector, ub::Vector, vtype::Vector, model_sense::String)
    # Convert to sparse matrix if not already sparse for memory efficiency
    A_sparse = issparse(A) ? A : sparse(A)
    
    # 0: initialize model with parameters settings
    env_p = Ref{Ptr{Cvoid}}()
    error = GRBloadenv(env_p, "")
    env = env_p[]

    GRBsetparam(env, "OutputFlag", "0"); # Update environment parameters
    model_p = Ref{Ptr{Cvoid}}()
    error = GRBnewmodel(env, model_p, "milp", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
    model = model_p[]

    # 1:  add variables
    NumVars = Ref{Cint}()
    NumConstrs = Ref{Cint}()
    (NumConstrs, NumVars) = size(A_sparse)
    
    # Convert variable types to Gurobi format
    vtype_converted = Vector{Int8}(undef, length(vtype))
    for i in eachindex(vtype)
        if vtype[i] == 'C' || vtype[i] == 'c'
            vtype_converted[i] = GRB_CONTINUOUS
        elseif vtype[i] == 'B' || vtype[i] == 'b'
            vtype_converted[i] = GRB_BINARY
        elseif vtype[i] == 'I' || vtype[i] == 'i'
            vtype_converted[i] = GRB_INTEGER
        else
            vtype_converted[i] = GRB_CONTINUOUS  # Default to continuous
        end
    end
    
    error = GRBaddvars(
        model, # model
        NumVars,      # : numvars
        0,      # : numnz
        C_NULL, # : *vbeg
        C_NULL, # : *vind
        C_NULL, # : *vval
        cobj,   # : *obj
        lb,    # : *lb
        ub,    # : *ub
        vtype_converted, # : *vtype
        C_NULL   # : **varnames
        )
    
    # 2:  add constraints - optimized for sparse matrices
    for i in 1 : NumConstrs
        # For sparse matrices, find non-zero elements in row i
        if issparse(A_sparse)
            row_indices = Int[]
            row_values = Float64[]
            
            # Iterate through the sparse matrix to find non-zeros in row i
            for j in 1:NumVars
                val = A_sparse[i, j]
                if val != 0
                    push!(row_indices, j - 1)  # Gurobi uses 0-based indexing
                    push!(row_values, val)
                end
            end
            
            numnz = length(row_indices)
            if numnz > 0
                nonzero_indices = convert(Vector{Cint}, row_indices)
                val = row_values
                
                # Determine constraint sense
                grb_sense = if sense[i] == '='
                    GRB_EQUAL
                elseif sense[i] == '<'
                    GRB_LESS_EQUAL
                elseif sense[i] == '>'
                    GRB_GREATER_EQUAL
                else
                    GRB_LESS_EQUAL  # Default
                end
                
                error = GRBaddconstr(
                    model,   # : *model
                    numnz,       # : numnz
                    nonzero_indices,   # : *cind
                    val,   # : *cval
                    grb_sense, # : sense
                    b[i],   # : rhs
                    C_NULL,    # : *constrname
                    )
            end
        else
            # Handle dense matrices
            nonzero_indices = findall(!iszero, A_sparse[i,:])
            numnz = length(nonzero_indices)
            if numnz > 0
                val = zeros(numnz)
                for j in 1:numnz
                    val[j] = A_sparse[i, nonzero_indices[j]]
                end
                for j in 1:numnz # Gurobi starts from 0 in the index
                    nonzero_indices[j] -= 1;
                end
                # convert nonzero_indices to Cint
                nonzero_indices = map(x -> Cint(x), nonzero_indices);
                
                # Determine constraint sense
                grb_sense = if sense[i] == '='
                    GRB_EQUAL
                elseif sense[i] == '<'
                    GRB_LESS_EQUAL
                elseif sense[i] == '>'
                    GRB_GREATER_EQUAL
                else
                    GRB_LESS_EQUAL  # Default
                end
                
                error = GRBaddconstr(
                    model,   # : *model
                    numnz,       # : numnz
                    nonzero_indices,   # : *cind
                    val,   # : *cval
                    grb_sense, # : sense
                    b[i],   # : rhs
                    C_NULL,    # : *constrname
                    )
            end
        end
    end
    
    # 3: update model parameters
    if cmp(model_sense, "min") == 0
        error = GRBsetintattr(model, "ModelSense", GRB_MINIMIZE)
    else
        error = GRBsetintattr(model, "ModelSense", GRB_MAXIMIZE)
    end
    error = GRBoptimize(model)

    # 4: obtain results
    optimstatus = Ref{Cint}()
    objval = Ref{Cdouble}()
    runtime = Ref{Cdouble}()
    mip_gap = Ref{Cdouble}()
    sol = ones(NumVars)

    error = GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, mip_gap);
    error = GRBgetdblattr(model, GRB_DBL_ATTR_RUNTIME, runtime);
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, optimstatus);
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, NumVars, sol);

    GRBfreemodel(model);
    GRBfreeenv(env);
    
    # Return solution using symbol indexing for consistency
    return Dict(:x => sol, :objval => objval[], :runtime => runtime[], :mip_gap => mip_gap[], :status => optimstatus[])
end

# Convenience method for problems with only inequality constraints (assumes all '<=' constraints)
function mixed_integer_linear_programming(cobj::Vector, A::AbstractMatrix, b::Vector, lb::Vector, ub::Vector, vtype::Vector, model_sense::String)
    sense = fill('<', length(b))  # Default to '<=' constraints
    return mixed_integer_linear_programming(cobj, A, b, sense, lb, ub, vtype, model_sense)
end

# Method for dictionary input (backwards compatibility)
function mixed_integer_linear_programming(problem::Dict)
    cobj = problem[:obj]
    A = problem[:A]
    b = problem[:rhs]
    lb = problem[:lb]
    ub = problem[:ub]
    vtype = problem[:vtype]
    sense = problem[:sense]
    model_sense = problem[:modelsense]
    
    return mixed_integer_linear_programming(cobj, A, b, sense, lb, ub, vtype, model_sense)
end

# Method for MILP struct input (using new data format)
function mixed_integer_linear_programming(milp::MILP)
    return mixed_integer_linear_programming(milp.c, milp.A, milp.b, milp.sense, milp.lb, milp.ub, milp.vtype, milp.model_sense)
end

# Method that returns MILP_Solution struct
function solve_milp(milp::MILP)
    result = mixed_integer_linear_programming(milp)
    return dict_to_milp_solution(result)
end

# # Type I: 
# lb = [0.0, 0.0, 0.0]
# ub = [1.0, 1.0, 1.0]
# cobj = [1.0, 1.0, 2.0]
# vtype = [GRB_BINARY, GRB_BINARY, GRB_BINARY]
# A = [1.0 2.0 3.0; -1.0 -1.0 0.0]
# b = [4.0; -1.0]
# sense = ['<'; '<']
# result =  @time mixed_integer_linear_programming(cobj, A, b, sense, lb, ub, vtype, "max")
# print(result)

# Type II:
problem = Dict(
    :A => [7 -2; 0 1; 2 -2],
    :rhs =>  [14; 3; 3],
    :sense => ['<' for _ in 1:3],
    :obj => [4.0, -1.0],
    :lb => zeros(2),
    :ub => fill(Inf, 2),
    :vtype => ['I', 'I'],
    :vtype => ['<', '<'],
    :modelsense => "max",
)

@time mixed_integer_linear_programming(problem)

# Type III: Sparse matrix test
# # Create a large sparse constraint matrix to test optimization
# n_vars = 100
# n_constrs = 50
# density = 0.1  # 10% of entries are non-zero

# # Generate random sparse matrix
# A_sparse = sprand(n_constrs, n_vars, density)
# b_sparse = rand(n_constrs) * 10
# sense_sparse = ['<' for _ in 1:n_constrs]
# cobj_sparse = rand(n_vars)
# lb_sparse = zeros(n_vars)
# ub_sparse = ones(n_vars)
# vtype_sparse = ['C' for _ in 1:n_vars]  # All continuous variables

# println("Testing with sparse matrix ($(n_constrs)x$(n_vars), $(round(100*density, digits=1))% density)")
# println("Non-zero elements: $(nnz(A_sparse))")

# result_sparse = @time mixed_integer_linear_programming(cobj_sparse, A_sparse, b_sparse, sense_sparse, lb_sparse, ub_sparse, vtype_sparse, "max")
# println("Sparse matrix test completed. Status: $(result_sparse[:status])")

# # Type I: 
# lb = [0.0, 0.0, 0.0]
# ub = [1.0, 1.0, 1.0]
# cobj = [1.0, 1.0, 2.0]
# vtype = [GRB_BINARY, GRB_BINARY, GRB_BINARY]
# A = [1.0 2.0 3.0; -1.0 -1.0 0.0]
# b = [4.0; -1.0]
# sense = ['<'; '<']
# result =  @time mixed_integer_linear_programming(cobj, A, b, sense, lb, ub, vtype, "max")
# print(result)

# Type II:
# problem = Dict(
#     :A => [7 -2; 0 1; 2 -2],
#     :rhs =>  [14; 3; 3],
#     :sense => ['<' for _ in 1:3],
#     :obj => [4.0, -1.0],
#     :lb => zeros(2),
#     :ub => fill(Inf, 2),
#     :vtype => ['I', 'I'],
#     :vtype => ['<', '<'],
#     :modelsense => "max",
# )

# @time mixed_integer_linear_programming(problem)

# Type III: Sparse matrix test
# # Create a large sparse constraint matrix to test optimization
# n_vars = 100
# n_constrs = 50
# density = 0.1  # 10% of entries are non-zero

# # Generate random sparse matrix
# A_sparse = sprand(n_constrs, n_vars, density)
# b_sparse = rand(n_constrs) * 10
# sense_sparse = ['<' for _ in 1:n_constrs]
# cobj_sparse = rand(n_vars)
# lb_sparse = zeros(n_vars)
# ub_sparse = ones(n_vars)
# vtype_sparse = ['C' for _ in 1:n_vars]  # All continuous variables

# println("Testing with sparse matrix ($(n_constrs)x$(n_vars), $(round(100*density, digits=1))% density)")
# println("Non-zero elements: $(nnz(A_sparse))")

# result_sparse = @time mixed_integer_linear_programming(cobj_sparse, A_sparse, b_sparse, sense_sparse, lb_sparse, ub_sparse, vtype_sparse, "max")
# println("Sparse matrix test completed. Status: $(result_sparse[:status])")
