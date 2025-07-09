using SparseArrays
using HiGHS
"""
    Standard mixed-integer linear programming (MILP) solver using HiGHS.
    The programming is to solve the following type of MILP problems:
    ⩟ min/max c'x
    ⩟ s.t. Ax ≤ b (inequality constraints)
    ⩟ Aeqx = beq (equality constraints)
    ⩟ lb ≤ x ≤ ub (bounds)
    ⩟ x ∈ {C, I, B} (continuous, integer, binary)
"""
# HiGHS status codes - defined at module level
const kHighsModelStatusNotset = 0
const kHighsModelStatusLoadError = 1
const kHighsModelStatusModelError = 2
const kHighsModelStatusPresolveError = 3
const kHighsModelStatusSolveError = 4
const kHighsModelStatusPostsolveError = 5
const kHighsModelStatusEmpty = 6
const kHighsModelStatusOptimal = 7
const kHighsModelStatusInfeasible = 8
const kHighsModelStatusPrimalInfeasible = 8
const kHighsModelStatusPrimalDualInfeasible = 8
const kHighsModelStatusDualInfeasible = 9
const kHighsModelStatusPrimalDualInfeasibleOrUnbounded = 9
const kHighsModelStatusUnbounded = 10
const kHighsModelStatusObjectiveBound = 11
const kHighsModelStatusObjectiveTarget = 12
const kHighsModelStatusTimeLimit = 13
const kHighsModelStatusIterationLimit = 14
const kHighsModelStatusUnknown = 15

function solve_milp_sparse(
    c::Vector{Float64},
    A::SparseMatrixCSC{Float64},
    b::Vector{Float64},
    Aeq::SparseMatrixCSC{Float64},
    beq::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64},
    vtype::Vector{Char};
    sense::String="min",
    verbose::Bool=false,
    debug::Bool=false
)
    # Validate input dimensions
    n = length(c)
    if size(A, 2) != n && !isempty(A)
        error("Inequality constraint matrix A has $(size(A, 2)) columns, but expected $n")
    end
    if size(Aeq, 2) != n && !isempty(Aeq)
        error("Equality constraint matrix Aeq has $(size(Aeq, 2)) columns, but expected $n")
    end
    if length(lb) != n
        error("Lower bounds vector has $(length(lb)) elements, but expected $n")
    end
    if length(ub) != n
        error("Upper bounds vector has $(length(ub)) elements, but expected $n")
    end
    if length(vtype) != n
        error("Variable type vector has $(length(vtype)) elements, but expected $n")
    end
    
    # Create a new HiGHS model
    model = Highs_create()
    
    # Set logging options
    if !verbose
        Highs_setBoolOptionValue(model, "output_flag", false)
        Highs_setBoolOptionValue(model, "log_to_console", false)
    else
        Highs_setBoolOptionValue(model, "output_flag", true)
        Highs_setBoolOptionValue(model, "log_to_console", true)
    end
    
    # Set solver precision parameters - use more relaxed tolerances
    Highs_setDoubleOptionValue(model, "primal_feasibility_tolerance", 1e-6)
    Highs_setDoubleOptionValue(model, "dual_feasibility_tolerance", 1e-6)
    
    # Set optimization direction
    if lowercase(sense) == "min"
        Highs_changeObjectiveSense(model, 1)  # 1 = minimize
        if verbose
            println("Set to MINIMIZE")
        end
    elseif lowercase(sense) == "max"
        Highs_changeObjectiveSense(model, -1)  # -1 = maximize
        if verbose
            println("Set to MAXIMIZE")
        end
    else
        error("Invalid optimization sense: $sense. Use 'min' or 'max'")
    end
    
    # Add variables with objective coefficients and bounds
    if verbose
        println("Variables:")
        println("  Objective coefficients: $c")
        println("  Lower bounds: $lb") 
        println("  Upper bounds: $ub")
        println("  Variable types: $vtype")
    end
    
    # Add variables one by one to ensure correct parameter order
    for i in 1:n
        # Add one variable at a time: cost, lower_bound, upper_bound
        status = Highs_addCol(model, c[i], lb[i], ub[i], 0, C_NULL, C_NULL)
        if verbose && i <= 3
            println("Added variable $i: cost=$(c[i]), lb=$(lb[i]), ub=$(ub[i]), status=$status")
        end
        
        # Set variable type
        if vtype[i] == 'I'
            status = Highs_changeColIntegrality(model, i-1, 1)  # Integer
            if verbose && i <= 3
                println("Set variable $(i-1) to integer, status: $status")
            end
        elseif vtype[i] == 'B'
            status = Highs_changeColIntegrality(model, i-1, 1)  # Binary
            if verbose && i <= 3
                println("Set variable $(i-1) to binary, status: $status")
            end
            # For binary variables, ensure bounds are [0,1]
            status = Highs_changeColBounds(model, i-1, 0.0, 1.0)
            if verbose && i <= 3
                println("Set binary bounds for variable $(i-1), status: $status")
            end
        end
    end
    
    # Add inequality constraints (Ax ≤ b)
    if !isempty(A) && !isempty(b)
        m_ineq = size(A, 1)
        if verbose
            println("Adding $m_ineq inequality constraints")
        end
        
        for i in 1:m_ineq
            # Extract non-zero elements from this row more efficiently
            row_vals = Float64[]
            row_inds = Int32[]
            
            # Use sparse matrix structure for better performance
            for j in 1:n
                val = A[i, j]
                if abs(val) > 1e-12  # Use smaller threshold
                    push!(row_vals, val)
                    push!(row_inds, j-1)  # 0-based indexing for HiGHS
                end
            end
            
            if !isempty(row_inds)
                # Add constraint: -∞ ≤ a'x ≤ b
                Highs_addRow(model, -1e20, b[i], length(row_inds), row_inds, row_vals)  # Use finite lower bound instead of -Inf
                
                if verbose && i <= 5  # Print first few constraints for debugging
                    println("Inequality constraint $i: indices=$(row_inds.+1), values=$row_vals, rhs=$(b[i])")
                end
            end
        end
    end
    
    # Add equality constraints (Aeq*x = beq)
    if !isempty(Aeq) && !isempty(beq)
        m_eq = size(Aeq, 1)
        if verbose
            println("Adding $m_eq equality constraints")
        end
        
        for i in 1:m_eq
            # Extract non-zero elements from this row
            row_vals = Float64[]
            row_inds = Int32[]
            
            for j in 1:n
                val = Aeq[i, j]
                if abs(val) > 1e-12  # Use smaller threshold
                    push!(row_vals, val)
                    push!(row_inds, j-1)  # 0-based indexing for HiGHS
                end
            end
            
            if !isempty(row_inds)
                # Add constraint: beq ≤ a'x ≤ beq
                Highs_addRow(model, beq[i], beq[i], length(row_inds), row_inds, row_vals)
                
                if verbose && i <= 5  # Print first few constraints for debugging
                    println("Equality constraint $i: indices=$(row_inds.+1), values=$row_vals, rhs=$(beq[i])")
                end
            end
        end
    end
    
    # Start timer for runtime calculation
    start_time = time()
    
    # Solve the model
    run_status = Highs_run(model)
    
    # Calculate runtime
    runtime = time() - start_time
    
    # Get model status
    model_status = Highs_getModelStatus(model)
    
    # Map HiGHS status to a string
    status_map = Dict(
        kHighsModelStatusNotset => "Not Set",
        kHighsModelStatusLoadError => "Load error",
        kHighsModelStatusModelError => "Model error",
        kHighsModelStatusPresolveError => "Presolve error",
        kHighsModelStatusSolveError => "Solve error",
        kHighsModelStatusPostsolveError => "Postsolve error",
        kHighsModelStatusEmpty => "Empty",
        kHighsModelStatusOptimal => "Optimal solution found",
        kHighsModelStatusInfeasible => "Problem is infeasible",
        kHighsModelStatusPrimalDualInfeasibleOrUnbounded => "Problem is infeasible or unbounded",
        kHighsModelStatusUnbounded => "Problem is unbounded",
        kHighsModelStatusObjectiveBound => "Objective bound reached",
        kHighsModelStatusObjectiveTarget => "Objective target reached",
        kHighsModelStatusTimeLimit => "Time limit reached",
        kHighsModelStatusIterationLimit => "Iteration limit reached",
        kHighsModelStatusUnknown => "Unknown"
    )
    
    status_name = get(status_map, model_status, "Unknown status ($model_status)")
    
    # Prepare solution array
    sol = zeros(n)
    objval = 0.0
    mip_gap = 0.0
    
    # Check if optimal solution was found
    if model_status == kHighsModelStatusOptimal
        # Get objective value
        objval = Highs_getObjectiveValue(model)
        if verbose
            println("Retrieved objective value: $objval")
        end
        
        # Get solution values - try different approach
        try
            col_value, col_dual, row_value, row_dual = Highs_getSolution(model)
            if verbose
                println("Solution array length: $(length(col_value))")
                println("Raw solution: $col_value")
            end
            
            if length(col_value) == n
                sol = col_value
            else
                if verbose
                    println("Warning: Solution array length mismatch, trying alternative method")
                end
                # Alternative: get solution one variable at a time
                sol = zeros(n)
                for i in 1:n
                    sol[i] = Highs_getColValue(model, i-1)
                end
            end
        catch e
            if verbose
                println("Error getting solution: $e")
                println("Trying alternative solution retrieval...")
            end
            # Fallback: get solution one variable at a time
            sol = zeros(n)
            for i in 1:n
                try
                    sol[i] = Highs_getColValue(model, i-1)
                catch
                    sol[i] = 0.0
                end
            end
        end
        
        if verbose
            println("Final solution: $sol")
            println("Solution length: $(length(sol))")
        end
        
        # Get MIP gap if it's a MIP problem - handle missing function gracefully
        has_discrete = any(vtype .!= 'C')
        if has_discrete
            try
                mip_gap = Highs_getMipGap(model)
            catch
                # MIP gap function not available in this version of HiGHS
                mip_gap = 0.0
                if verbose
                    println("Note: MIP gap not available in this HiGHS version")
                end
            end
        end
        
        # Verify constraints if debug is enabled
        if debug
            println("\nSolution verification:")
            
            # Check inequality constraints
            if !isempty(A) && !isempty(b)
                max_violation = 0.0
                violation_count = 0
                
                for i in axes(A, 1)
                    lhs = 0.0
                    # Manually compute dot product to avoid potential sparse matrix issues
                    for j in 1:n
                        lhs += A[i, j] * sol[j]
                    end
                    
                    violation = max(0.0, lhs - b[i])
                    max_violation = max(max_violation, violation)
                    if violation > 1e-6
                        violation_count += 1
                        if violation_count <= 10  # Limit output to first 10 violations
                            println("Inequality constraint $i violated: $lhs > $(b[i]), violation = $violation")
                        end
                    end
                end
                
                if violation_count > 10
                    println("... and $(violation_count - 10) more inequality violations")
                end
                
                println("Max inequality violation: $max_violation (total violations: $violation_count)")
            end
            
            # Check equality constraints
            if !isempty(Aeq) && !isempty(beq)
                max_violation = 0.0
                violation_count = 0
                
                for i in axes(Aeq, 1)
                    lhs = 0.0
                    # Manually compute dot product to avoid potential sparse matrix issues
                    for j in 1:n
                        lhs += Aeq[i, j] * sol[j]
                    end
                    
                    violation = abs(lhs - beq[i])
                    max_violation = max(max_violation, violation)
                    if violation > 1e-6
                        violation_count += 1
                        if violation_count <= 10  # Limit output to first 10 violations
                            println("Equality constraint $i violated: $lhs != $(beq[i]), violation = $violation")
                        end
                    end
                end
                
                if violation_count > 10
                    println("... and $(violation_count - 10) more equality violations")
                end
                
                println("Max equality violation: $max_violation (total violations: $violation_count)")
            end
        end
    else
        objval = NaN
        mip_gap = NaN
        
        # Print a message about the non-optimal status
        if verbose
            println("No optimal solution found. Status: $(status_name)")
            println("Model status code: $model_status")
        end
    end
    
    # Free resources
    Highs_destroy(model)
    
    # Return results
    return Dict(
        :x => sol, 
        :objval => objval, 
        :runtime => runtime, 
        :mip_gap => mip_gap,
        :status => model_status,
        :status_name => status_name
    )
end

# Define the HiGHS wrapper functions
function Highs_create()
    ccall((:Highs_create, "libhighs"), Ptr{Cvoid}, ())
end

function Highs_destroy(highs)
    ccall((:Highs_destroy, "libhighs"), Cvoid, (Ptr{Cvoid},), highs)
end

function Highs_run(highs)
    ccall((:Highs_run, "libhighs"), Cint, (Ptr{Cvoid},), highs)
end

function Highs_setBoolOptionValue(highs, option, value::Bool)
    ccall((:Highs_setBoolOptionValue, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cstring, Cint), highs, option, value ? 1 : 0)
end

function Highs_setDoubleOptionValue(highs, option, value)
    ccall((:Highs_setDoubleOptionValue, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cstring, Cdouble), highs, option, value)
end

function Highs_changeObjectiveSense(highs, sense)
    ccall((:Highs_changeObjectiveSense, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cint), highs, sense)
end

function Highs_addVars(highs, num_vars, costs, lower_bounds, upper_bounds)
    ccall((:Highs_addVars, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
          highs, num_vars, costs, lower_bounds, upper_bounds)
end

function Highs_changeColIntegrality(highs, col, integrality)
    ccall((:Highs_changeColIntegrality, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cint, Cint), highs, col, integrality)
end

function Highs_changeColBounds(highs, col, lower, upper)
    ccall((:Highs_changeColBounds, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, col, lower, upper)
end

function Highs_addRow(highs, lower, upper, num_entries, indices, values)
    ccall((:Highs_addRow, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), 
          highs, lower, upper, num_entries, indices, values)
end

function Highs_addCol(highs, cost, lower, upper, num_entries, indices, values)
    ccall((:Highs_addCol, "libhighs"), Cint, 
          (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), 
          highs, cost, lower, upper, num_entries, indices, values)
end

function Highs_getModelStatus(highs)
    ccall((:Highs_getModelStatus, "libhighs"), Cint, (Ptr{Cvoid},), highs)
end

function Highs_getObjectiveValue(highs)
    value_ref = Ref{Cdouble}(0.0)
    ccall((:Highs_getObjectiveValue, "libhighs"), Cint, 
          (Ptr{Cvoid}, Ptr{Cdouble}), highs, value_ref)
    return value_ref[]
end

function Highs_getSolution(highs)
    n_col = Ref{Cint}(0)
    n_row = Ref{Cint}(0)
    
    # Get model dimensions
    ccall((:Highs_getNumCol, "libhighs"), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, n_col)
    ccall((:Highs_getNumRow, "libhighs"), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, n_row)
    
    # Allocate memory for solution vectors
    col_value = zeros(Cdouble, n_col[])
    col_dual = zeros(Cdouble, n_col[])
    row_value = zeros(Cdouble, n_row[])
    row_dual = zeros(Cdouble, n_row[])
    
    # Get solution
    ccall((:Highs_getSolution, "libhighs"), Cint, 
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
          highs, col_value, col_dual, row_value, row_dual)
    
    return col_value, col_dual, row_value, row_dual
end

function Highs_getMipGap(highs)
    # Try to get MIP gap, but handle gracefully if function doesn't exist
    try
        gap_ref = Ref{Cdouble}(0.0)
        ccall((:Highs_getMipGap, "libhighs"), Cint, (Ptr{Cvoid}, Ptr{Cdouble}), highs, gap_ref)
        return gap_ref[]
    catch
        # Function not available in this version
        return 0.0
    end
end

function Highs_getColValue(highs, col)
    value_ref = Ref{Cdouble}(0.0)
    status = ccall((:Highs_getColValue, "libhighs"), Cint, 
                   (Ptr{Cvoid}, Cint, Ptr{Cdouble}), highs, col, value_ref)
    if status != 0
        error("Failed to get column value for column $col")
    end
    return value_ref[]
end

function test_mixed_integer()
    println("=== Test Case: Mixed Integer Problem with HiGHS ===")
    # max 2x1 + 3x2
    # s.t. x1 + x2 <= 3.5
    #      2x1 + x2 <= 5
    #      x1 integer, x2 continuous
    #      x1, x2 >= 0
    
    c = [2.0, 3.0]
    A = sparse([1.0 1.0; 2.0 1.0])
    b = [3.5, 5.0]
    Aeq = spzeros(0, 2)
    beq = Float64[]
    lb = [0.0, 0.0]
    ub = [100.0, 100.0]  # Use reasonable finite upper bounds
    vtype = ['I', 'C']
    
    # Debug: print the problem setup
    println("Problem setup:")
    println("  Objective: max $(c[1])*x1 + $(c[2])*x2")
    println("  Constraints:")
    println("    $(A[1,1])*x1 + $(A[1,2])*x2 <= $(b[1])")
    println("    $(A[2,1])*x1 + $(A[2,2])*x2 <= $(b[2])")
    println("  Bounds: x1 >= $(lb[1]), x2 >= $(lb[2])")
    println("  Types: x1 = $(vtype[1]), x2 = $(vtype[2])")
    
    result = solve_milp_sparse(c, A, b, Aeq, beq, lb, ub, vtype, sense="max", verbose=false)
    
    println("Solution: ", result[:x])
    println("Objective: ", result[:objval])
    println("Status: ", result[:status_name])
    println("Runtime: ", result[:runtime])
    println("MIP Gap: ", result[:mip_gap])
    println("Expected solution: [1.0, 2.5] (approximately)")
    println("Expected objective: 9.5")
end

# # Additional test cases from Gurobi version
# function test_sparse_matrix()
#     println("\n=== Test Case: Sparse Matrix Problem with HiGHS ===")
#     # Create a smaller, more manageable sparse constraint matrix
#     n_vars = 20
#     n_constrs = 10
#     density = 0.3  # Higher density for more realistic constraints
    
#     # Generate random sparse matrix with guaranteed feasible solution
#     A_sparse = sprand(n_constrs, n_vars, density)
#     # Make RHS generous to ensure feasibility
#     b_sparse = rand(n_constrs) * 20 .+ 10
#     Aeq_sparse = spzeros(0, n_vars)
#     beq_sparse = Float64[]
#     cobj_sparse = rand(n_vars)
#     lb_sparse = zeros(n_vars)
#     ub_sparse = fill(10.0, n_vars)  # Reasonable upper bounds
#     vtype_sparse = ['C' for _ in 1:n_vars]  # All continuous variables
    
#     println("Testing with sparse matrix ($(n_constrs)x$(n_vars), $(round(100*density, digits=1))% density)")
#     println("Non-zero elements: $(nnz(A_sparse))")
    
#     result_sparse = @time solve_milp_sparse(cobj_sparse, A_sparse, b_sparse, Aeq_sparse, beq_sparse, lb_sparse, ub_sparse, vtype_sparse, sense="max")
#     println("Sparse matrix test completed. Status: $(result_sparse[:status_name])")
#     println("Objective: $(result_sparse[:objval])")
#     println("Runtime: $(result_sparse[:runtime])")
# end

# function test_binary_problem()
#     println("\n=== Test Case: Binary Problem with HiGHS ===")
#     # max x1 + x2 + 2*x3
#     # s.t. x1 + x2 + x3 <= 3  (allows [1,1,1])
#     #      x1, x2, x3 binary
    
#     cobj = [1.0, 1.0, 2.0]
#     A = sparse([1.0 1.0 1.0; -1.0 -1.0 0.0])  # Single constraint that allows [1,1,1]
#     b = [3.0; -1.0]  # Right-hand side that allows [1,1,1]
#     Aeq = spzeros(0, 3)
#     beq = Float64[]
#     lb = [0.0, 0.0, 0.0]
#     ub = [1.0, 1.0, 1.0]
#     vtype = ['B', 'B', 'B']  # All binary variables
    
#     println("Problem setup:")
#     println("  Objective: max $(cobj[1])*x1 + $(cobj[2])*x2 + $(cobj[3])*x3")
#     println("  Constraints:")
#     println("    $(A[1,1])*x1 + $(A[1,2])*x2 + $(A[1,3])*x3 <= $(b[1])")
#     println("    $(-A[2,1])*x1 - $(A[2,2])*x2 <= $(b[2])")
#     println("  Bounds: 0 <= x1,x2,x3 <= 1")
#     println("  Types: all binary")
#     println("  Expected optimal solution: [1, 1, 1] with objective = 4")
    
#     result = @time solve_milp_sparse(cobj, A, b, Aeq, beq, lb, ub, vtype, sense="max", verbose=true)
    
#     println("Solution: ", result[:x])
#     println("Objective: ", result[:objval])
#     println("Status: ", result[:status_name])
#     println("Runtime: ", result[:runtime])
#     println("MIP Gap: ", result[:mip_gap])
    
#     # Manual verification
#     if result[:status] == kHighsModelStatusOptimal
#         x = result[:x]
#         manual_obj = cobj[1]*x[1] + cobj[2]*x[2] + cobj[3]*x[3]
#         println("Manual objective calculation: $manual_obj")
#         println("Constraint check: $(A[1,1]*x[1] + A[1,2]*x[2] + A[1,3]*x[3]) <= $(b[1])")
        
#         # Check if solution matches expected
#         expected = [1.0, 1.0, 1.0]
#         if isapprox(x, expected, atol=1e-6)
#             println("✓ Solution matches expected [1, 1, 1]")
#         else
#             println("✗ Solution differs from expected [1, 1, 1]")
#         end
#     end
# end

# function test_integer_problem()
#     println("\n=== Test Case: Integer Problem with HiGHS ===")
#     # Type II from Gurobi tests - but ensure feasibility
#     A = sparse([7.0 -2.0; 0.0 1.0; 2.0 -2.0])
#     b = [14.0, 3.0, 10.0]  # Make third constraint less restrictive
#     Aeq = spzeros(0, 2)
#     beq = Float64[]
#     cobj = [4.0, -1.0]
#     lb = zeros(2)
#     ub = [10.0, 10.0]  # Use finite upper bounds
#     vtype = ['I', 'I']  # All integer variables
    
#     result = @time solve_milp_sparse(cobj, A, b, Aeq, beq, lb, ub, vtype, sense="max", verbose=true)
    
#     println("Solution: ", result[:x])
#     println("Objective: ", result[:objval])
#     println("Status: ", result[:status_name])
#     println("Runtime: ", result[:runtime])
#     println("MIP Gap: ", result[:mip_gap])
# end

# Run all tests
test_mixed_integer()
# test_binary_problem()
# test_sparse_matrix()
# test_integer_problem()
