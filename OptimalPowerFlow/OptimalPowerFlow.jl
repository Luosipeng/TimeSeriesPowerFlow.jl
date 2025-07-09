"""
OptimalPowerFlow.jl

A unified, high-performance Optimal Power Flow module for Julia.
Supports both AC and DC formulations with a modular, extensible design.

Features:
- Type-safe formulation system using Julia's type hierarchy
- Zero-cost abstractions through traits and multiple dispatch
- Unified configuration and validation
- Consistent APIs across all formulations
- High performance through compile-time specialization
"""
module OptimalPowerFlow

using LinearAlgebra, SparseArrays, Statistics

# Import PowerFlow module if it exists
try
    using ..PowerFlow
catch
    # If PowerFlow is not available as a parent module, try to include it
    if isfile("../src/PowerFlow.jl")
        include("../src/PowerFlow.jl")
        using .PowerFlow
    end
end

# Core exports - main API
export solve_opf, OPFOptions, OPFResult
export DCFormulation, ACFormulation
export validate_problem, get_problem_info

# Formulation trait exports
export has_voltage_magnitudes, has_reactive_power, uses_complex_matrices
export is_linear_constraints, get_variable_count, get_constraint_count

# Advanced exports for extension
export AbstractOPFFormulation, AbstractOPFProblem, AbstractOPFSolver
export create_problem, create_solver, extract_solution

# Utility exports
export profile_solve, compare_formulations, validate_solution, export_results
export from_legacy_options, to_legacy_options, rundcopf, runacopf

# Check if components are already loaded, if not include them
if !@isdefined(AbstractOPFFormulation)
    include("types.jl")
end

if !@isdefined(has_voltage_magnitudes)
    include("traits.jl")
end

if !@isdefined(OPFOptions)
    include("config.jl")
end

if !@isdefined(validate_problem)
    include("validation.jl")
end

if !@isdefined(create_problem)
    # Include solver dependencies first
    if isfile("../solvers/types.jl")
        include("../solvers/types.jl")
    end
    if isfile("../solvers/interior_point_method.jl")
        include("../solvers/interior_point_method.jl")
    end
    
    include("formulations/dc_formulation.jl")
    include("formulations/ac_formulation.jl")
end

if !@isdefined(solve_opf)
    include("solvers/unified_solver.jl")
end

if !@isdefined(profile_solve)
    include("utils.jl")
end

end # module OptimalPowerFlow
