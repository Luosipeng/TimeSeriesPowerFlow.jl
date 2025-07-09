"""
Trait system for OPF formulations.
Enables compile-time specialization and zero-cost abstractions.
"""

# ============================================================================
# Formulation Traits
# ============================================================================

"""
Check if formulation includes voltage magnitude variables.
"""
has_voltage_magnitudes(::Type{DCFormulation}) = false
has_voltage_magnitudes(::Type{ACFormulation}) = true
has_voltage_magnitudes(::Type{F}) where {F<:AbstractOPFFormulation} = error("Trait not implemented for $F")

"""
Check if formulation includes reactive power variables.
"""
has_reactive_power(::Type{DCFormulation}) = false
has_reactive_power(::Type{ACFormulation}) = true
has_reactive_power(::Type{F}) where {F<:AbstractOPFFormulation} = error("Trait not implemented for $F")

"""
Check if formulation uses complex matrix operations.
"""
uses_complex_matrices(::Type{DCFormulation}) = false
uses_complex_matrices(::Type{ACFormulation}) = true
uses_complex_matrices(::Type{F}) where {F<:AbstractOPFFormulation} = error("Trait not implemented for $F")

"""
Check if formulation has linear constraints only.
"""
is_linear_constraints(::Type{DCFormulation}) = true
is_linear_constraints(::Type{ACFormulation}) = false
is_linear_constraints(::Type{F}) where {F<:AbstractOPFFormulation} = error("Trait not implemented for $F")

"""
Check if formulation requires complex voltage variables.
"""
requires_complex_voltage(::Type{DCFormulation}) = false
requires_complex_voltage(::Type{ACFormulation}) = true
requires_complex_voltage(::Type{F}) where {F<:AbstractOPFFormulation} = error("Trait not implemented for $F")

# ============================================================================
# Variable Count Traits (Multiple Dispatch)
# ============================================================================

"""
Get total number of optimization variables for DC formulation.
Variables: [va; pij; pg] where va excludes reference bus
"""
function get_variable_count(::Type{DCFormulation}, nb::Int, ng::Int, nl::Int)
    return nb + nl + ng  # angles + line_flows + gen_active_power
end

"""
Get total number of optimization variables for AC formulation.
Variables: [va; vm; pg; qg] where va excludes reference bus
"""
function get_variable_count(::Type{ACFormulation}, nb::Int, ng::Int, nl::Int)
    return nb + nb + ng + ng  # angles + magnitudes + gen_active + gen_reactive
end

"""
Get variable count from problem instance.
"""
get_variable_count(problem::OPFProblem{F}) where {F} = 
    get_variable_count(F, problem.nb, problem.ng, problem.nl)

# ============================================================================
# Constraint Count Traits
# ============================================================================

"""
Get number of equality constraints for DC formulation.
"""
function get_equality_constraint_count(::Type{DCFormulation}, nb::Int, ng::Int, nl::Int)
    return nb + nl  # power_balance + line_flow_def
end

"""
Get number of equality constraints for AC formulation.
"""
function get_equality_constraint_count(::Type{ACFormulation}, nb::Int, ng::Int, nl::Int)
    return 2 * nb  # active_power_balance + reactive_power_balance
end

"""
Get number of inequality constraints (bounds + operational limits).
"""
function get_inequality_constraint_count(::Type{F}, nb::Int, ng::Int, nl::Int) where {F}
    nvar = get_variable_count(F, nb, ng, nl)
    return 2 * nvar + nl  # variable_bounds + line_flow_limits
end

# ============================================================================
# Numeric Type Traits
# ============================================================================

"""
Get preferred numeric type for formulation.
"""
preferred_numeric_type(::Type{DCFormulation}) = Float64
preferred_numeric_type(::Type{ACFormulation}) = Float64

"""
Get matrix element type for formulation.
"""
matrix_element_type(::Type{DCFormulation}) = Float64
matrix_element_type(::Type{ACFormulation}) = ComplexF64

# ============================================================================
# Variable Indexing Traits
# ============================================================================

"""
Get variable indices for DC formulation.
Returns NamedTuple with symbolic variable ranges.
"""
function get_variable_indices(::Type{DCFormulation}, nb::Int, ng::Int, nl::Int)
    va_end = nb
    pij_end = va_end + nl
    pg_end = pij_end + ng
    
    return (
        va = 1:va_end,
        pij = (va_end + 1):pij_end,
        pg = (pij_end + 1):pg_end,
        total = pg_end
    )
end

"""
Get variable indices for AC formulation.
"""
function get_variable_indices(::Type{ACFormulation}, nb::Int, ng::Int, nl::Int)
    va_end = nb
    vm_end = va_end + nb
    pg_end = vm_end + ng
    qg_end = pg_end + ng
    
    return (
        va = 1:va_end,
        vm = (va_end + 1):vm_end,
        pg = (vm_end + 1):pg_end,
        qg = (pg_end + 1):qg_end,
        total = qg_end
    )
end

# ============================================================================
# Solver Compatibility Traits
# ============================================================================

"""
Check if formulation is compatible with linear solvers.
"""
is_linear_compatible(::Type{DCFormulation}) = true
is_linear_compatible(::Type{ACFormulation}) = false

"""
Check if formulation requires nonlinear solver.
"""
requires_nonlinear_solver(::Type{DCFormulation}) = false
requires_nonlinear_solver(::Type{ACFormulation}) = true

"""
Get recommended solver tolerance for formulation.
"""
default_tolerance(::Type{DCFormulation}) = 1e-6
default_tolerance(::Type{ACFormulation}) = 1e-5

"""
Get recommended maximum iterations for formulation.
"""
default_max_iterations(::Type{DCFormulation}) = 100
default_max_iterations(::Type{ACFormulation}) = 200
