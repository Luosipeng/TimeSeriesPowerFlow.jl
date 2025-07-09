"""
Core type definitions for the OptimalPowerFlow module.
"""

# ============================================================================
# Abstract Type Hierarchy
# ============================================================================

"""
Abstract base type for all OPF formulations.
Enables type-safe dispatch and compile-time specialization.
"""
abstract type AbstractOPFFormulation end

"""
Abstract base type for OPF problems.
Parameterized by formulation type and numeric type for performance.
"""
abstract type AbstractOPFProblem{F<:AbstractOPFFormulation, T<:Number} end

"""
Abstract base type for OPF solvers.
"""
abstract type AbstractOPFSolver end

# ============================================================================
# Concrete Formulation Types
# ============================================================================

"""
DC Optimal Power Flow formulation.
Linear approximation with voltage angles only.
"""
struct DCFormulation <: AbstractOPFFormulation end

"""
AC Optimal Power Flow formulation.
Full nonlinear AC power flow equations.
"""
struct ACFormulation <: AbstractOPFFormulation end

# ============================================================================
# Problem Structure
# ============================================================================

"""
Unified OPF problem structure.
Contains all problem data in a type-safe, performance-oriented format.
"""
struct OPFProblem{F<:AbstractOPFFormulation, T<:Number} <: AbstractOPFProblem{F, T}
    formulation::F
    
    # System dimensions
    nb::Int                     # number of buses
    ng::Int                     # number of generators  
    nl::Int                     # number of lines
    
    # Base MVA for per-unit conversion
    baseMVA::Float64
    
    # Network matrices (type-parameterized for performance)
    Ybus::SparseMatrixCSC{Complex{T}, Int}     # Bus admittance matrix
    Yf::SparseMatrixCSC{Complex{T}, Int}       # From-bus admittance matrix
    Yt::SparseMatrixCSC{Complex{T}, Int}       # To-bus admittance matrix
    Cg::SparseMatrixCSC{T, Int}                # Generator connection matrix
    
    # System data
    bus_data::Matrix{T}
    gen_data::Matrix{T}
    branch_data::Matrix{T}
    
    # Reference indices
    ref_buses::Vector{Int}
    pv_buses::Vector{Int}
    pq_buses::Vector{Int}
    active_gens::Vector{Int}
    
    # Variable bounds
    var_bounds::NamedTuple
    
    # Cost coefficients
    cost_coeffs::NamedTuple
    
    # Internal mapping
    int_to_ext::Vector{Int}
end

# ============================================================================
# Solution Structure
# ============================================================================

"""
Unified OPF solution structure.
Type-safe container for solution data across all formulations.
"""
struct OPFSolution{F<:AbstractOPFFormulation, T<:Number}
    formulation::F
    
    # Solution status
    success::Bool
    objective::T
    iterations::Int
    solve_time::Float64
    
    # Solution variables (formulation-dependent)
    variables::NamedTuple
    
    # Constraint violations
    max_violation::T
    constraint_violations::NamedTuple
    
    # Additional solver information
    solver_info::Dict{String, Any}
end

# ============================================================================
# Result Structure (Public API)
# ============================================================================

"""
Public result structure returned by solve_opf.
Contains solution and updated case data.
"""
struct OPFResult{F<:AbstractOPFFormulation}
    formulation::F
    success::Bool
    objective::Float64
    iterations::Int
    solve_time::Float64
    
    # Updated case data
    case_data::Dict{String, Any}
    
    # Solution details
    solution::OPFSolution{F}
    
    # Performance metrics
    memory_usage::Int
    constraint_violation::Float64
end
