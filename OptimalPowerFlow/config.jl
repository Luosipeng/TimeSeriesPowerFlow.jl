"""
Unified configuration system for OptimalPowerFlow module.
Provides type-safe, validated configuration with sensible defaults.
"""

# ============================================================================
# Configuration Structures
# ============================================================================

"""
Unified OPF options structure.
Uses keyword arguments with sensible defaults for all formulations.
"""
Base.@kwdef struct OPFOptions
    # Solver settings
    max_iterations::Int = 200
    tolerance::Float64 = 1e-5
    verbose::Int = 1
    
    # Physical limits
    angle_limit_deg::Float64 = 60.0
    voltage_min_default::Float64 = 0.95
    voltage_max_default::Float64 = 1.05
    
    # Validation settings
    enable_validation::Bool = true
    validation_tolerance::Float64 = 1e-8
    bounds_margin::Float64 = 1e-6
    
    # Performance settings
    use_sparse_matrices::Bool = true
    preallocate_memory::Bool = true
    
    # Line flow settings
    line_capacity_factor::Float64 = 2.0
    default_line_limit_fraction::Float64 = 0.9
    
    # Generator settings
    min_power_gap::Float64 = 0.01
    default_q_limit_fraction::Float64 = 0.5
    
    # Solver-specific settings
    interior_point_mu::Float64 = 0.1
    barrier_parameter::Float64 = 0.995
    feasibility_tolerance::Float64 = 1e-8
    
    # Advanced options
    formulation_specific::Dict{String, Any} = Dict{String, Any}()
end

# ============================================================================
# Configuration Validation
# ============================================================================

"""
Validate and normalize OPF options.
"""
function validate_options(options::OPFOptions, formulation::Type{F}) where {F<:AbstractOPFFormulation}
    validated_opts = deepcopy(options)
    
    # Validate basic constraints
    validated_opts.max_iterations < 1 && throw(ArgumentError("max_iterations must be positive"))
    validated_opts.tolerance <= 0 && throw(ArgumentError("tolerance must be positive"))
    validated_opts.angle_limit_deg <= 0 && throw(ArgumentError("angle_limit_deg must be positive"))
    
    # Voltage bounds validation
    if validated_opts.voltage_min_default >= validated_opts.voltage_max_default
        throw(ArgumentError("voltage_min_default must be less than voltage_max_default"))
    end
    
    # Formulation-specific validation
    if F == DCFormulation
        # DC formulation doesn't use voltage magnitude limits
        if haskey(validated_opts.formulation_specific, "ignore_voltage_limits")
            validated_opts.formulation_specific["ignore_voltage_limits"] = true
        end
    elseif F == ACFormulation
        # Ensure reasonable voltage bounds for AC
        validated_opts.voltage_min_default < 0.5 && @warn "Very low voltage minimum: $(validated_opts.voltage_min_default)"
        validated_opts.voltage_max_default > 1.5 && @warn "Very high voltage maximum: $(validated_opts.voltage_max_default)"
    end
    
    # Adjust tolerance based on formulation
    if validated_opts.tolerance < default_tolerance(F) / 10
        @warn "Tolerance $(validated_opts.tolerance) is very strict for $F, consider using $(default_tolerance(F))"
    end
    
    return validated_opts
end

# ============================================================================
# Solver Configuration
# ============================================================================

"""
Convert unified options to solver-specific format.
"""
function get_solver_options(options::OPFOptions, formulation::Type{F}) where {F<:AbstractOPFFormulation}
    base_opts = Dict{String, Any}(
        "max_iterations" => options.max_iterations,
        "tolerance" => options.tolerance,
        "verbose" => options.verbose,
        "feasibility_tolerance" => options.feasibility_tolerance,
        "barrier_parameter" => options.barrier_parameter
    )
    
    # Add formulation-specific options
    if F == DCFormulation
        merge!(base_opts, Dict{String, Any}(
            "angle_limit_deg" => options.angle_limit_deg,
            "line_capacity_factor" => options.line_capacity_factor,
            "linear_solver" => true
        ))
    elseif F == ACFormulation
        merge!(base_opts, Dict{String, Any}(
            "angle_limit_deg" => options.angle_limit_deg,
            "voltage_min_default" => options.voltage_min_default,
            "voltage_max_default" => options.voltage_max_default,
            "nonlinear_solver" => true,
            "complex_matrices" => uses_complex_matrices(F)
        ))
    end
    
    # Add any custom formulation-specific options
    merge!(base_opts, options.formulation_specific)
    
    return base_opts
end

# ============================================================================
# Performance Configuration
# ============================================================================

"""
Get performance-optimized configuration for given problem size.
"""
function get_performance_config(nb::Int, ng::Int, nl::Int, formulation::Type{F}) where {F}
    if nb < 50  # Small system
        return OPFOptions(
            max_iterations = default_max_iterations(F),
            tolerance = default_tolerance(F),
            use_sparse_matrices = false,  # Dense might be faster for small systems
            preallocate_memory = true
        )
    elseif nb < 500  # Medium system
        return OPFOptions(
            max_iterations = default_max_iterations(F) + 50,
            tolerance = default_tolerance(F),
            use_sparse_matrices = true,
            preallocate_memory = true
        )
    else  # Large system
        return OPFOptions(
            max_iterations = default_max_iterations(F) + 100,
            tolerance = default_tolerance(F) * 2,  # Slightly relaxed for large systems
            use_sparse_matrices = true,
            preallocate_memory = true,
            verbose = 0  # Reduce output for large systems
        )
    end
end

# ============================================================================
# Legacy Compatibility
# ============================================================================

"""
Convert legacy option dictionary to unified OPFOptions.
"""
function from_legacy_options(legacy_opts::Dict{String, Any})
    options = OPFOptions()
    
    # Map common legacy options
    options = OPFOptions(
        max_iterations = get(legacy_opts, "OPF_MAX_IT", options.max_iterations),
        tolerance = get(legacy_opts, "OPF_VIOLATION", options.tolerance),
        verbose = get(legacy_opts, "VERBOSE", options.verbose),
        angle_limit_deg = get(legacy_opts, "ANGLE_LIMIT_DEG", options.angle_limit_deg),
        voltage_min_default = get(legacy_opts, "DEFAULT_VMIN", options.voltage_min_default),
        voltage_max_default = get(legacy_opts, "DEFAULT_VMAX", options.voltage_max_default),
        line_capacity_factor = get(legacy_opts, "AVG_LINE_FACTOR", options.line_capacity_factor)
    )
    
    return options
end

"""
Convert unified options to legacy format for backward compatibility.
"""
function to_legacy_options(options::OPFOptions, formulation::Type{F}) where {F}
    legacy_opts = Dict{String, Any}(
        "OPF_MAX_IT" => options.max_iterations,
        "OPF_VIOLATION" => options.tolerance,
        "VERBOSE" => options.verbose,
        "ANGLE_LIMIT_DEG" => options.angle_limit_deg
    )
    
    if F == DCFormulation
        legacy_opts["AVG_LINE_FACTOR"] = options.line_capacity_factor
    elseif F == ACFormulation
        legacy_opts["DEFAULT_VMIN"] = options.voltage_min_default
        legacy_opts["DEFAULT_VMAX"] = options.voltage_max_default
    end
    
    return legacy_opts
end
