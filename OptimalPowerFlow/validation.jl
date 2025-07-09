"""
Unified validation system for OPF problems.
Provides comprehensive bounds checking and problem feasibility analysis.
"""

# ============================================================================
# Validation Errors
# ============================================================================

abstract type OPFValidationError <: Exception end

struct InfeasibleBoundsError <: OPFValidationError
    message::String
    component::String
    indices::Vector{Int}
end

struct InvalidDataError <: OPFValidationError
    message::String
    component::String
end

struct TopologyError <: OPFValidationError
    message::String
end

# ============================================================================
# Problem Validation
# ============================================================================

"""
Comprehensive validation of OPF problem data.
"""
function validate_problem(jpc::Dict{String, Any}, options::OPFOptions, formulation::Type{F}) where {F}
    if !options.enable_validation
        return jpc  # Skip validation if disabled
    end
    
    validated_jpc = deepcopy(jpc)
    
    try
        # Basic data structure validation
        validate_data_structure!(validated_jpc)
        
        # Network topology validation
        validate_topology!(validated_jpc)
        
        # Bounds validation (formulation-specific)
        validate_bounds!(validated_jpc, options, formulation)
        
        # Feasibility checks
        validate_feasibility!(validated_jpc, options, formulation)
        
        return validated_jpc
        
    catch e
        if isa(e, OPFValidationError)
            rethrow(e)
        else
            throw(InvalidDataError("Unexpected validation error: $(e)", "unknown"))
        end
    end
end

"""
Validate basic data structure integrity.
"""
function validate_data_structure!(jpc::Dict{String, Any})
    # Check required fields
    required_fields = ["baseMVA", "busAC", "genAC", "branchAC"]
    for field in required_fields
        if !haskey(jpc, field)
            throw(InvalidDataError("Missing required field: $field", "structure"))
        end
    end
    
    # Check data types and dimensions
    bus = jpc["busAC"]
    gen = jpc["genAC"]
    branch = jpc["branchAC"]
    
    if !isa(bus, Matrix) || size(bus, 2) < 13
        throw(InvalidDataError("Invalid bus matrix format", "bus"))
    end
    
    if !isa(gen, Matrix) || size(gen, 2) < 21
        throw(InvalidDataError("Invalid generator matrix format", "gen"))
    end
    
    if !isa(branch, Matrix) || size(branch, 2) < 13
        throw(InvalidDataError("Invalid branch matrix format", "branch"))
    end
    
    # Check for NaN or Inf values
    for (name, data) in [("bus", bus), ("gen", gen), ("branch", branch)]
        if any(isnan, data) || any(isinf, data)
            throw(InvalidDataError("NaN or Inf values detected in $name data", name))
        end
    end
end

"""
Validate network topology.
"""
function validate_topology!(jpc::Dict{String, Any})
    bus = jpc["busAC"]
    gen = jpc["genAC"]
    branch = jpc["branchAC"]
    
    nb = size(bus, 1)
    
    # Check bus numbering
    bus_numbers = Int.(bus[:, 1])  # BUS_I
    if !allunique(bus_numbers)
        throw(TopologyError("Duplicate bus numbers detected"))
    end
    
    # Check generator bus connections
    gen_buses = Int.(gen[:, 1])  # GEN_BUS
    invalid_gen_buses = setdiff(gen_buses, bus_numbers)
    if !isempty(invalid_gen_buses)
        throw(TopologyError("Generators connected to non-existent buses: $invalid_gen_buses"))
    end
    
    # Check branch connections
    from_buses = Int.(branch[:, 1])  # F_BUS
    to_buses = Int.(branch[:, 2])    # T_BUS
    
    invalid_from = setdiff(from_buses, bus_numbers)
    invalid_to = setdiff(to_buses, bus_numbers)
    
    if !isempty(invalid_from) || !isempty(invalid_to)
        throw(TopologyError("Branches connected to non-existent buses"))
    end
    
    # Check for reference bus
    bus_types = Int.(bus[:, 2])  # BUS_TYPE
    ref_buses = findall(bus_types .== 3)  # REF type
    
    if isempty(ref_buses)
        throw(TopologyError("No reference bus found"))
    end
    
    if length(ref_buses) > 1
        @warn "Multiple reference buses found, using first one"
    end
end

"""
Validate and fix variable bounds (formulation-specific).
"""
function validate_bounds!(jpc::Dict{String, Any}, options::OPFOptions, formulation::Type{F}) where {F}
    bus = jpc["busAC"]
    gen = jpc["genAC"]
    branch = jpc["branchAC"]
    baseMVA = jpc["baseMVA"]
    
    # Validate voltage bounds (AC formulation only)
    if has_voltage_magnitudes(F)
        validate_voltage_bounds!(bus, options)
    end
    
    # Validate generator bounds
    validate_generator_bounds!(gen, baseMVA, options)
    
    # Validate line limits
    validate_line_limits!(branch, baseMVA, options)
    
    # Update jpc with validated data
    jpc["busAC"] = bus
    jpc["genAC"] = gen
    jpc["branchAC"] = branch
end

"""
Validate voltage magnitude bounds.
"""
function validate_voltage_bounds!(bus::Matrix, options::OPFOptions)
    nb = size(bus, 1)
    
    # Get voltage limits (columns 12, 13 are VMIN, VMAX)
    if size(bus, 2) >= 13
        vmin = bus[:, 12]
        vmax = bus[:, 13]
        
        # Set defaults for missing limits
        zero_min = findall(vmin .<= 0)
        zero_max = findall(vmax .<= 0)
        
        if !isempty(zero_min)
            bus[zero_min, 12] .= options.voltage_min_default
        end
        
        if !isempty(zero_max)
            bus[zero_max, 13] .= options.voltage_max_default
        end
        
        # Check bounds consistency
        vmin = bus[:, 12]
        vmax = bus[:, 13]
        
        inconsistent = findall(vmin .>= vmax)
        if !isempty(inconsistent)
            # Fix by adding margin
            margin = options.bounds_margin
            bus[inconsistent, 13] = bus[inconsistent, 12] .+ margin
            @warn "Fixed $(length(inconsistent)) inconsistent voltage bounds"
        end
    else
        # Add voltage limit columns if missing
        bus_new = zeros(nb, 13)
        bus_new[:, 1:size(bus, 2)] = bus
        bus_new[:, 12] .= options.voltage_min_default  # VMIN
        bus_new[:, 13] .= options.voltage_max_default  # VMAX
        bus = bus_new
    end
    
    return bus
end

"""
Validate generator power bounds.
"""
function validate_generator_bounds!(gen::Matrix, baseMVA::Float64, options::OPFOptions)
    ng = size(gen, 1)
    
    # Active power bounds (columns 9, 10 are PMIN, PMAX)
    if size(gen, 2) >= 10
        pmin = gen[:, 9]
        pmax = gen[:, 10]
        
        # Set defaults for missing bounds
        zero_min = findall(pmin .== 0)
        zero_max = findall(pmax .== 0)
        
        if !isempty(zero_max)
            gen[zero_max, 10] .= baseMVA  # Default 1 p.u.
        end
        
        # Check consistency
        inconsistent = findall(pmin .>= pmax)
        if !isempty(inconsistent)
            gen[inconsistent, 10] = gen[inconsistent, 9] .+ options.min_power_gap * baseMVA
        end
    end
    
    # Reactive power bounds (columns 4, 5 are QMIN, QMAX)
    if size(gen, 2) >= 5
        qmin = gen[:, 4]
        qmax = gen[:, 5]
        
        zero_min = findall(qmin .== 0)
        zero_max = findall(qmax .== 0)
        
        if !isempty(zero_min)
            gen[zero_min, 4] .= -options.default_q_limit_fraction * baseMVA
        end
        
        if !isempty(zero_max)
            gen[zero_max, 5] .= options.default_q_limit_fraction * baseMVA
        end
        
        # Check consistency
        qmin = gen[:, 4]
        qmax = gen[:, 5]
        inconsistent = findall(qmin .>= qmax)
        if !isempty(inconsistent)
            gen[inconsistent, 5] = gen[inconsistent, 4] .+ options.min_power_gap * baseMVA
        end
    end
end

"""
Validate transmission line limits.
"""
function validate_line_limits!(branch::Matrix, baseMVA::Float64, options::OPFOptions)
    nl = size(branch, 1)
    
    if size(branch, 2) >= 6  # RATE_A column
        rate_a = branch[:, 6]
        
        # Set defaults for missing limits
        zero_limits = findall(rate_a .<= 0)
        if !isempty(zero_limits)
            # Calculate intelligent defaults based on impedance
            for i in zero_limits
                x = branch[i, 4]  # BR_X
                if x > 0
                    # Thermal limit based on impedance
                    thermal_limit = options.line_capacity_factor * baseMVA / x
                    branch[i, 6] = min(thermal_limit, 10 * baseMVA)  # Cap at 10 p.u.
                else
                    branch[i, 6] = baseMVA  # Conservative default
                end
            end
        end
    else
        # Add RATE_A column if missing
        branch_new = zeros(nl, max(6, size(branch, 2)))
        branch_new[:, 1:size(branch, 2)] = branch
        branch_new[:, 6] .= baseMVA  # Default 1 p.u. limit
        branch = branch_new
    end
    
    return branch
end

"""
Check overall problem feasibility.
"""
function validate_feasibility!(jpc::Dict{String, Any}, options::OPFOptions, formulation::Type{F}) where {F}
    bus = jpc["busAC"]
    gen = jpc["genAC"]
    baseMVA = jpc["baseMVA"]
    
    # Active power balance check
    total_load = sum(bus[:, 3]) / baseMVA  # PD column
    active_gens = findall(gen[:, 8] .> 0)  # GEN_STATUS
    total_gen_capacity = sum(gen[active_gens, 10]) / baseMVA  # PMAX
    
    if total_load > total_gen_capacity * 0.95  # 5% margin
        throw(InfeasibleBoundsError(
            "Insufficient generation capacity: Load=$(total_load) > Capacity=$(total_gen_capacity)",
            "generation", active_gens
        ))
    end
    
    # For AC formulation, also check reactive power
    if has_reactive_power(F)
        total_q_load = sum(bus[:, 4]) / baseMVA  # QD column
        total_q_capacity = sum(gen[active_gens, 5]) / baseMVA  # QMAX
        
        if total_q_load > total_q_capacity * 0.95
            @warn "Tight reactive power capacity: Load=$(total_q_load), Capacity=$(total_q_capacity)"
        end
    end
end

"""
Get comprehensive problem information.
"""
function get_problem_info(jpc::Dict{String, Any}, formulation::Type{F}) where {F}
    bus = jpc["busAC"]
    gen = jpc["genAC"] 
    branch = jpc["branchAC"]
    
    nb = size(bus, 1)
    ng = size(gen, 1)
    nl = size(branch, 1)
    
    active_gens = sum(gen[:, 8] .> 0)  # GEN_STATUS
    ref_buses = sum(bus[:, 2] .== 3)   # REF type
    
    info = Dict{String, Any}(
        "formulation" => string(F),
        "dimensions" => Dict(
            "buses" => nb,
            "generators" => ng,
            "active_generators" => active_gens,
            "branches" => nl,
            "reference_buses" => ref_buses
        ),
        "variables" => get_variable_count(F, nb, active_gens, nl),
        "equality_constraints" => get_equality_constraint_count(F, nb, active_gens, nl),
        "inequality_constraints" => get_inequality_constraint_count(F, nb, active_gens, nl),
        "properties" => Dict(
            "has_voltage_magnitudes" => has_voltage_magnitudes(F),
            "has_reactive_power" => has_reactive_power(F),
            "uses_complex_matrices" => uses_complex_matrices(F),
            "is_linear" => is_linear_constraints(F)
        )
    )
    
    return info
end
