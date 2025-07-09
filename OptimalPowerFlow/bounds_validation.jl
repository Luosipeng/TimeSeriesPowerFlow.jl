"""
Bounds Validation Module for OPF Problems

This module provides robust handling of boundary conditions for OPF problems,
ensuring feasibility even when input data is incomplete or inconsistent.
"""

using ..PowerFlow: REF, PV, PQ, BUS_TYPE, VM, VMIN, VMAX
using ..PowerFlow: GEN_STATUS, PMIN, PMAX, QMIN, QMAX, PG, QG
using ..PowerFlow: RATE_A, BR_X

"""
    validate_voltage_bounds!(bus, opt)

Validate and correct voltage magnitude bounds to ensure feasibility.
"""
function validate_voltage_bounds!(bus, opt)
    nb = size(bus, 1)
    
    # Default voltage bounds if not specified
    default_vmin = get(opt, "DEFAULT_VMIN", 0.95)
    default_vmax = get(opt, "DEFAULT_VMAX", 1.05)
    
    # Check if VMIN/VMAX columns exist
    if size(bus, 2) < VMAX
        # Expand bus matrix if needed
        bus_new = zeros(size(bus, 1), VMAX)
        bus_new[:, 1:size(bus, 2)] = bus
        bus = bus_new
    end
    
    # Set default bounds where missing or invalid
    invalid_min = (bus[:, VMIN] .<= 0) .| (bus[:, VMIN] .>= bus[:, VMAX])
    invalid_max = (bus[:, VMAX] .<= 0) .| (bus[:, VMAX] .<= bus[:, VMIN])
    
    bus[invalid_min, VMIN] .= default_vmin
    bus[invalid_max, VMAX] .= default_vmax
    
    # Ensure VMIN < VMAX with minimum gap
    min_gap = get(opt, "MIN_VOLTAGE_GAP", 0.05)
    tight_bounds = (bus[:, VMAX] - bus[:, VMIN]) .< min_gap
    if any(tight_bounds)
        bus[tight_bounds, VMAX] = bus[tight_bounds, VMIN] .+ min_gap
    end
    
    return bus
end

"""
    validate_generator_bounds!(gen, baseMVA, opt)

Validate and correct generator power bounds to ensure feasibility.
"""
function validate_generator_bounds!(gen, baseMVA, opt)
    ng = size(gen, 1)
    on = findall(gen[:, GEN_STATUS] .> 0)
    
    if isempty(on)
        error("No generators are online")
    end
    
    # Default generator bounds as fraction of base MVA
    default_pmin_fraction = get(opt, "DEFAULT_PMIN_FRACTION", 0.0)
    default_pmax_fraction = get(opt, "DEFAULT_PMAX_FRACTION", 1.0)
    default_qmin_fraction = get(opt, "DEFAULT_QMIN_FRACTION", -0.5)
    default_qmax_fraction = get(opt, "DEFAULT_QMAX_FRACTION", 0.5)
    
    # Validate P bounds
    invalid_pmin = (gen[on, PMIN] .< 0) .| (gen[on, PMIN] .>= gen[on, PMAX])
    invalid_pmax = gen[on, PMAX] .<= 0
    
    gen[on[invalid_pmin], PMIN] .= default_pmin_fraction * baseMVA
    gen[on[invalid_pmax], PMAX] .= default_pmax_fraction * baseMVA
    
    # Validate Q bounds
    invalid_qmin = gen[on, QMIN] .>= gen[on, QMAX]
    invalid_qmax = gen[on, QMAX] .<= gen[on, QMIN]
    
    gen[on[invalid_qmin], QMIN] .= default_qmin_fraction * baseMVA
    gen[on[invalid_qmax], QMAX] .= default_qmax_fraction * baseMVA
    
    # Ensure minimum gap between bounds
    min_gap = get(opt, "MIN_POWER_GAP", 0.01) * baseMVA
    
    # P bounds gap
    tight_p = (gen[on, PMAX] - gen[on, PMIN]) .< min_gap
    if any(tight_p)
        gen[on[tight_p], PMAX] = gen[on[tight_p], PMIN] .+ min_gap
    end
    
    # Q bounds gap
    tight_q = (gen[on, QMAX] - gen[on, QMIN]) .< min_gap
    if any(tight_q)
        gen[on[tight_q], QMAX] = gen[on[tight_q], QMIN] .+ min_gap
    end
    
    return gen
end

"""
    validate_line_limits!(branch, baseMVA, opt)

Validate and correct line flow limits to ensure reasonable bounds.
"""
function validate_line_limits!(branch, baseMVA, opt)
    nl = size(branch, 1)
    
    # Ensure RATE_A column exists
    if size(branch, 2) < RATE_A
        branch_new = zeros(size(branch, 1), RATE_A)
        branch_new[:, 1:size(branch, 2)] = branch
        branch = branch_new
    end
    
    # Default line limits based on impedance
    default_rate_multiplier = get(opt, "DEFAULT_RATE_MULTIPLIER", 5.0)
    
    for i in 1:nl
        if branch[i, RATE_A] <= 0
            # Estimate capacity based on impedance
            x = branch[i, BR_X]
            if x > 0
                # Simple estimate: higher capacity for lower impedance
                estimated_rate = default_rate_multiplier * baseMVA / x
                branch[i, RATE_A] = min(estimated_rate, 10 * baseMVA)  # Cap at 10x base
            else
                branch[i, RATE_A] = baseMVA  # Default to base MVA
            end
        end
    end
    
    return branch
end

"""
    validate_feasible_initialization(x0, xmin, xmax, opt)

Ensure initial point is feasible and well-conditioned.
"""
function validate_feasible_initialization(x0, xmin, xmax, opt)
    margin = get(opt, "BOUNDS_MARGIN", 1e-6)
    
    # Identify equality constraints (where xmin == xmax)
    equality_constraints = abs.(xmax - xmin) .< 1e-12
    
    # Check for inconsistent bounds (xmin > xmax) excluding equality constraints
    inconsistent = (xmin .> xmax) .& .!equality_constraints
    if any(inconsistent)
        inconsistent_indices = findall(inconsistent)
        error("Inconsistent bounds detected at variables: $(inconsistent_indices)")
    end
    
    # For equality constraints, set x0 to the constraint value
    x0_feasible = copy(x0)
    x0_feasible[equality_constraints] = xmin[equality_constraints]
    
    # For inequality constraints, project to feasible region with margin
    inequality_constraints = .!equality_constraints
    if any(inequality_constraints)
        x0_feasible[inequality_constraints] = max.(
            xmin[inequality_constraints] .+ margin,
            min.(xmax[inequality_constraints] .- margin, 
                 x0_feasible[inequality_constraints])
        )
    end
    
    # Check if projection changed the point significantly (excluding equality constraints)
    if any(inequality_constraints)
        max_change = maximum(abs.(x0_feasible[inequality_constraints] - x0[inequality_constraints]))
        if max_change > get(opt, "INIT_PROJECTION_WARN", 0.1)
            @warn "Initial point projected significantly (max change: $max_change)"
        end
    end
    
    return x0_feasible
end

"""
    check_problem_feasibility(bus, gen, branch, baseMVA, opt)

Perform comprehensive feasibility checks before optimization.
"""
function check_problem_feasibility(bus, gen, branch, baseMVA, opt)
    # Check power balance feasibility
    on = findall(gen[:, GEN_STATUS] .> 0)
    total_load_p = sum(bus[:, PD]) / baseMVA
    total_load_q = sum(bus[:, QD]) / baseMVA
    
    max_gen_p = sum(gen[on, PMAX]) / baseMVA
    min_gen_p = sum(gen[on, PMIN]) / baseMVA
    max_gen_q = sum(gen[on, QMAX]) / baseMVA
    min_gen_q = sum(gen[on, QMIN]) / baseMVA
    
    # Check if generation can meet load
    if max_gen_p < total_load_p * 0.99  # Allow 1% tolerance
        error("Insufficient generation capacity: max=$(max_gen_p), load=$(total_load_p)")
    end
    
    if min_gen_p > total_load_p * 1.01  # Allow 1% tolerance
        error("Minimum generation exceeds load: min=$(min_gen_p), load=$(total_load_p)")
    end
    
    # Reactive power feasibility (less strict due to voltage control)
    if max_gen_q < total_load_q - abs(total_load_q) * 2  # Allow large margin
        @warn "Potential reactive power shortage: max=$(max_gen_q), load=$(total_load_q)"
    end
    
    # Check reference bus
    ref_buses = findall(bus[:, BUS_TYPE] .== REF)
    if isempty(ref_buses)
        error("No reference bus found")
    elseif length(ref_buses) > 1
        @warn "Multiple reference buses found, using first one"
    end
    
    return true
end
