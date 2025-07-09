"""
Utility functions for OptimalPowerFlow module.
"""

# ============================================================================
# Performance Monitoring
# ============================================================================

"""
Monitor and profile OPF solve performance.
"""
function profile_solve(case_data::Dict{String, Any}, 
                      formulation::Type{F},
                      options::OPFOptions = OPFOptions()) where {F}
    
    # Compile-time overhead measurement
    compile_time = @elapsed begin
        # Force compilation by solving a minimal problem first
        minimal_case = create_minimal_case()
        solve_opf(minimal_case, formulation, OPFOptions(verbose=0))
    end
    
    # Actual solve timing
    solve_time = @elapsed result = solve_opf(case_data, formulation, options)
    
    # Memory allocation
    memory_used = @allocated solve_opf(case_data, formulation, options)
    
    # Performance metrics
    info = get_problem_info(case_data, formulation)
    
    perf_info = Dict{String, Any}(
        "compile_time" => compile_time,
        "solve_time" => solve_time,
        "memory_mb" => memory_used / 1024 / 1024,
        "variables" => info["variables"],
        "constraints" => info["equality_constraints"] + info["inequality_constraints"],
        "success" => result.success,
        "iterations" => result.iterations,
        "objective" => result.objective,
        "efficiency" => result.success ? result.objective / solve_time : 0.0
    )
    
    return result, perf_info
end

"""
Create minimal test case for compilation timing.
"""
function create_minimal_case()
    return Dict{String, Any}(
        "baseMVA" => 100.0,
        "busAC" => [1 3 0 0 0 0 1 1 0 0 1 0.95 1.05],  # Simple 1-bus system
        "genAC" => [1 0 0 -100 100 1 100 1 100 0 0 0 0 0 0 0 0 0 0 0 2 0 100],
        "branchAC" => zeros(0, 13)  # No branches
    )
end

# ============================================================================
# Problem Scaling and Conditioning
# ============================================================================

"""
Analyze problem conditioning and suggest scaling.
"""
function analyze_conditioning(case_data::Dict{String, Any}, 
                             formulation::Type{F}) where {F}
    
    bus = case_data["busAC"]
    gen = case_data["genAC"]
    branch = case_data["branchAC"]
    baseMVA = case_data["baseMVA"]
    
    analysis = Dict{String, Any}()
    
    # Load analysis
    if !isempty(bus)
        loads = bus[:, 3]  # PD
        analysis["load_range"] = (minimum(loads), maximum(loads))
        analysis["load_std"] = std(loads)
        analysis["total_load"] = sum(loads)
    end
    
    # Generation analysis
    if !isempty(gen)
        pmax = gen[:, 10]  # PMAX
        analysis["gen_range"] = (minimum(pmax), maximum(pmax))
        analysis["gen_std"] = std(pmax)
        analysis["total_capacity"] = sum(pmax)
    end
    
    # Impedance analysis
    if !isempty(branch)
        x_values = branch[:, 4]  # BR_X
        nonzero_x = x_values[x_values .> 0]
        if !isempty(nonzero_x)
            analysis["impedance_range"] = (minimum(nonzero_x), maximum(nonzero_x))
            analysis["impedance_ratio"] = maximum(nonzero_x) / minimum(nonzero_x)
        end
    end
    
    # Conditioning recommendations
    recommendations = String[]
    
    if haskey(analysis, "impedance_ratio") && analysis["impedance_ratio"] > 1000
        push!(recommendations, "Consider impedance scaling - high ratio detected")
    end
    
    if haskey(analysis, "load_std") && analysis["load_std"] > analysis["total_load"] * 0.5
        push!(recommendations, "Consider load balancing - high variance detected")
    end
    
    analysis["recommendations"] = recommendations
    
    return analysis
end

# ============================================================================
# Comparison and Benchmarking
# ============================================================================

"""
Compare multiple formulations on the same case.
"""
function compare_formulations(case_data::Dict{String, Any}, 
                             formulations::Vector{Type{<:AbstractOPFFormulation}},
                             options::OPFOptions = OPFOptions())
    
    results = Dict{String, Any}()
    
    for formulation in formulations
        form_name = string(formulation)
        
        try
            result, perf = profile_solve(case_data, formulation, options)
            
            results[form_name] = Dict{String, Any}(
                "success" => result.success,
                "objective" => result.objective,
                "solve_time" => perf["solve_time"],
                "memory_mb" => perf["memory_mb"],
                "iterations" => result.iterations,
                "variables" => perf["variables"],
                "violation" => result.constraint_violation
            )
            
        catch e
            results[form_name] = Dict{String, Any}(
                "success" => false,
                "error" => string(e)
            )
        end
    end
    
    return results
end

# ============================================================================
# Validation Utilities
# ============================================================================

"""
Validate solution against physical constraints.
"""
function validate_solution(result::OPFResult{F}, tolerance::Float64 = 1e-6) where {F}
    
    if !result.success
        return Dict("valid" => false, "reason" => "Solve failed")
    end
    
    validation = Dict{String, Any}("valid" => true, "warnings" => String[])
    
    # Check constraint violation
    if result.constraint_violation > tolerance * 10
        validation["valid"] = false
        validation["reason"] = "Large constraint violation: $(result.constraint_violation)"
        return validation
    end
    
    # Check physical feasibility
    case_data = result.case_data
    
    # Power balance check
    if haskey(case_data, "busAC") && haskey(case_data, "genAC")
        bus = case_data["busAC"]
        gen = case_data["genAC"]
        
        total_gen = sum(gen[gen[:, 8] .> 0, 2])  # Active generators PG
        total_load = sum(bus[:, 3])  # PD
        
        power_balance_error = abs(total_gen - total_load) / total_load
        
        if power_balance_error > tolerance
            push!(validation["warnings"], "Power balance error: $(power_balance_error)")
        end
    end
    
    # Voltage magnitude checks (AC only)
    if F == ACFormulation && haskey(case_data, "busAC")
        bus = case_data["busAC"]
        vm = bus[:, 8]  # VM
        
        if any(vm .< 0.8) || any(vm .> 1.2)
            push!(validation["warnings"], "Voltage magnitudes outside reasonable range")
        end
    end
    
    return validation
end

# ============================================================================
# Export Utilities
# ============================================================================

"""
Export results to various formats.
"""
function export_results(result::OPFResult{F}, format::Symbol = :dict) where {F}
    
    if format == :dict
        return result.case_data
        
    elseif format == :summary
        return Dict{String, Any}(
            "formulation" => string(F),
            "success" => result.success,
            "objective" => result.objective,
            "solve_time" => result.solve_time,
            "iterations" => result.iterations,
            "constraint_violation" => result.constraint_violation
        )
        
    elseif format == :csv_buses
        if haskey(result.case_data, "busAC")
            bus = result.case_data["busAC"]
            headers = ["Bus", "Type", "PD", "QD", "GS", "BS", "Area", "VM", "VA", "BaseKV", "Zone", "VMin", "VMax"]
            return (headers, bus)
        else
            return nothing
        end
        
    else
        throw(ArgumentError("Unsupported export format: $format"))
    end
end
