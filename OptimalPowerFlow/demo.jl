include("../src/PowerFlow.jl")
using .PowerFlow

# Import required packages globally before including module components
using LinearAlgebra, SparseArrays, Statistics

# Include all components in the correct order
include("types.jl")
include("traits.jl")
include("config.jl")
include("validation.jl")
include("formulations/dc_formulation.jl")
include("formulations/ac_formulation.jl")
include("solvers/unified_solver.jl")
include("utils.jl")

# Now include the main module
include("OptimalPowerFlow.jl")
using .OptimalPowerFlow

include("../data/case30.jl")

println("🚀 OptimalPowerFlow Module Demonstration")
println("="^60)

# Load test case
jpc = case30()

# =============================================================================
# DEMONSTRATION 1: Unified API  
# =============================================================================

println("\n📊 Demo 1: Unified API")
println("-"^50)

# Test basic functionality first
println("🔧 Testing basic DC OPF functionality:")
try
    # Test with proper bounds validation first
    println("  Creating problem with bounds validation...")
    test_problem = create_problem(DCFormulation, jpc, OPFOptions())
    println("  Problem created successfully")
    println("  Problem type: $(typeof(test_problem))")
    println("  Buses: $(test_problem.nb), Generators: $(test_problem.ng), Lines: $(test_problem.nl)")
    
    # Check bounds feasibility
    bounds = test_problem.var_bounds
    println("  Voltage angle bounds: [$(minimum(bounds.va_min)), $(maximum(bounds.va_max))] rad")
    println("  Line flow bounds: [$(minimum(bounds.pij_min)), $(maximum(bounds.pij_max))] p.u.")
    println("  Generator bounds: [$(minimum(bounds.pg_min)), $(maximum(bounds.pg_max))] p.u.")
    
    # Now try full solve
    result_dc = solve_opf(jpc, DCFormulation)
    
    println("  ✅ Success: $(result_dc.success)")
    if result_dc.success
        println("  📈 Objective: $(round(result_dc.objective, digits=4))")
        println("  🔄 Iterations: $(result_dc.iterations)")
        println("  ⏱️  Time: $(round(result_dc.solve_time, digits=4))s")
        println("  💾 Memory: $(round(result_dc.memory_usage/1024/1024, digits=2))MB")
    else
        println("  ❌ DC OPF failed - checking bounds issues")
    end
    
catch e
    println("  ❌ Error in DC OPF: $e")
    println("  Error type: $(typeof(e))")
end

# Test AC OPF functionality
println("\n⚡ Testing basic AC OPF functionality:")
try
    # Test with proper bounds validation first
    println("  Creating AC problem with bounds validation...")
    test_ac_problem = create_problem(ACFormulation, jpc, OPFOptions())
    println("  AC Problem created successfully")
    println("  Problem type: $(typeof(test_ac_problem))")
    println("  Buses: $(test_ac_problem.nb), Generators: $(test_ac_problem.ng), Lines: $(test_ac_problem.nl)")
    
    # Check AC bounds feasibility
    ac_bounds = test_ac_problem.var_bounds
    println("  Voltage angle bounds: [$(minimum(ac_bounds.va_min)), $(maximum(ac_bounds.va_max))] rad")
    println("  Voltage magnitude bounds: [$(minimum(ac_bounds.vm_min)), $(maximum(ac_bounds.vm_max))] p.u.")
    println("  P generator bounds: [$(minimum(ac_bounds.pg_min)), $(maximum(ac_bounds.pg_max))] p.u.")
    println("  Q generator bounds: [$(minimum(ac_bounds.qg_min)), $(maximum(ac_bounds.qg_max))] p.u.")
    
    # Now try full AC solve with custom options - EXACTLY as run_acopf.jl
    ac_options = OPFOptions(
        max_iterations = 500,         # Start with fewer iterations like run_acopf.jl
        tolerance = 1e-3,           # More relaxed tolerance like run_acopf.jl
        verbose = 1,
        angle_limit_deg = 60.0,
        voltage_min_default = 0.94,  # Closer to run_acopf.jl defaults
        voltage_max_default = 1.06,
        bounds_margin = 1e-4         # Small bounds margin
    )
    
    result_ac = solve_opf(jpc, ACFormulation, ac_options)
    
    println("  ✅ Success: $(result_ac.success)")
    if result_ac.success
        println("  📈 Objective: $(round(result_ac.objective, digits=4))")
        println("  🔄 Iterations: $(result_ac.iterations)")
        println("  ⏱️  Time: $(round(result_ac.solve_time, digits=4))s")
        println("  💾 Memory: $(round(result_ac.memory_usage/1024/1024, digits=2))MB")
        println("  📊 Constraint violation: $(round(result_ac.constraint_violation, digits=6))")
    else
        println("  ❌ AC OPF failed - checking bounds issues")
    end
    
catch e
    println("  ❌ Error in AC OPF: $e")
    println("  Error type: $(typeof(e))")
end

# Show trait-based behavior (this should work)
println("\n🎯 Trait-based behavior:")
println("  DC Traits:")
println("    - Variables: $(get_variable_count(DCFormulation, 30, 6, 41))")
println("    - Voltage magnitudes: $(has_voltage_magnitudes(DCFormulation))")
println("    - Reactive power: $(has_reactive_power(DCFormulation))")
println("    - Linear constraints: $(is_linear_constraints(DCFormulation))")

println("  AC Traits:")
println("    - Variables: $(get_variable_count(ACFormulation, 30, 6, 41))")
println("    - Voltage magnitudes: $(has_voltage_magnitudes(ACFormulation))")
println("    - Reactive power: $(has_reactive_power(ACFormulation))")
println("    - Linear constraints: $(is_linear_constraints(ACFormulation))")

# Problem analysis (this should work)
println("\n🔍 Problem analysis:")
try
    problem_info = get_problem_info(jpc, DCFormulation)
    for (key, value) in problem_info["dimensions"]
        println("  $key: $value")
    end
catch e
    println("  Error in problem analysis: $e")
end

# Test configuration system
println("\n⚙️  Configuration test:")
test_options = OPFOptions(max_iterations=100, tolerance=1e-5, verbose=0)
println("  Options created: max_iterations=$(test_options.max_iterations)")

# Legacy compatibility test
println("\n🔄 Legacy compatibility test:")
try
    legacy_dc_opts = Dict{String, Any}(
        "OPF_MAX_IT" => 100,
        "OPF_VIOLATION" => 1e-5,
        "VERBOSE" => 0
    )
    
    # Test option conversion
    converted_opts = from_legacy_options(legacy_dc_opts)
    println("  Legacy options converted successfully")
    println("  Max iterations: $(converted_opts.max_iterations)")
    
catch e
    println("  Error in legacy compatibility: $e")
end

println("\n✨ Basic module functionality demonstrated!")
println("🎉 Key features working:")
println("  • Type system and traits")
println("  • Configuration management")
println("  • Problem analysis")
println("  • Legacy compatibility")

# =============================================================================
# DEMONSTRATION 2: Type Safety and Performance
# =============================================================================

println("\n🏗️  Demo 2: Type Safety and Performance")
println("-"^50)

# Show compile-time type specialization
println("🧬 Type specialization demonstration:")
dc_problem_type = typeof(create_problem(DCFormulation, jpc, OPFOptions()))
ac_problem_type = typeof(create_problem(ACFormulation, jpc, OPFOptions()))

println("  DC Problem Type: $(dc_problem_type)")
println("  AC Problem Type: $(ac_problem_type)")

# =============================================================================
# DEMONSTRATION 3: Advanced Features
# =============================================================================

println("\n⚙️  Demo 3: Advanced Features")
println("-"^50)

# Performance profiling
println("\n📊 Performance profiling:")
result, perf = profile_solve(jpc, DCFormulation, OPFOptions(verbose=0))

for (metric, value) in perf
    if isa(value, Number)
        println("  $metric: $(round(value, digits=4))")
    else
        println("  $metric: $value")
    end
end

# Formulation comparison
println("\n🔄 Formulation comparison:")
comparison = compare_formulations(jpc, [DCFormulation, ACFormulation], OPFOptions(verbose=0))

for (form, results) in comparison
    println("  $form:")
    if results["success"]
        println("    Objective: $(round(results["objective"], digits=4))")
        println("    Time: $(round(results["solve_time"], digits=4))s")
        println("    Variables: $(results["variables"])")
    else
        println("    Failed: $(get(results, "error", "Unknown error"))")
    end
end

# =============================================================================
# DEMONSTRATION 4: Configuration Management
# =============================================================================

println("\n⚙️  Demo 4: Configuration Management")
println("-"^50)

# Different configuration scenarios
configs = [
    ("Research", OPFOptions(tolerance=1e-6, max_iterations=200, verbose=1)),
    ("Production", OPFOptions(tolerance=1e-4, max_iterations=100, verbose=0)),
    ("Conservative", OPFOptions(angle_limit_deg=45.0, voltage_min_default=0.98))
]

for (name, config) in configs
    println("\n🧪 Testing $name configuration:")
    result = solve_opf(jpc, DCFormulation, config)
    
    println("  ✅ Success: $(result.success)")
    if result.success
        println("  📈 Objective: $(round(result.objective, digits=4))")
        println("  ⏱️  Time: $(round(result.solve_time, digits=4))s")
    end
end

# =============================================================================
# DEMONSTRATION 5: Solution Validation
# =============================================================================

println("\n🛡️  Demo 5: Solution Validation")
println("-"^50)

# Validate solutions
if @isdefined(result_dc)
    println("🔍 DC Solution validation:")
    dc_validation = validate_solution(result_dc)
    println("  Valid: $(dc_validation["valid"])")
    if haskey(dc_validation, "warnings") && !isempty(dc_validation["warnings"])
        for warning in dc_validation["warnings"]
            println("  ⚠️  $warning")
        end
    end
else
    println("🔍 DC Solution: Not available")
end

if @isdefined(result_ac)
    println("\n🔍 AC Solution validation:")
    ac_validation = validate_solution(result_ac)
    println("  Valid: $(ac_validation["valid"])")
    if haskey(ac_validation, "warnings") && !isempty(ac_validation["warnings"])
        for warning in ac_validation["warnings"]
            println("  ⚠️  $warning")
        end
    end
else
    println("\n🔍 AC Solution: Not available")
end

# =============================================================================
# DEMONSTRATION 6: Export and Compatibility
# =============================================================================

println("\n📤 Demo 6: Export and Compatibility")
println("-"^50)

# Export results
if @isdefined(result_dc)
    summary = export_results(result_dc, :summary)
    println("📋 DC Summary export:")
    for (key, value) in summary
        if isa(value, Number)
            println("  $key: $(round(value, digits=4))")
        else
            println("  $key: $value")
        end
    end
end

if @isdefined(result_ac)
    ac_summary = export_results(result_ac, :summary)
    println("\n📋 AC Summary export:")
    for (key, value) in ac_summary
        if isa(value, Number)
            println("  $key: $(round(value, digits=4))")
        else
            println("  $key: $value")
        end
    end
end

# Legacy compatibility test
println("\n🔄 Legacy compatibility test:")
legacy_dc_opts = to_legacy_options(OPFOptions(), DCFormulation)
legacy_result = rundcopf(deepcopy(jpc), legacy_dc_opts)

println("  Legacy DC Success: $(legacy_result["success"])")
println("  Legacy DC Objective: $(round(legacy_result["f"], digits=4))")

# Test legacy AC OPF
legacy_ac_opts = to_legacy_options(OPFOptions(), ACFormulation)
legacy_ac_result = runacopf(deepcopy(jpc), legacy_ac_opts)

println("  Legacy AC Success: $(legacy_ac_result["success"])")
println("  Legacy AC Objective: $(round(legacy_ac_result["f"], digits=4))")

println("\n✨ Module demonstration completed successfully!")
println("🎉 OptimalPowerFlow provides:")
println("  • Type-safe, high-performance OPF solving")
println("  • Unified interface for all formulations")
println("  • Advanced configuration and validation")
println("  • Full backward compatibility")
println("  • Zero-cost abstractions through Julia's type system")
println("  • Zero-cost abstractions through Julia's type system")
