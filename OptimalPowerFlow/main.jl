include("../data/case30.jl")
include("../solvers/types.jl")
include("../solvers/interior_point_method.jl")
include("../src/PowerFlow.jl")

using .PowerFlow
using Statistics  # Add this import for mean function
include("run_dcopf.jl")
include("run_acopf.jl")


jpc = case30()

# ROBUST DC OPF OPTIONS
opt_dc = Dict(
    "solver" => "interior_point_method",
    "OPF_MAX_IT" => 150,
    "OPF_VIOLATION" => 1e-5,
    "VERBOSE" => 1,
    # Bounds validation options
    "ANGLE_LIMIT_DEG" => 60.0,
    "AVG_LINE_FACTOR" => 2.0,
    "MIN_POWER_GAP" => 0.01,
    "BOUNDS_MARGIN" => 1e-6,
    "INIT_PROJECTION_WARN" => 0.1
)

# ROBUST AC OPF OPTIONS  
opt_ac = Dict(
    "solver" => "interior_point_method",
    "OPF_MAX_IT" => 200,
    "OPF_VIOLATION" => 5e-5,
    "VERBOSE" => 0,  # Reduce verbosity for comparison
    # Voltage bounds
    "DEFAULT_VMIN" => 0.95,
    "DEFAULT_VMAX" => 1.05,
    "MIN_VOLTAGE_GAP" => 0.05,
    # Generator bounds
    "DEFAULT_PMIN_FRACTION" => 0.0,
    "DEFAULT_PMAX_FRACTION" => 1.0,
    "DEFAULT_QMIN_FRACTION" => -0.5,
    "DEFAULT_QMAX_FRACTION" => 0.5,
    # Angle limits
    "ANGLE_LIMIT_DEG" => 60.0,
    # Other options
    "MIN_POWER_GAP" => 0.01,
    "BOUNDS_MARGIN" => 1e-6
)

println("Running DC OPF with robust bounds validation...")
jpc_dc = rundcopf(deepcopy(jpc), opt_dc)

println("\n" * "="^60)
println("COMPARING AC OPF IMPLEMENTATIONS")
println("="^60)

# Method 1: Original monolithic runacopf
println("\n📊 Method 1: Original runacopf function")
println("-"^40)
jpc_original = deepcopy(jpc)
time_original = @elapsed result_original = runacopf(jpc_original, opt_ac)
memory_original = @allocated runacopf(deepcopy(jpc), opt_ac)
