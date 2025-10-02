using  Dates
using  XLSX
using  DataFrames
using  Base.Threads
using  HyDistFlow
file_path = joinpath(pwd(), "data", "test_case2.xlsx")

case = load_julia_power_data(file_path)

# Topology processing
results, new_case = topology_analysis(case, output_file="topology_results.xlsx")

# View results
println("Found ", nrow(results["cycles"]), " cycles")
println("Network is divided into ", length(unique(results["nodes"].Partition)), " partitions")

# empty!(new_case.storageetap)
# new_case.converters[3].control_mode = "Droop_Udc_Us"
# new_case.converters[2].control_mode = "Droop_Udc_Us"
# new_case.converters[1].control_mode = "Droop_Udc_Qs"

jpc = JuliaPowerCase2Jpc(new_case)

opt = options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab_gpu";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;
opt["PF"]["GPU_ACCELERATION"] = 1; 

jpc_list, isolated = extract_islands_acdc(jpc)
n_islands = length(jpc_list)
println("Extracted $(n_islands) islands in total")

# Create results array
results_array = Vector{Any}(undef, n_islands)

println("Starting multi-threaded calculation...")
t_start = time()

# Use multi-threading to calculate power flow for each island
@threads for i in 1:n_islands
    results_array[i] = runhpf(jpc_list[i], opt)
end

t_end = time()
elapsed = t_end - t_start

# Construct result similar to @timed return
results = (value=results_array, time=elapsed)

# # Get voltage results for all nodes
voltage_results = get_bus_voltage_results_acdc(results, new_case)

# # Compare results with reference file
result_file = joinpath(pwd(), "data", "test_case2_result.xlsx")
analyze_voltage_results(results, case, result_file, output_dir="./analysis_results", save_pdf=false)

# println("Calculation completed, time elapsed: $(results.time) seconds")
# process_result(results, isolated, "PowerFlow_report.txt")
