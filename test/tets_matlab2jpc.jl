using Base.Threads
using HyDistFlow
using PowerModels
using PyCall
using JuMP
using Ipopt
using Statistics
# include("../../data/output_case.jl")

# Convert a MATPOWER case file to a Julia PowerCase
# Option1 - Convert a MATPOWER case file to a Julia PowerCase
jpc = convert_matpower_case("E:/luosipeng/Package/matpower8.0/data/case30.m", "data/output_case.jl")
# Option2 - Use the predefined case data
# jpc = case_data()

opt = PowerFlow.options() # The initial settings
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 0;
# Extract islands from the JPC
jpc_list, isolated = PowerFlow.extract_islands(jpc)
n_islands = length(jpc_list)
println("Extracted $(n_islands) islands in total")

# Create results array
results_array = Vector{Any}(undef, n_islands)

println("Starting multi-threaded computation...")
t_start = time()

# Use multi-threading to compute power flow for each island
@threads for i in 1:n_islands
    results_array[i] = PowerFlow.runpf(jpc_list[i], opt)
end

t_end = time()
elapsed = t_end - t_start

# Construct results similar to @timed return format
results = (value=results_array, time=elapsed)

# Get voltage results for all buses
Vm_tspf = results[1][1].busAC[:,VM]
Va_tspf = results[1][1].busAC[:,VA]

## PowerModels PF results

nlp_solver = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "tol"=>1e-6,
    "print_level"=>0,
)
result = solve_ac_pf("E:/luosipeng/Package/matpower8.0/data/case30.m", nlp_solver)

bus_dict = result["solution"]["bus"]

# 获取所有节点编号并排序
bus_numbers = sort([parse(Int, k) for k in keys(bus_dict)])

# 提取电压幅值
Vm_powermodels = [bus_dict[string(bus_num)]["vm"] for bus_num in bus_numbers]
Va_powermodels = [bus_dict[string(bus_num)]["va"] for bus_num in bus_numbers]   

## pandapower PF results
pp = pyimport("pandapower")
pp_networks = pyimport("pandapower.networks")

py"""
import pandapower as pp
import pandapower.networks as pp_networks

net = pp_networks.case30()
pp.runpp(net)

# Access complete results
bus_voltages = net.res_bus[['vm_pu', 'va_degree']]
line_loading = net.res_line[['p_from_mw', 'q_from_mvar', 'loading_percent']]
generator_output = net.res_gen[['p_mw', 'q_mvar']]

Vm = bus_voltages['vm_pu'].values
Va = bus_voltages['va_degree'].values
"""
Vm_pandapower = py"Vm"
Va_pandapower = py"Va"

## comparison
Vm_error_with_powermodels = abs.(Vm_tspf - Vm_powermodels)
Vm_error_with_pandapower = abs.(Vm_tspf - Vm_pandapower)
Va_error_with_powermodels = abs.(Va_tspf - Va_powermodels*180/pi)
Va_error_with_pandapower = abs.(Va_tspf - Va_pandapower)

max_vm_error_with_powermodels = maximum(Vm_error_with_powermodels)
max_vm_error_with_pandapower = maximum(Vm_error_with_pandapower)
max_va_error_with_powermodels = maximum(Va_error_with_powermodels)
max_va_error_with_pandapower = maximum(Va_error_with_pandapower)

min_vm_error_with_powermodels = minimum(Vm_error_with_powermodels)
min_vm_error_with_pandapower = minimum(Vm_error_with_pandapower)
min_va_error_with_powermodels = minimum(Va_error_with_powermodels)
min_va_error_with_pandapower = minimum(Va_error_with_pandapower)

average_vm_error_with_powermodels = mean(Vm_error_with_powermodels)
average_vm_error_with_pandapower = mean(Vm_error_with_pandapower)
average_va_error_with_powermodels = mean(Va_error_with_powermodels)
average_va_error_with_pandapower = mean(Va_error_with_pandapower)

standard_deviation_vm_error_with_powermodels = std(Vm_error_with_powermodels)
standard_deviation_vm_error_with_pandapower = std(Vm_error_with_pandapower)
standard_deviation_va_error_with_powermodels = std(Va_error_with_powermodels)
standard_deviation_va_error_with_pandapower = std(Va_error_with_pandapower) 

