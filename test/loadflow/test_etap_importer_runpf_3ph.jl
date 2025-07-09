using  Dates
using  XLSX
using  DataFrames
using  PowerFlow
using  Base.Threads


# file_path = joinpath(pwd(), "data", "etap_runpf_acdc.xlsx")
file_path = "C:/Users/13733/Desktop/etap-main/parameters.xlsx"
# file_path = joinpath(pwd(), "data", "石桥F12草河F27交直流.xlsx")

case = PowerFlow.load_julia_power_data(file_path)

#拓扑处理

results, new_case = PowerFlow.topology_analysis(case, output_file="topology_results.xlsx")

# 查看结果
println("发现了 ", nrow(results["cycles"]), " 个环路")
println("网络被分为 ", length(unique(results["nodes"].Partition)), " 个分区")

jpc_3ph, gs_eg, bs_eg = PowerFlow.JuliaPowerCase2Jpc_3ph(new_case)

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "\\";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;

jpc_3ph = PowerFlow.runupf(case, jpc_3ph ,gs_eg,bs_eg ,opt)


