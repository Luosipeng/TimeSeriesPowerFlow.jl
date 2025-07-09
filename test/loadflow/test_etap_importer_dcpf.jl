using  Dates
using  XLSX
using  DataFrames
using  Base.Threads
using  PowerFlow

# file_path = joinpath(pwd(), "data", "etap_runpf.xlsx")
# file_path = "C:/Users/13733/Desktop/PV测试.xlsx"
# file_path = "C:/Users/13733/Desktop/etap-main/parameters.xlsx"
file_path = joinpath(pwd(), "data", "石桥F12草河F27交直流.xlsx")

case = PowerFlow.load_julia_power_data(file_path)

#拓扑处理

results, new_case = PowerFlow.topology_analysis(case, output_file="topology_results.xlsx")

# 查看结果
println("发现了 ", nrow(results["cycles"]), " 个环路")
println("网络被分为 ", length(unique(results["nodes"].Partition)), " 个分区")

jpc = PowerFlow.JuliaPowerCase2Jpc(new_case)

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "\\";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;

# jpc = PowerFlow.runhpf(jpc, opt)
jpc = PowerFlow.rundcpf(jpc, opt)
println("潮流计算完成")

# 比较结果与参考文件
# result_file = joinpath(pwd(), "data", "etap_runpf_result.xlsx")
# result_file = "C:/Users/13733/Desktop/石桥F12草河F27result.xlsx"
# result_file = "C:/Users/13733/Desktop/石桥F12草河F27连接result.xlsx"
# PowerFlow.analyze_voltage_results(results, case, result_file, output_dir="./analysis_results")

# println("计算完成，耗时: $(results.time) 秒")
# PowerFlow.process_result(results, isolated, "powerflow_report.txt")
