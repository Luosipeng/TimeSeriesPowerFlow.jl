using  Dates
using  XLSX
using  DataFrames
using  Base.Threads

using  PowerFlow
# file_path = joinpath(pwd(), "data", "etap_runpf_acdc.xlsx")
file_path = joinpath(pwd(), "data", "control_mode_test.xlsx")
# file_path = joinpath(pwd(), "data", "石桥F12草河F27交直流.xlsx")

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

jpc_list, isolated = PowerFlow.extract_islands_acdc(jpc)
n_islands = length(jpc_list)
println("共提取出 $(n_islands) 个孤岛")

# 创建结果数组
results_array = Vector{Any}(undef, n_islands)

println("开始多线程计算...")
t_start = time()

# 使用多线程计算每个孤岛的潮流
@threads for i in 1:n_islands
    results_array[i] = PowerFlow.runhpf(jpc_list[i], opt)
end

t_end = time()
elapsed = t_end - t_start

# 构造类似@timed返回的结果
results = (value=results_array, time=elapsed)


# # 获取所有节点的电压结果
voltage_results = PowerFlow.get_bus_voltage_results_acdc(results, new_case)

# # 比较结果与参考文件
# result_file = joinpath(pwd(), "data", "石桥F12草河F27交直流result.xlsx")
# result_file = "C:/Users/13733/Desktop/etap-main/result.xlsx"
# PowerFlow.analyze_voltage_results(results, case, result_file, output_dir="./analysis_results")

# println("计算完成，耗时: $(results.time) 秒")
# PowerFlow.process_result(results, isolated, "powerflow_report.txt")
