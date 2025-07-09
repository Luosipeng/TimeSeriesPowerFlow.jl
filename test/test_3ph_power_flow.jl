using  Dates
using  XLSX
using  DataFrames
using  PowerFlow
using  Base.Threads

file_path = "C:/Users/13733/Desktop/etap-main/parameters.xlsx"

case = PowerFlow.load_julia_power_data(file_path)

#拓扑处理

results, new_case = PowerFlow.topology_analysis(case, output_file="topology_results.xlsx")

jpc_3ph = PowerFlow.JuliaPowerCase2Jpc_3ph(new_case)