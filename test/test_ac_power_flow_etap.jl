push!(LOAD_PATH, pwd()*"/src/")
include(pwd()*"/data/case118.jl")
include(pwd()*"/data/case33bw.jl")
include(pwd()*"/data/case_reliability.jl")

using PowerFlow
# using MATLAB
using Base.Threads
using Plots

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
# opt["PF"]["DC_PREPROCESS"] = 1;
#ETAP ac system parameters file path
# file_path = "C:/Users/13733/Desktop/one_line_test/九龙变01/九龙变data.xlsx"
file_path = "C:/Users/13733/Desktop/system_test/ac_power_flow/ac_system_data.xlsx"
#Process the ac system data
mpc, dict_bus, node_mapping, pv_curves = PowerFlow.excel2jpc(file_path)

mpc = PowerFlow.runpf(mpc,opt)

etap_file_path = "C:/Users/13733/Desktop/system_test/ac_power_flow/ac_power_flow_result.xlsx"
# etap_file_path = "C:/Users/13733/Desktop/one_line_test/九龙变01/九龙变result.xlsx"

voltage_errors, angle_errors, avg_v_error, avg_a_error =PowerFlow.ac_power_flow_compared(mpc, dict_bus, node_mapping, etap_file_path)
