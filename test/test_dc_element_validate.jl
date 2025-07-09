"""
    Main function for the AC power flow
"""

push!(LOAD_PATH, pwd()*"/src/")
include(pwd()*"/data/case118.jl")
include(pwd()*"/data/case33bw.jl")
include(pwd()*"/data/case_reliability.jl")

using PowerFlow
# using MATLAB
# using Base.Threads
using Plots

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;


file_path = "C:/Users/13733/Desktop/model_test/battery_test/battery_data.xlsx"

mpc, Dict_busdc, pv_curves = PowerFlow.excel2jpc(file_path)

mpc = PowerFlow.rundcpf(mpc, opt)

result_path = "C:/Users/13733/Desktop/model_test/battery_test/battery_result.xlsx"

PowerFlow.dc_element_validate("DCLUMPLOAD", "RatedKW", Dict_busdc, result_path, opt, mpc; connectedname=2)
