"""
    Main function for the AC power flow
"""

push!(LOAD_PATH, pwd()*"/src/")
# include(pwd()*"/data/case118.jl")
# include(pwd()*"/data/case33bw.jl")
# include(pwd()*"/data/case_reliability.jl")

using PowerFlow
# using MATLAB
# using Base.Threads
using Plots

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;


file_path = joinpath(pwd(), "data", "ac_dc_power_flow_data.xlsx")


mpc, dict_bus, node_mapping, pv_curves, Dict_busdc = PowerFlow.excel2jpc(file_path)

mpc = PowerFlow.process_inverter(mpc)

mpc=PowerFlow.runhpf(mpc,opt)

result_path = joinpath(pwd(), "data", "ac_dc_power_flow_result.xlsx")

PowerFlow.acdc_power_flow_compared(mpc, dict_bus,Dict_busdc,node_mapping,result_path)
