"""
    Main function for the AC power flow
"""

push!(LOAD_PATH, pwd()*"/src/")
# include(pwd()*"/data/case118.jl")
include(pwd()*"/data/case33mg_new.jl")
include(pwd()*"/data/case69_new.jl")
include(pwd()*"/data/case56_new.jl")
include(pwd()*"/data/case33mg_new_dist.jl")
include(pwd()*"/data/case69_new_dist.jl")
# include(pwd()*"/data/case_reliability.jl")

using PowerFlow
# using MATLAB
# using Base.Threads
using Plots

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;

# mpc=case33mg_new()
# mpc=case33mg_new_dist()
# mpc=case69_new()
mpc=case69_new_dist()
# mpc=case56_new()

# mpc = PowerFlow.runpf(mpc,opt)
# file_path = joinpath(pwd(), "data", "lindistflow_parameters_2.xlsx")

# mpc, dict_bus, node_mapping, pv_curves, Dict_busdc = PowerFlow.excel2jpc(file_path)

lindist_result = PowerFlow.run_lindistflow(mpc)

mpc = PowerFlow.process_inverter(lindist_result)
# mpc["genAC"][1, 6] = 1.00935
mpc["busAC"][:,8].= 1.0
mpc["busAC"][:, 9] .= 0.0
@time mpc=PowerFlow.runhpf(mpc,opt)

linVac = lindist_result["busAC"][:,8]
linVdc = lindist_result["busDC"][:,8]
Vac = mpc["busAC"][:,8]
Vdc = mpc["busDC"][:,8]

plot(Vac, label = "AC Power Flow", color = :blue, title = "Voltage Comparison", xlabel = "Bus Number", ylabel = "Voltage (p.u.)")
plot!(linVac, label = "Lindistflow", color = :red, title = "Voltage Comparison", xlabel = "Bus Number", ylabel = "Voltage (p.u.)")
# plot(linVdc, label = "DC Power Flow", color = :green, title = "Voltage Comparison", xlabel = "Bus Number", ylabel = "Voltage (p.u.)")
# plot!(Vdc, label = "AC Power Flow", color = :blue, title = "Voltage Comparison", xlabel = "Bus Number", ylabel = "Voltage (p.u.)")
# result_path = joinpath(pwd(), "data", "lindistflow_result.xlsx")
# PowerFlow.acdc_power_flow_compared(mpc, dict_bus,Dict_busdc,node_mapping,result_path)