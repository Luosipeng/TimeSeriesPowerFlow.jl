using  Dates
using  XLSX
using  DataFrames
using  PowerFlow
using  Base.Threads
include(pwd()*"/data/case3_integrated_energy_systems.jl")

jpc = case3()
opt = PowerFlow.options() # The initial settings
opt["PF"]["NR_ALG"] = "\\";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 0;

jpc_results = PowerFlow.runpf(jpc, opt)
