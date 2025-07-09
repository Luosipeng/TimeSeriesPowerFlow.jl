"""
    Main function for the AC power flow
"""

# Detect the current working operating system
if Sys.iswindows()
    # Add the path to the data folder
    push!(LOAD_PATH, pwd()*"\\src\\")
    include(pwd()*"\\data\\case118.jl")
else
    # Add the path to the data folder
    push!(LOAD_PATH, pwd()*"/src/")
    include(pwd()*"/data/case118.jl")
    using AppleAccelerate
end
# push!(LOAD_PATH, pwd()*"\\data\\");
using PowerFlow

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "gmres";
opt["PF"]["ENFORCE_Q_LIMS"]=0
mpc=case118()
@time mpc = PowerFlow.rundcpf(mpc, opt)