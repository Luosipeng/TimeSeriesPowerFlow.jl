function runprepf(mpc, opt)
    #
    opt["PF"]["DC"] = 1
    # Run the power flow
    mpc = PowerFlow.runpf(mpc, opt)
    opt["PF"]["DC"] = 0
    return mpc
    
end