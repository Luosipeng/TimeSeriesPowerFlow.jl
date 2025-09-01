"""
   Main function to call the DC power flow function
    Input: case file
    Output: results of the power flow as a dictionary
    Example:
    bus, gen, branch = rundcpf(casefile)
"""
function rundcpf(mpc, opt::Dict{String})
   # Step 2.1: Define the data structures
    # Define named indices into bus, gen, branch matrices

    # Step 2.2: Convert the data into the required format
    baseMVA = mpc["baseMVA"];
    bus =  mpc["busDC"];
    gen = mpc["genDC"];
    branch = mpc["branchDC"];   
    load = mpc["loadDC"];
    pvarray = mpc["pv"]
    gpu_flag = opt["PF"]["GPU_ACCELERATION"];  ## use GPU acceleration?
    success = false;
    # convert the external data to internal data 
    (bus, gen, branch, load, pvarray, i2e) = PowerFlow.ext2int(bus, gen, branch, load, pvarray);
    ## get bus index lists of each type of bus
    (ref, p) = PowerFlow.dcbustypes(bus, gen);
    ## generator info
    on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
    gbus = gen[on, GEN_BUS]  # what buses are they at?
    # Step 2.3: Run the power flow
    ##-----  run the power flow  ----- 
    alg = opt["PF"]["PF_ALG"];
    its = 0;            ## total iterations
    ## initialize
    V0  = bus[:, VM]
    ## build admittance matrices
    (Ybus, Yf, Yt) = PowerFlow.makeYbus(baseMVA, bus, branch)
    repeat=1;
    while (repeat>0)
        ## function for computing V dependent complex bus power injections
        if alg == "NR"
            if gpu_flag == 1 
                V, success, iterations = newtondcpf_gpu(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, opt["PF"]["PF_DC_TOL"], opt["PF"]["PF_DC_MAX_IT"], opt["PF"]["NR_ALG"])
            else
                V, success, iterations = newtondcpf(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, opt["PF"]["PF_DC_TOL"], opt["PF"]["PF_DC_MAX_IT"], opt["PF"]["NR_ALG"]);
            end
            # V, success, iterations =  PowerFlow.adaptive_damped_newton(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, opt["PF"]["PF_DC_TOL"], opt["PF"]["PF_DC_MAX_IT"], opt["PF"]["NR_ALG"]);
            its += iterations;
        end
        bus, gen, branch = dcpfsoln(baseMVA, bus, gen, branch, load, Ybus, Yf, Yt, V, ref, p)
        repeat = 0;
    end
    bus, gen, branch, load, pvarray, areas=PowerFlow.int2ext(i2e, bus, gen, branch, load, pvarray);
    mpc["busDC"] = bus
    mpc["genDC"] = gen
    mpc["branchDC"] = branch
    mpc["loadDC"] = load
    mpc["pv"] = pvarray
    mpc["iterationsDC"] = its
    mpc["success"] = success
    return mpc
end