function case3()
    jpc = Dict{String, Any}();
    jpc["version"] = "2";
    jpc["baseMVA"] = 100.0;

    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    jpc["busAC"] = [
        1	3	100	50	0	0	1	1.0 	0       500	1	1.1	0.9;
	    2	1	0	0	0	0	1	1.0 	0       220	1	1.1	0.9;
        3	1	150	75	0	0	1	1.0 	0       110	1	1.1	0.9;
        ];
    ## generator data
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q_apf apf
    jpc["genAC"] = [
        1	0.0     0.0     5000	-5000	1.0 	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
    ];

    ## branch data
    # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax MVSC1 MVSC2 BRANCHMODE ETCR ETCI PHI
    jpc["branchAC"] = [
        1	2	0.05	0.1	 0	80	150	250	0	0	1	-360	360;
        2	3	0.05	0.1	 0	80	150	250	0	0	1	-360	360;
    ];

    ## load data
    # LOAD_I bus_i status Pd Qd load_z_percent load_i_percent load_i_percent
    jpc["loadAC"] = [
        1	1	1	100	50	1.0	0.0	0.0;
        2	3	1	150	75	1.0	0.0	0.0;
    ];

    # PV as an empty array, as we are not using PV generators in this case
    jpc["pv"] =   Array{Float64}(undef, 0, 8);

    return jpc
end