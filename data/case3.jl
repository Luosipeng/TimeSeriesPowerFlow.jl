function case3()
    mpc = Dict{String, Any}();
    mpc["version"] = "2";
    mpc["baseMVA"] = 10.0;
    #TODO:
    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    mpc["busAC"] = [
        1	3	0	0	0	0	1	1.0 	0       0.4	1	1.1	0.9;
	    2	1	0.07	0.0554	0	0	1	1.0 	0       0.4	1	1.1	0.9;
        ];
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    mpc["busDC"] = [
        1	1	0.3111 0.0	0   0	1	1.0 0 0.824	1	1.1	0.9;
        2	1	0   0.0	0   0	1	1.0 0 0.824	1	1.1	0.9;
        3	2	0   0.0	0   0	1	1.0 0 0.824	1	1.1	0.9;
    ];
    ## generator data
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q_apf apf
    mpc["genAC"] = [
        1	0.0     0.0     300	-300	1.0 	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
    ];

    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q_apf apf
    mpc["genDC"] = [
        3	0.0     0.0     300	-300	1.0 	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
    ];

    ## branch data
    # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax MVSC1 MVSC2 BRANCHMODE ETCR ETCI PHI
    mpc["branchAC"] = [
        1	2	0.5	0.5	0	250	150	150	0	0	1	-360	360;
    ];
    # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax MVSC1 MVSC2 BRANCHMODE ETCR ETCI PHI
    mpc["branchDC"] = [
        1	2	5.8912	0	0	250	150	150	0	0	1	-360	360	;
        2	3	0.3594	0	0	250	150	150	0	0	1	-360	360	;
    ];
    return mpc
end