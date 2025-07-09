function case3ac()
    mpc = Dict{String, Any}();
    mpc["version"] = "2";
    mpc["baseMVA"] = 100.0;
    
    ## bus data
    #	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
    mpc["bus"] = [
        1  3  0.0    0.0        0.0  0.0  0.0  1.0  0.0  10.0     1.0  1.05  0.8
        2  1  0.036  0.0174356  0.0  0.0  0.0  1.0  0.0   0.3996  1.0  1.05  0.8
    ];
    
    ## generator data
    #	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    mpc["gen"] = [
        1.0  0.0  0.0  700.0  -20.0  1.0  100.0  1.0  1000.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
    ];
    
    ## branch data
    #	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
    mpc["branch"] = [
         1.0  2.0  9.98752  199.75  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0
    ];
    return mpc
end