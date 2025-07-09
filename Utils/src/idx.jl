function idx_bus()
    # define bus types
    PQ      = 1;
    PV      = 2;
    REF     = 3;
    NONE    = 4;

    # define the indices
    BUS_I       = 1;    ## bus number (1 to 29997)
    BUS_TYPE    = 2;    ## bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
    PD          = 3;    ## Pd, real power demand (MW)
    QD          = 4;    ## Qd, reactive power demand (MVAr)
    GS          = 5;    ## Gs, shunt conductance (MW at V = 1.0 p.u.)
    BS          = 6;    ## Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
    BUS_AREA    = 7;    ## area number, 1-100
    VM          = 8;    ## Vm, voltage magnitude (p.u.)
    VA          = 9;    ## Va, voltage angle (degrees)
    BASE_KV     = 10;   ## baseKV, base voltage (kV)
    ZONE        = 11;   ## zone, loss zone (1-999)
    VMAX        = 12;   ## maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
    VMIN        = 13;   ## minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    LAM_P       = 14;   ## Lagrange multiplier on real power mismatch (u/MW)
    LAM_Q       = 15;   ## Lagrange multiplier on reactive power mismatch (u/MVAr)
    MU_VMAX     = 16;   ## Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    MU_VMIN     = 17;   ## Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
    PER_CONSUMER=18;    ## PER_CONSUMER, PER_CONSUMER

return PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER
end

#-------------------------------------------------
## define the indices
function idx_brch()
    F_BUS       = 1;    ## f, from bus number
    T_BUS       = 2;    ## t, to bus number
    BR_R        = 3;    ## r, resistance (p.u.)
    BR_X        = 4;    ## x, reactance (p.u.)
    BR_B        = 5;    ## b, total line charging susceptance (p.u.)
    RATE_A      = 6;    ## rateA, MVA rating A (long term rating)
    RATE_B      = 7;    ## rateB, MVA rating B (short term rating)
    RATE_C      = 8;    ## rateC, MVA rating C (emergency rating)
    TAP         = 9;    ## ratio, transformer off nominal turns ratio
    SHIFT       = 10;   ## angle, transformer phase shift angle (degrees)
    BR_STATUS   = 11;   ## initial branch status, 1 - in service, 0 - out of service
    ANGMIN      = 12;   ## minimum angle difference, angle(Vf) - angle(Vt) (degrees)
    ANGMAX      = 13;   ## maximum angle difference, angle(Vf) - angle(Vt) (degrees)
    DICTKEY     = 14;   ## dictionnary key (not in PTI format)
    ## included in power flow solution, not necessarily in input
    PF          = 15;   ## real power injected at "from" bus end (MW)       (not in PTI format)
    QF          = 16;   ## reactive power injected at "from" bus end (MVAr) (not in PTI format)
    PT          = 17;   ## real power injected at "to" bus end (MW)         (not in PTI format)
    QT          = 18;   ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    MU_SF       = 19;   ## Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
    MU_ST       = 20;   ## Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
    MU_ANGMIN   = 21;   ## Kuhn-Tucker multiplier lower angle difference limit (u/degree)
    MU_ANGMAX   = 22;   ## Kuhn-Tucker multiplier upper angle difference limit (u/degree)
    ##
    LAMBDA      = 23;   ## Lagrange multiplier on real power mismatch at "from" bus (u/MW)
    SW_TIME = 24;   ## time of last switch [s]
    RP_TIME = 25;
    BR_TYPE = 26;
    BR_AREA = 27;

    return F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA
end

#-------------------------------------------------
function idx_gen()
    ## define the indices
    GEN_BUS     = 1;    ## bus number
    PG          = 2;    ## Pg, real power output (MW)
    QG          = 3;    ## Qg, reactive power output (MVAr)
    QMAX        = 4;    ## Qmax, maximum reactive power output at Pmin (MVAr)
    QMIN        = 5;    ## Qmin, minimum reactive power output at Pmin (MVAr)
    VG          = 6;    ## Vg, voltage magnitude setpoint (p.u.)
    MBASE       = 7;    ## mBase, total MVA base of this machine, defaults to baseMVA
    GEN_STATUS  = 8;    ## status, 1 - machine in service, 0 - machine out of service
    PMAX        = 9;    ## Pmax, maximum real power output (MW)
    PMIN        = 10;   ## Pmin, minimum real power output (MW)
    PC1         = 11;   ## Pc1, lower real power output of PQ capability curve (MW)
    PC2         = 12;   ## Pc2, upper real power output of PQ capability curve (MW)
    QC1MIN      = 13;   ## Qc1min, minimum reactive power output at Pc1 (MVAr)
    QC1MAX      = 14;   ## Qc1max, maximum reactive power output at Pc1 (MVAr)
    QC2MIN      = 15;   ## Qc2min, minimum reactive power output at Pc2 (MVAr)
    QC2MAX      = 16;   ## Qc2max, maximum reactive power output at Pc2 (MVAr)
    RAMP_AGC    = 17;   ## ramp rate for load following/AGC (MW/min)
    RAMP_10     = 18;   ## ramp rate for 10 minute reserves (MW)
    RAMP_30     = 19;   ## ramp rate for 30 minute reserves (MW)
    RAMP_Q      = 20;   ## ramp rate for reactive power (2 sec timescale) (MVAr/min)
    APF         = 21;   ## area participation factor
    PW_LINEAR   = 22;   ## piecewise linear cost data, 2n pairs of real (MW) and positive real (u$/hr) (n must be less than 10)
    POLYNOMIAL = 23;   ## polynomial cost data, p+1 coefficients (p must be less than 10)
    MODEL       = 24;   ## generator model, 0 - model 1, 1 - model 2
    STARTUP     = 25;   ## startup cost in US dollars
    SHUTDOWN    = 26;   ## shutdown cost in US dollars
    NCOST       = 27;   ## number of cost coefficients for polynomial cost function, or number of PW pairs
    COST        = 28;   ## parameters defining total cost function f(Pg) for general model

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    MU_PMAX     = 29;   ## Kuhn-Tucker multiplier on upper Pg limit (u/MW)
    MU_PMIN     = 30;   ## Kuhn-Tucker multiplier on lower Pg limit (u/MW)
    MU_QMAX     = 31;   ## Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
    MU_QMIN     = 32;   ## Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

    GEN_AREA = 33;   ## area number, 1-100
    ## Note: When a generator's PQ capability curve is not simply a box and the
    ## upper Qg limit is binding, the multiplier on this constraint is split into
    ## it's P and Q components and combined with the appropriate MU_Pxxx and
    ## MU_Qxxx values. Likewise for the lower Q limits.
    return GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, GEN_AREA
end

#-------------------------------------------------
function idx_dcbus()
    # define bus types
    P      = 1;
    REF     = 2;
    NONE    = 3;

    # define the indices
    BUS_I       = 1;    ## bus number (1 to 29997)
    BUS_TYPE    = 2;    ## bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
    PD          = 3;    ## Pd, real power demand (MW)
    QD          = 4;    ## Qd, reactive power demand (MVAr)
    GS          = 5;    ## Gs, shunt conductance (MW at V = 1.0 p.u.)
    BS          = 6;    ## Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
    BUS_AREA    = 7;    ## area number, 1-100
    VM          = 8;    ## Vm, voltage magnitude (p.u.)
    VA          = 9;    ## Va, voltage angle (degrees)
    BASE_KV     = 10;   ## baseKV, base voltage (kV)
    ZONE        = 11;   ## zone, loss zone (1-999)
    VMAX        = 12;   ## maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
    VMIN        = 13;   ## minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    LAM_P       = 14;   ## Lagrange multiplier on real power mismatch (u/MW)
    LAM_Q       = 15;   ## Lagrange multiplier on reactive power mismatch (u/MVAr)
    MU_VMAX     = 16;   ## Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    MU_VMIN     = 17;   ## Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
    PER_CONSUMER = 18;    ## PER_CONSUMER, PER_CONSUMER
    return P, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN,PER_CONSUMER
end

function idx_ld()
    LOAD_I=1;
    LOAD_CND=2;
    LOAD_STATUS=3;
    LOAD_PD=4;
    LOAD_QD=5;
    LOADZ_PERCENT=6;
    LOADI_PERCENT=7;
    LOADP_PERCENT=8;


    return LOAD_I,LOAD_CND,LOAD_STATUS,LOAD_PD,LOAD_QD,LOADZ_PERCENT,LOADI_PERCENT,LOADP_PERCENT
end

function idx_hvcb()
    HVCB_ID=1;
    HVCB_FROM_ELEMENT=2;
    HVCB_TO_ELEMENT=3;
    HVCB_INSERVICE=4;
    HVCB_STATUS=5;

    return HVCB_ID,HVCB_FROM_ELEMENT,HVCB_TO_ELEMENT,HVCB_INSERVICE,HVCB_STATUS
end

function idx_microgrid()
    MG_ID=1;
    MG_CAPACITY=2;
    MG_PEAK_LOAD=3;
    MG_DURATION=4;
    MG_AREA=5;
    return MG_ID,MG_CAPACITY,MG_PEAK_LOAD,MG_DURATION,MG_AREA
end

function idx_pv()
    PV_ID=1;
    PV_BUS=2;
    PV_VOC=3;
    PV_VMPP=4;
    PV_ISC =5;
    PV_IMPP=6;
    PV_IRRADIANCE=7; # PV irradiance in W/m^2
    PV_AREA=8;
    PV_IN_SERVICE=9;
    return PV_ID,PV_BUS,PV_VOC,PV_VMPP,PV_ISC,PV_IMPP,PV_IRRADIANCE,PV_AREA,PV_IN_SERVICE
end

function idx_pv_acsystem()
    PV_AC_ID=1;
    PV_AC_BUS=2;
    PV_AC_VOC=3;
    PV_AC_VMPP=4;
    PV_AC_ISC =5;
    PV_AC_IMPP=6;
    PV_AC_IRRADIANCE=7; # PV irradiance in W/m^2
    PV_AC_INVERTER_EFFICIENCY=8; # Inverter efficiency, typically a value between 0 and 1
    PV_AC_INVERTER_MODE=9; 
    PV_AC_INVERTER_PAC=10; # Active power output of the inverter in MW
    PV_AC_INVERTER_QAC=11; # Reactive power output of the inverter in MVAr
    PV_AC_INVERTER_QAC_MAX=12; # Maximum reactive power output of the inverter in MVAr
    PV_AC_INVERTER_QAC_MIN=13; # Minimum reactive power output of the inverter in MVAr
    PV_AC_AREA=14;
    PV_AC_IN_SERVICE=15;
    return PV_AC_ID, PV_AC_BUS, PV_AC_VOC, PV_AC_VMPP, PV_AC_ISC, PV_AC_IMPP, PV_AC_IRRADIANCE, PV_AC_INVERTER_EFFICIENCY, PV_AC_INVERTER_MODE, PV_AC_INVERTER_PAC, PV_AC_INVERTER_QAC, PV_AC_INVERTER_QAC_MAX, PV_AC_INVERTER_QAC_MIN, PV_AC_AREA, PV_AC_IN_SERVICE
end

function idx_conv()
    CONV_ACBUS=1;
    CONV_DCBUS=2;
    CONV_INSERVICE=3;
    CONV_P_AC=4;
    CONV_Q_AC=5;
    CONV_P_DC=6;
    CONV_EFF =7;
    CONV_MODE=8; # 1-constant δs、Us, 2-constant Ps、Qs, 3-constant Ps、Us, 4-constant Udc、Qs, 5-constant Udc、Us, 6-Droop Udc、Constant Qs, 7-Droop Udc、Constant Us
    CONV_DROOP_KP = 9; # Droop gain, only for mode 6 and 7

    return CONV_ACBUS,CONV_DCBUS,CONV_INSERVICE,CONV_P_AC,CONV_Q_AC,CONV_P_DC,CONV_EFF, CONV_MODE, CONV_DROOP_KP
end
function idx_ess()
    ESS_BUS=1;
    ESS_POWER_CAPACITY=2;
    ESS_ENERGY_CAPACITY=3;
    ESS_SOC_INIT=4;
    ESS_SOC_MIN=5;
    ESS_SOC_MAX=6;
    ESS_EFFICIENCY = 7;
    ESS_STATUS=8;
    return ESS_BUS,ESS_POWER_CAPACITY,ESS_ENERGY_CAPACITY,ESS_SOC_INIT,ESS_SOC_MIN,ESS_SOC_MAX,ESS_EFFICIENCY, ESS_STATUS
end

function idx_ev()
    EV_ID=1;
    EV_CAPACITY=2;
    EV_FLEX_CAPACITY=3;
    EV_AREA=4;
    return EV_ID,EV_CAPACITY,EV_FLEX_CAPACITY,EV_AREA
end

function idx_jpc_3ph()
    RES_3PH_BUS = 1;    ## bus number (1 to 29997)
    RES_3PH_VM_A = 2;    ## Vm, voltage magnitude (p.u.)
    RES_3PH_VM_B = 3;    ## Vm, voltage magnitude (p.u.)
    RES_3PH_VM_C = 4;    ## Vm, voltage magnitude (p.u.)
    RES_3PH_VA_A = 5;    ## Va, voltage angle (degrees)
    RES_3PH_VA_B = 6;    ## Va, voltage angle (degrees)
    RES_3PH_VA_C = 7;    ## Va, voltage angle (degrees)
    RES_3PH_UNBALANCED = 8; ## unbalanced flag, 1 - unbalanced, 0 - balanced
    RES_3PH_PA_MW = 9;    ## real power injected at "from" bus end (MW)       (not in PTI format)
    RES_3PH_PB_MW = 10;   ## real power injected at "to" bus end (MW)         (not in PTI format)
    RES_3PH_PC_MW = 11;   ## real power injected at "to" bus end (MW)         (not in PTI format)
    RES_3PH_QA_MVAR = 12; ## reactive power injected at "from" bus end (MVAr) (not in PTI format)
    RES_3PH_QB_MVAR = 13; ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)
    RES_3PH_QC_MVAR = 14; ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)

    return RES_3PH_BUS, RES_3PH_VM_A, RES_3PH_VM_B, RES_3PH_VM_C, RES_3PH_VA_A, RES_3PH_VA_B, RES_3PH_VA_C, RES_3PH_UNBALANCED, RES_3PH_PA_MW, RES_3PH_PB_MW, RES_3PH_PC_MW, RES_3PH_QA_MVAR, RES_3PH_QB_MVAR, RES_3PH_QC_MVAR
end
