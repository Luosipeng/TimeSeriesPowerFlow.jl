"""
TimeDomainPowerFlow module for power system time-domain analysis

This module provides tools and functions for time-series power flow analysis,
including dynamic dispatch, voltage profile analysis, and renewable integration.
"""
module TimeDomainPowerFlow

    using Printf
    using SparseArrays
    using LinearAlgebra
    using PrettyTables
    using AMD
    using SuiteSparse
    using LinearOperators
    # using AlgebraicMultigrid
    using XLSX
    using DataFrames
    using StatsBase
    using Graphs
    using DataStructures  
    using Dates
    using Base.Threads
    using Plots
    using JuMP
    using Ipopt
    
    # ...other files in Utils directory...
    include("../../PowerFlow/src/PowerFlow.jl")

    using .PowerFlow

    const Utils = PowerFlow.Utils
    const ComponentModel = PowerFlow.Utils.ComponentModel
    using .PowerFlow.Utils: load_julia_power_data, JuliaPowerCase,JPC,topology_analysis, extract_islands_acdc,JuliaPowerCase2Jpc,JuliaPowerCase2Jpc_3ph
    using .PowerFlow.ComponentModel: Storage
    using .PowerFlow.Utils: idx_brch, idx_bus, idx_gen, idx_dcbus, idx_ld, idx_hvcb, idx_microgrid, idx_pv,idx_conv, idx_ess, idx_jpc_3ph,idx_ev, idx_pv_acsystem
    using .PowerFlow: options
    

    include(joinpath(@__DIR__,"renumber_hybrid_system.jl"))
    include(joinpath(@__DIR__,"run_dynamic_dispatch.jl"))
    include(joinpath(@__DIR__,"runtdpf.jl"))
    include(joinpath(@__DIR__,"run_single_day.jl"))
    include(joinpath(@__DIR__,"plot_voltage_time_series.jl"))
    include(joinpath(@__DIR__,"plot_PD_time_series.jl"))
    include(joinpath(@__DIR__,"record_voltage_violation.jl"))
    include(joinpath(@__DIR__,"plot_losses_time_series.jl"))
    include(joinpath(@__DIR__,"plot_flow_violations.jl"))

    # ... export other indexes defined in idx.jl ...
    # Export bus indices
    const PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, 
        LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER = idx_bus()
    export PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, 
         LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER
  
    # Export branch indices
    const F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX, 
        DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA = idx_brch()
    export F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX, 
         DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA
  
    # Export generator indices
    const GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, 
        RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST, 
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, GEN_AREA = idx_gen()
    export GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, 
         RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST, 
         MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, GEN_AREA
  
    # Export DC bus indices
    const P, DC_REF, DC_NONE, DC_BUS_I, DC_BUS_TYPE, DC_PD, DC_QD, DC_GS, DC_BS, DC_BUS_AREA, DC_VM, DC_VA, 
        DC_BASE_KV, DC_ZONE, DC_VMAX, DC_VMIN, DC_LAM_P, DC_LAM_Q, DC_MU_VMAX, DC_MU_VMIN, DC_PER_CONSUMER = idx_dcbus()
    export P, DC_REF, DC_NONE, DC_BUS_I, DC_BUS_TYPE, DC_PD, DC_QD, DC_GS, DC_BS, DC_BUS_AREA, DC_VM, DC_VA, 
         DC_BASE_KV, DC_ZONE, DC_VMAX, DC_VMIN, DC_LAM_P, DC_LAM_Q, DC_MU_VMAX, DC_MU_VMIN, DC_PER_CONSUMER
  
    # Export load indices
    const LOAD_I, LOAD_CND, LOAD_STATUS, LOAD_PD, LOAD_QD, LOADZ_PERCENT, LOADI_PERCENT, LOADP_PERCENT = idx_ld()
    export LOAD_I, LOAD_CND, LOAD_STATUS, LOAD_PD, LOAD_QD, LOADZ_PERCENT, LOADI_PERCENT, LOADP_PERCENT
  
    # Export HVCB indices
    const HVCB_ID, HVCB_FROM_ELEMENT, HVCB_TO_ELEMENT, HVCB_INSERVICE, HVCB_STATUS = idx_hvcb()
    export HVCB_ID, HVCB_FROM_ELEMENT, HVCB_TO_ELEMENT, HVCB_INSERVICE, HVCB_STATUS
  
    # Export microgrid indices
    const MG_ID, MG_CAPACITY, MG_PEAK_LOAD, MG_DURATION, MG_AREA = idx_microgrid()
    export MG_ID, MG_CAPACITY, MG_PEAK_LOAD, MG_DURATION, MG_AREA
  
    # Export PV indices
    const PV_ID,PV_BUS,PV_VOC,PV_VMPP,PV_ISC,PV_IMPP,PV_IRRADIANCE,PV_AREA,PV_IN_SERVICE = idx_pv()
    export PV_ID,PV_BUS,PV_VOC,PV_VMPP,PV_ISC,PV_IMPP,PV_IRRADIANCE,PV_AREA,PV_IN_SERVICE

    # Export PV AC system indices
    const PV_AC_ID, PV_AC_BUS, PV_AC_VOC, PV_AC_VMPP, PV_AC_ISC, PV_AC_IMPP, PV_AC_IRRADIANCE, PV_AC_INVERTER_EFFICIENCY, PV_AC_INVERTER_MODE, PV_AC_INVERTER_PAC, PV_AC_INVERTER_QAC, PV_AC_INVERTER_QAC_MAX, PV_AC_INVERTER_QAC_MIN, PV_AC_AREA, PV_AC_IN_SERVICE = idx_pv_acsystem()
    export PV_AC_ID, PV_AC_BUS, PV_AC_VOC, PV_AC_VMPP, PV_AC_ISC, PV_AC_IMPP, PV_AC_IRRADIANCE, PV_AC_INVERTER_EFFICIENCY, PV_AC_INVERTER_MODE, PV_AC_INVERTER_PAC, PV_AC_INVERTER_QAC, PV_AC_INVERTER_QAC_MAX, PV_AC_INVERTER_QAC_MIN, PV_AC_AREA, PV_AC_IN_SERVICE


    # Export Converter indices
    const CONV_ACBUS,CONV_DCBUS,CONV_INSERVICE,CONV_P_AC,CONV_Q_AC,CONV_P_DC,CONV_EFF,CONV_MODE,CONV_DROOP_KP = idx_conv()
    export CONV_ACBUS,CONV_DCBUS,CONV_INSERVICE,CONV_P_AC,CONV_Q_AC,CONV_P_DC,CONV_EFF,CONV_MODE,CONV_DROOP_KP
  
    # Export ESS indices
    const ESS_BUS,ESS_POWER_CAPACITY,ESS_ENERGY_CAPACITY,ESS_SOC_INIT,ESS_SOC_MIN,ESS_SOC_MAX,ESS_EFFICIENCY, ESS_STATUS = idx_ess()
    export ESS_BUS,ESS_POWER_CAPACITY,ESS_ENERGY_CAPACITY,ESS_SOC_INIT,ESS_SOC_MIN,ESS_SOC_MAX,ESS_EFFICIENCY, ESS_STATUS

    # Export JPC_3ph indices
    const RES_3PH_BUS, RES_3PH_VM_A, RES_3PH_VM_B, RES_3PH_VM_C, RES_3PH_VA_A, RES_3PH_VA_B, RES_3PH_VA_C, RES_3PH_UNBALANCED, RES_3PH_PA_MW, RES_3PH_PB_MW, RES_3PH_PC_MW, RES_3PH_QA_MVAR, RES_3PH_QB_MVAR, RES_3PH_QC_MVAR = idx_jpc_3ph()
    export RES_3PH_BUS, RES_3PH_VM_A, RES_3PH_VM_B, RES_3PH_VM_C, RES_3PH_VA_A, RES_3PH_VA_B, RES_3PH_VA_C, RES_3PH_UNBALANCED, RES_3PH_PA_MW, RES_3PH_PB_MW, RES_3PH_PC_MW, RES_3PH_QA_MVAR, RES_3PH_QB_MVAR, RES_3PH_QC_MVAR

    # Export EV indices
    const EV_ID, EV_CAPACITY, EV_FLEX_CAPACITY, EV_AREA = idx_ev()
    export EV_ID, EV_CAPACITY, EV_FLEX_CAPACITY, EV_AREA

    # export  read_load_data,read_price_data,read_irradiance_data
    export runtdpf, run_dynamic_dispatch, run_single_day,plot_voltage_time_series,plot_PD_time_series,create_time_series_loads
    export record_voltage_violation, plot_losses_time_series, plot_flow_violations

    export topology_analysis, options
    export load_julia_power_data
    export create_time_series_storage_profile
    
end
