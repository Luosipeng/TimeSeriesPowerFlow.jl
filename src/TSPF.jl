"""
TSPF

Time Series Power Flow module for power system analysis.
"""
module TSPF
    using TimeSeriesPowerFlow
    using XLSX
    using DataFrames
    using Base.Threads
    using Plots
    using DataStructures
    using Dates

    const Utils = TimeSeriesPowerFlow.Utils
    const ComponentModel = TimeSeriesPowerFlow.ComponentModel
    const PowerFlow = TimeSeriesPowerFlow.PowerFlow

    using Utils: load_julia_power_data, JuliaPowerCase,JPC,topology_analysis, extract_islands_acdc,JuliaPowerCase2Jpc,JuliaPowerCase2Jpc_3ph,extract_islands_acdc,extract_islands
    using TimeSeriesPowerFlow: read_load_data, read_price_data, read_irradiance_data,create_time_series_loads,create_time_series_prices,create_time_series_irradiance,extract_load_matrix_by_islands
    using Utils:idx_brch, idx_bus, idx_gen, idx_dcbus, idx_ld, idx_hvcb, idx_microgrid, idx_pv,idx_conv, idx_ess, idx_jpc_3ph,idx_ev, idx_pv_acsystem
    using PowerFlow: options,runhpf,runpf,rundcpf,get_bus_voltage_results_acdc,analyze_voltage_results,process_result
    using ComponentModel: AbstractComponent, Bus, Line, LineDC, Switch, BusDC, HighVoltageCircuitBreaker,Transformer2W, Transformer3W, Transformer2Wetap
    using ComponentModel: StaticGenerator, Generator, StaticGeneratorDC,Load, AsymmetricLoad, LoadDC
    using ComponentModel: Storage, MobileStorage, Storageetap, StorageAC
    using ComponentModel: ExternalGrid, Converter
    using ComponentModel: EquipmentCarbon, CarbonTimeSeries, CarbonScenario
    using ComponentModel: VirtualPowerPlant, FlexLoad
    using ComponentModel: ChargingStation, Charger, EVAggregator, V2GService
    using ComponentModel: Microgrid,PVArray, ACPVSystem


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

    export  read_load_data,read_price_data,read_irradiance_data
    export runtdpf, run_dynamic_dispatch, run_single_day,plot_voltage_time_series,plot_PD_time_series
    export record_voltage_violation, plot_losses_time_series, plot_flow_violations, 
        get_bus_voltage_results_acdc,analyze_voltage_results,process_result

    export topology_analysis, options, runhpf, runpf, rundcpf
    export load_julia_power_data
    export JuliaPowerCase2Jpc, JuliaPowerCase2Jpc_3ph,extract_islands_acdc,extract_islands
    export create_time_series_loads,create_time_series_prices,create_time_series_irradiance,extract_load_matrix_by_islands

    export AbstractComponent, Bus, Line, LineDC, Switch, BusDC, HighVoltageCircuitBreaker,Transformer2W, Transformer3W, Transformer2Wetap
    export StaticGenerator, Generator, StaticGeneratorDC,Load, AsymmetricLoad, LoadDC
    export Storage, MobileStorage, Storageetap, StorageAC
    export ExternalGrid, Converter
    export EquipmentCarbon, CarbonTimeSeries, CarbonScenario
    export VirtualPowerPlant, FlexLoad
    export ChargingStation, Charger, EVAggregator, V2GService
    export Microgrid,PVArray, ACPVSystem

    export ComponentModel
    export Utils
    export PowerFlow
    export TimeSeriesPowerFlow
end

