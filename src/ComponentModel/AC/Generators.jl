# Static Generator Structure
"""
    StaticGenerator <: AbstractComponent

Structure representing a static generator, used for modeling simplified generator units.

# Fields
- `index::Int`: Unique identifier for the generator
- `name::String`: Generator name
- `bus::Int`: Connected bus number
- `p_mw::Float64`: Active power output (MW)
- `q_mvar::Float64`: Reactive power output (MVar)
- `scaling::Float64`: Power scaling factor
- `max_p_mw::Float64`: Maximum active power limit (MW)
- `min_p_mw::Float64`: Minimum active power limit (MW)
- `max_q_mvar::Float64`: Maximum reactive power limit (MVar)
- `min_q_mvar::Float64`: Minimum reactive power limit (MVar)
- `k::Float64`: Coefficient k (typically used for control equations)
- `rx::Float64`: Resistance to reactance ratio
- `in_service::Bool`: Operating status flag
- `type::String`: Generator type (WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.)
- `controllable::Bool`: Controllability flag
"""
mutable struct StaticGenerator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    k::Float64
    rx::Float64
    in_service::Bool
    type::String  # WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.
    controllable::Bool
    
    # Constructor
    """
        StaticGenerator(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                        max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable)

    Create a new static generator instance.

    # Parameters
    - `index`: Unique identifier for the generator
    - `name`: Generator name
    - `bus`: Connected bus number
    - `p_mw`: Active power output (MW)
    - `q_mvar`: Reactive power output (MVar)
    - `scaling`: Power scaling factor
    - `max_p_mw`: Maximum active power limit (MW)
    - `min_p_mw`: Minimum active power limit (MW)
    - `max_q_mvar`: Maximum reactive power limit (MVar)
    - `min_q_mvar`: Minimum reactive power limit (MVar)
    - `k`: Coefficient k (typically used for control equations)
    - `rx`: Resistance to reactance ratio
    - `in_service`: Operating status flag
    - `type`: Generator type (WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.)
    - `controllable`: Controllability flag
    """
    function StaticGenerator(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                            max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable)
        return new(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                  max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable)
    end
end

# Conventional Generator Structure
"""
    Generator <: AbstractComponent

Structure representing a conventional generator unit, including detailed generator parameters and operational characteristics.

# Fields
- `index::Int`: Unique identifier for the generator
- `name::String`: Generator name
- `bus::Int`: Connected bus number
- `p_mw::Float64`: Active power output (MW)
- `vm_pu::Float64`: Voltage magnitude (per unit)
- `sn_mva::Float64`: Rated capacity (MVA)
- `scaling::Float64`: Power scaling factor
- `max_p_mw::Float64`: Maximum active power limit (MW)
- `min_p_mw::Float64`: Minimum active power limit (MW)
- `max_q_mvar::Float64`: Maximum reactive power limit (MVar)
- `min_q_mvar::Float64`: Minimum reactive power limit (MVar)
- `vn_kv::Float64`: Rated voltage (kV)
- `xdss_pu::Float64`: Subtransient direct-axis reactance (per unit)
- `rdss_pu::Float64`: Subtransient resistance (per unit)
- `cos_phi::Float64`: Power factor
- `controllable::Bool`: Controllability flag
- `in_service::Bool`: Operating status flag
- `type::String`: Generator type
- `generator_type::String`: Generator technology type
- `fuel_type::String`: Fuel type
- `startup_time_cold_h::Float64`: Cold startup time (hours)
- `startup_time_warm_h::Float64`: Warm startup time (hours)
- `startup_time_hot_h::Float64`: Hot startup time (hours)
- `min_up_time_h::Float64`: Minimum up time (hours)
- `min_down_time_h::Float64`: Minimum down time (hours)
- `ramp_up_rate_mw_per_min::Float64`: Ramp-up rate (MW/minute)
- `ramp_down_rate_mw_per_min::Float64`: Ramp-down rate (MW/minute)
- `efficiency_percent::Float64`: Efficiency percentage
- `heat_rate_mmbtu_per_mwh::Float64`: Heat rate (MMBtu/MWh)
- `co2_emission_rate_kg_per_mwh::Float64`: CO2 emission rate (kg/MWh)
"""
mutable struct Generator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    vm_pu::Float64
    sn_mva::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    vn_kv::Float64
    xdss_pu::Float64
    rdss_pu::Float64
    cos_phi::Float64
    controllable::Bool
    in_service::Bool
    type::String
    generator_type::String
    fuel_type::String
    
    # Startup and operation parameters
    startup_time_cold_h::Float64
    startup_time_warm_h::Float64
    startup_time_hot_h::Float64
    min_up_time_h::Float64
    min_down_time_h::Float64
    ramp_up_rate_mw_per_min::Float64
    ramp_down_rate_mw_per_min::Float64
    
    # Efficiency and emission parameters
    efficiency_percent::Float64
    heat_rate_mmbtu_per_mwh::Float64
    co2_emission_rate_kg_per_mwh::Float64
    
    # Constructor
    """
        Generator(index, name, bus, p_mw, vm_pu, sn_mva, scaling, 
                 max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                 vn_kv, xdss_pu, rdss_pu, cos_phi, controllable, in_service,
                 type, generator_type, fuel_type; kwargs...)

    Create a new conventional generator unit instance.

    # Parameters
    - `index`: Unique identifier for the generator
    - `name`: Generator name
    - `bus`: Connected bus number
    - `p_mw`: Active power output (MW)
    - `vm_pu`: Voltage magnitude (per unit)
    - `sn_mva`: Rated capacity (MVA)
    - `scaling`: Power scaling factor
    - `max_p_mw`: Maximum active power limit (MW)
    - `min_p_mw`: Minimum active power limit (MW)
    - `max_q_mvar`: Maximum reactive power limit (MVar)
    - `min_q_mvar`: Minimum reactive power limit (MVar)
    - `vn_kv`: Rated voltage (kV)
    - `xdss_pu`: Subtransient direct-axis reactance (per unit)
    - `rdss_pu`: Subtransient resistance (per unit)
    - `cos_phi`: Power factor
    - `controllable`: Controllability flag
    - `in_service`: Operating status flag
    - `type`: Generator type
    - `generator_type`: Generator technology type
    - `fuel_type`: Fuel type
    - `kwargs...`: Optional parameters including startup, operation, efficiency and emission parameters

    # Optional Parameters
    - `startup_time_cold_h`: Cold startup time (hours)
    - `startup_time_warm_h`: Warm startup time (hours)
    - `startup_time_hot_h`: Hot startup time (hours)
    - `min_up_time_h`: Minimum up time (hours)
    - `min_down_time_h`: Minimum down time (hours)
    - `ramp_up_rate_mw_per_min`: Ramp-up rate (MW/minute)
    - `ramp_down_rate_mw_per_min`: Ramp-down rate (MW/minute)
    - `efficiency_percent`: Efficiency percentage
    - `heat_rate_mmbtu_per_mwh`: Heat rate (MMBtu/MWh)
    - `co2_emission_rate_kg_per_mwh`: CO2 emission rate (kg/MWh)
    """
    function Generator(index, name, bus, p_mw, vm_pu, sn_mva, scaling, 
                      max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                      vn_kv, xdss_pu, rdss_pu, cos_phi, controllable, in_service,
                      type, generator_type, fuel_type; kwargs...)
        self = new(index, name, bus, p_mw, vm_pu, sn_mva, scaling, 
                  max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                  vn_kv, xdss_pu, rdss_pu, cos_phi, controllable, in_service,
                  type, generator_type, fuel_type)
        
        # Set other parameters
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

# External Grid Structure
"""
    ExternalGrid <: AbstractComponent

Structure representing an external grid, used for modeling equivalent representations of external power systems.

# Fields
- `index::Int`: Unique identifier for the external grid
- `name::String`: External grid name
- `bus::Int`: Connected bus number
- `vm_pu::Float64`: Voltage magnitude (per unit)
- `va_degree::Float64`: Voltage angle (degrees)
- `in_service::Bool`: Operating status flag
- `s_sc_max_mva::Float64`: Maximum short-circuit capacity (MVA)
- `s_sc_min_mva::Float64`: Minimum short-circuit capacity (MVA)
- `rx_max::Float64`: Maximum R/X ratio
- `rx_min::Float64`: Minimum R/X ratio
- `r0x0_max::Float64`: Maximum zero-sequence R/X ratio
- `x0x_max::Float64`: Maximum ratio of zero-sequence to positive-sequence reactance
- `controllable::Bool`: Controllability flag
"""
mutable struct ExternalGrid <: AbstractComponent
    index::Int
    name::String
    bus::Int
    vm_pu::Float64
    va_degree::Float64
    in_service::Bool
    s_sc_max_mva::Float64
    s_sc_min_mva::Float64
    rx_max::Float64
    rx_min::Float64
    r0x0_max::Float64
    x0x_max::Float64
    controllable::Bool
    
    # Constructor
    """
        ExternalGrid(index, name, bus, vm_pu, va_degree, in_service, 
                    s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable)

    Create a new external grid instance.

    # Parameters
    - `index`: Unique identifier for the external grid
    - `name`: External grid name
    - `bus`: Connected bus number
    - `vm_pu`: Voltage magnitude (per unit)
    - `va_degree`: Voltage angle (degrees)
    - `in_service`: Operating status flag
    - `s_sc_max_mva`: Maximum short-circuit capacity (MVA)
    - `s_sc_min_mva`: Minimum short-circuit capacity (MVA)
    - `rx_max`: Maximum R/X ratio
    - `rx_min`: Minimum R/X ratio
    - `r0x0_max`: Maximum zero-sequence R/X ratio
    - `x0x_max`: Maximum ratio of zero-sequence to positive-sequence reactance
    - `controllable`: Controllability flag
    """
    function ExternalGrid(index, name, bus, vm_pu, va_degree, in_service, 
                         s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable)
        return new(index, name, bus, vm_pu, va_degree, in_service, 
                  s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable)
    end
end

# AC PV System
"""
    ACPVSystem <: AbstractComponent

Structure representing an AC photovoltaic system, including PV module and inverter parameters.

# Fields
- `index::Int`: Unique identifier for the PV system
- `name::String`: PV system name
- `bus::Int`: Connected bus number
- `p_mw::Float64`: Active power output (MW)
- `q_mvar::Float64`: Reactive power output (MVar)
- `vm_ac_pu::Float64`: AC-side voltage magnitude (per unit)
- `vm_dc_pu::Float64`: DC-side voltage magnitude (per unit)
- `loss_percent::Float64`: Loss percentage
- `loss_mw::Float64`: Power loss (MW)
- `max_p_mw::Float64`: Maximum active power limit (MW)
- `min_p_mw::Float64`: Minimum active power limit (MW)
- `max_q_mvar::Float64`: Maximum reactive power limit (MVar)
- `min_q_mvar::Float64`: Minimum reactive power limit (MVar)
- `numpanelseries::Int`: Number of panels in series
- `numpanelparallel::Int`: Number of panels in parallel
- `vmpp::Float64`: Maximum power point voltage
- `impp::Float64`: Maximum power point current
- `voc::Float64`: Open-circuit voltage
- `isc::Float64`: Short-circuit current
- `pvanumcells::Int`: Number of PV cells
- `irradiance::Float64`: Irradiance (W/m²)
- `temperature::Float64`: Temperature (°C)
- `α_isc::Float64`: Short-circuit current temperature coefficient
- `β_voc::Float64`: Open-circuit voltage temperature coefficient
- `control_mode::String`: Control mode
- `controllable::Bool`: Controllability flag
- `in_service::Bool`: Operating status flag
"""
mutable struct ACPVSystem <: AbstractComponent
    index::Int
    name::String
    bus::Int
    # Inverter parameters
    p_mw::Float64
    q_mvar::Float64
    vm_ac_pu::Float64
    vm_dc_pu::Float64
    loss_percent::Float64
    loss_mw::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    # PV module parameters
    numpanelseries::Int
    numpanelparallel::Int
    vmpp::Float64
    impp::Float64
    voc::Float64
    isc::Float64
    pvanumcells::Int
    irradiance::Float64
    temperature::Float64
    α_isc::Float64
    β_voc::Float64
    control_mode::String
    controllable::Bool
    in_service::Bool
    
    # Constructor
    """
        ACPVSystem(index, name, bus, p_mw, q_mvar, vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                   max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                   numpanelseries, numpanelparallel, vmpp, impp, voc, isc,
                   pvanumcells, irradiance, temperature, α_isc, β_voc,
                   control_mode, controllable, in_service)

    Create a new AC photovoltaic system instance.

    # Parameters
    - `index`: Unique identifier for the PV system
    - `name`: PV system name
    - `bus`: Connected bus number
    - `p_mw`: Active power output (MW)
    - `q_mvar`: Reactive power output (MVar)
    - `vm_ac_pu`: AC-side voltage magnitude (per unit)
    - `vm_dc_pu`: DC-side voltage magnitude (per unit)
    - `loss_percent`: Loss percentage
    - `loss_mw`: Power loss (MW)
    - `max_p_mw`: Maximum active power limit (MW)
    - `min_p_mw`: Minimum active power limit (MW)
    - `max_q_mvar`: Maximum reactive power limit (MVar)
    - `min_q_mvar`: Minimum reactive power limit (MVar)
    - `numpanelseries`: Number of panels in series
    - `numpanelparallel`: Number of panels in parallel
    - `vmpp`: Maximum power point voltage
    - `impp`: Maximum power point current
    - `voc`: Open-circuit voltage
    - `isc`: Short-circuit current
    - `pvanumcells`: Number of PV cells
    - `irradiance`: Irradiance (W/m²)
    - `temperature`: Temperature (°C)
    - `α_isc`: Short-circuit current temperature coefficient
    - `β_voc`: Open-circuit voltage temperature coefficient
    - `control_mode`: Control mode
    - `controllable`: Controllability flag
    - `in_service`: Operating status flag
    """
    function ACPVSystem(index, name, bus, p_mw, q_mvar, vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                        max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                        numpanelseries, numpanelparallel, vmpp, impp, voc, isc,
                        pvanumcells, irradiance, temperature, α_isc, β_voc,
                        control_mode, controllable, in_service)
        return new(index, name, bus, p_mw, q_mvar, vm_ac_pu, vm_dc_pu,
                  loss_percent, loss_mw, max_p_mw, min_p_mw,
                  max_q_mvar, min_q_mvar,
                  numpanelseries, numpanelparallel,
                  vmpp, impp, voc, isc,
                  pvanumcells,
                  irradiance, temperature,
                  α_isc, β_voc,
                  control_mode,
                  controllable,
                  in_service)
    end

end

# Storage Structure
"""
    StorageAC <: AbstractComponent

Structure representing an AC storage system.

# Fields
- `index::Int`: Unique identifier for the storage system
- `name::String`: Storage system name
- `bus::Int`: Connected bus number
- `power_capacity_mw::Float64`: Power capacity (MW)
- `energy_capacity_mwh::Float64`: Energy capacity (MWh)
- `soc_init::Float64`: Initial state of charge
- `min_soc::Float64`: Minimum state of charge
- `max_soc::Float64`: Maximum state of charge
- `efficiency::Float64`: Charging/discharging efficiency
- `in_service::Bool`: Operating status flag
- `type::String`: Storage type
- `controllable::Bool`: Controllability flag
"""
mutable struct StorageAC <: AbstractComponent
    index::Int
    name::String
    bus::Int
    power_capacity_mw::Float64
    energy_capacity_mwh::Float64
    soc_init::Float64
    min_soc::Float64
    max_soc::Float64
    efficiency::Float64
    in_service::Bool
    type::String
    controllable::Bool
    
    # Constructor
    """
        StorageAC(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                 min_soc, max_soc, efficiency, in_service, type, controllable)

    Create a new AC storage system instance.

    # Parameters
    - `index`: Unique identifier for the storage system
    - `name`: Storage system name
    - `bus`: Connected bus number
    - `power_capacity_mw`: Power capacity (MW)
    - `energy_capacity_mwh`: Energy capacity (MWh)
    - `soc_init`: Initial state of charge
    - `min_soc`: Minimum state of charge
    - `max_soc`: Maximum state of charge
    - `efficiency`: Charging/discharging efficiency
    - `in_service`: Operating status flag
    - `type`: Storage type
    - `controllable`: Controllability flag
    """
    function StorageAC(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                     min_soc, max_soc, efficiency, in_service, type, controllable)
        return new(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                  min_soc, max_soc, efficiency, in_service, type, controllable)
    end
end
