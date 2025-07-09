"""
    Transformer2W <: AbstractComponent

Represents a two-winding transformer in power systems.

# Fields
- `index::Int`: Unique identifier for the transformer
- `name::String`: Transformer name
- `std_type::String`: Standard type
- `hv_bus::Int`: High voltage bus number
- `lv_bus::Int`: Low voltage bus number
- `sn_mva::Float64`: Rated power (MVA)
- `vn_hv_kv::Float64`: Rated voltage at HV side (kV)
- `vn_lv_kv::Float64`: Rated voltage at LV side (kV)
- `vk_percent::Float64`: Short-circuit impedance in percent
- `vkr_percent::Float64`: Short-circuit resistance in percent
- `pfe_kw::Float64`: Iron losses (kW)
- `i0_percent::Float64`: No-load current in percent
- `shift_degree::Float64`: Phase shift angle in degrees

# Tap changer parameters
- `tap_side::String`: Tap changer side (hv/lv)
- `tap_neutral::Int`: Neutral tap position
- `tap_min::Int`: Minimum tap position
- `tap_max::Int`: Maximum tap position
- `tap_step_percent::Float64`: Tap step size in percent
- `tap_step_degree::Float64`: Tap step angle in degrees
- `tap_pos::Int`: Current tap position
- `tap_phase_shifter::Bool`: Whether it's a phase shifter

# Other parameters
- `parallel::Int`: Number of parallel transformers
- `max_loading_percent::Float64`: Maximum loading in percent
- `df::Float64`: Distribution factor
- `in_service::Bool`: Operating status flag
- `oltc::Bool`: Whether it has on-load tap changer
- `power_station_unit::Bool`: Whether it's a power station transformer

# Technical parameters
- `vector_group::String`: Vector group
- `hv_connection::String`: HV side connection type (Y/D)
- `lv_connection::String`: LV side connection type (Y/D)
- `thermal_capacity_mw::Float64`: Thermal capacity (MW)
- `cooling_type::String`: Cooling method
- `oil_volume_liters::Float64`: Oil volume (liters)
- `winding_material::String`: Winding material
"""
mutable struct Transformer2W <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    lv_bus::Int
    sn_mva::Float64
    vn_hv_kv::Float64
    vn_lv_kv::Float64
    vk_percent::Float64
    vkr_percent::Float64
    pfe_kw::Float64
    i0_percent::Float64
    shift_degree::Float64
    
    # Tap changer parameters
    tap_side::String
    tap_neutral::Int
    tap_min::Int
    tap_max::Int
    tap_step_percent::Float64
    tap_step_degree::Float64
    tap_pos::Int
    tap_phase_shifter::Bool
    
    parallel::Int
    max_loading_percent::Float64
    df::Float64
    in_service::Bool
    oltc::Bool
    power_station_unit::Bool
    
    # Technical parameters
    vector_group::String
    hv_connection::String
    lv_connection::String
    thermal_capacity_mw::Float64
    cooling_type::String
    oil_volume_liters::Float64
    winding_material::String
    
    # Constructor
    """
        Transformer2W(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                      vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree; kwargs...)

    Create a new two-winding transformer instance.

    # Parameters
    - `index`: Unique identifier for the transformer
    - `name`: Transformer name
    - `std_type`: Standard type
    - `hv_bus`: High voltage bus number
    - `lv_bus`: Low voltage bus number
    - `sn_mva`: Rated power (MVA)
    - `vn_hv_kv`: Rated voltage at HV side (kV)
    - `vn_lv_kv`: Rated voltage at LV side (kV)
    - `vk_percent`: Short-circuit impedance in percent
    - `vkr_percent`: Short-circuit resistance in percent
    - `pfe_kw`: Iron losses (kW)
    - `i0_percent`: No-load current in percent
    - `shift_degree`: Phase shift angle in degrees
    - `kwargs...`: Optional parameters

    # Optional parameters
    - `tap_side`: Tap changer side (hv/lv)
    - `tap_neutral`: Neutral tap position
    - `tap_min`: Minimum tap position
    - `tap_max`: Maximum tap position
    - `tap_step_percent`: Tap step size in percent
    - `tap_step_degree`: Tap step angle in degrees
    - `tap_pos`: Current tap position
    - `tap_phase_shifter`: Whether it's a phase shifter
    - `parallel`: Number of parallel transformers
    - `max_loading_percent`: Maximum loading in percent
    - `df`: Distribution factor
    - `in_service`: Operating status flag
    - `oltc`: Whether it has on-load tap changer
    - `power_station_unit`: Whether it's a power station transformer
    - `vector_group`: Vector group
    - `hv_connection`: HV side connection type (Y/D)
    - `lv_connection`: LV side connection type (Y/D)
    - `thermal_capacity_mw`: Thermal capacity (MW)
    - `cooling_type`: Cooling method
    - `oil_volume_liters`: Oil volume (liters)
    - `winding_material`: Winding material
    """
    function Transformer2W(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                          vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree; kwargs...)
        self = new(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                  vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree)
        
        # Set other parameters
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

"""
    Transformer2Wetap <: AbstractComponent

Represents a two-winding transformer with ETAP parameters.

# Fields
- `index::Int`: Unique identifier for the transformer
- `name::String`: Transformer name
- `std_type::String`: Standard type
- `hv_bus::Int`: High voltage bus number
- `lv_bus::Int`: Low voltage bus number
- `sn_mva::Float64`: Rated power (MVA)
- `vn_hv_kv::Float64`: Rated voltage at HV side (kV)
- `vn_lv_kv::Float64`: Rated voltage at LV side (kV)
- `z_percent::Float64`: Impedance in percent
- `x_r::Float64`: X/R ratio
- `z0_percent::Float64`: Zero sequence impedance in percent
- `x0_r0::Float64`: Zero sequence X/R ratio

# Tap changer parameters
- `tap_neutral::Int`: Neutral tap position
- `prim_tap::Float64`: Primary side tap
- `prim_tap_min::Int`: Primary side minimum tap position
- `prim_tap_max::Int`: Primary side maximum tap position
- `sec_tap::Float64`: Secondary side tap
- `sec_tap_min::Int`: Secondary side minimum tap position
- `sec_tap_max::Int`: Secondary side maximum tap position
- `vectororwinding::String`: Vector or winding type
- `phaseshifthl::Float64`: Phase shift angle between HV and LV sides
- `phaseshiftps`: Phase shifter angle

# Other parameters
- `parallel::Int`: Number of parallel transformers
- `max_loading_percent::Float64`: Maximum loading in percent
- `df::Float64`: Distribution factor
- `in_service::Bool`: Operating status flag
- `oltc::Bool`: Whether it has on-load tap changer
- `power_station_unit::Bool`: Whether it's a power station transformer
- `vector_group::String`: Vector group
"""
mutable struct Transformer2Wetap <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    lv_bus::Int
    sn_mva::Float64
    vn_hv_kv::Float64
    vn_lv_kv::Float64
    z_percent::Float64
    x_r::Float64
    z0_percent::Float64
    x0_r0::Float64
   
    
    # Tap changer parameters
    tap_neutral::Int
    prim_tap::Float64
    prim_tap_min::Int
    prim_tap_max::Int
    sec_tap::Float64
    sec_tap_min::Int
    sec_tap_max::Int
    vectororwinding::String
    phaseshifthl::Float64
    phaseshiftps::Float64
    
    parallel::Int
    max_loading_percent::Float64
    df::Float64
    in_service::Bool
    oltc::Bool
    power_station_unit::Bool
    
    # Technical parameters
    vector_group::String
    
    # Constructor
    """
        Transformer2Wetap(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                         z_percent, x_r, z0_percent, x0_r0, tap_neutral; kwargs...)

    Create a new two-winding transformer instance with ETAP parameters.

    # Parameters
    - `index`: Unique identifier for the transformer
    - `name`: Transformer name
    - `std_type`: Standard type
    - `hv_bus`: High voltage bus number
    - `lv_bus`: Low voltage bus number
    - `sn_mva`: Rated power (MVA)
    - `vn_hv_kv`: Rated voltage at HV side (kV)
    - `vn_lv_kv`: Rated voltage at LV side (kV)
    - `z_percent`: Impedance in percent
    - `x_r`: X/R ratio
    - `z0_percent`: Zero sequence impedance in percent
    - `x0_r0`: Zero sequence X/R ratio
    - `tap_neutral`: Neutral tap position
    - `kwargs...`: Optional parameters

    # Optional parameters
    - `prim_tap`: Primary side tap
    - `prim_tap_min`: Primary side minimum tap position
    - `prim_tap_max`: Primary side maximum tap position
    - `sec_tap`: Secondary side tap
    - `sec_tap_min`: Secondary side minimum tap position
    - `sec_tap_max`: Secondary side maximum tap position
    - `vectororwinding`: Vector or winding type
    - `phaseshifthl`: Phase shift angle between HV and LV sides
    - `phaseshiftps`: Phase shifter angle
    - `parallel`: Number of parallel transformers
    - `max_loading_percent`: Maximum loading in percent
    - `df`: Distribution factor
    - `in_service`: Operating status flag
    - `oltc`: Whether it has on-load tap changer
    - `power_station_unit`: Whether it's a power station transformer
    - `vector_group`: Vector group
    """
    function Transformer2Wetap(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
        z_percent, x_r, z0_percent, x0_r0, tap_neutral; kwargs...)
        self = new(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
            z_percent, x_r, z0_percent, x0_r0, tap_neutral)
        
        # Set other parameters
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

# Three-winding transformer structure
"""
    Transformer3W <: AbstractComponent

Represents a three-winding transformer in power systems.

# Fields
- `index::Int`: Unique identifier for the transformer
- `name::String`: Transformer name
- `std_type::String`: Standard type
- `hv_bus::Int`: High voltage bus number
- `mv_bus::Int`: Medium voltage bus number
- `lv_bus::Int`: Low voltage bus number
- `sn_hv_mva::Float64`: HV side rated power (MVA)
- `sn_mv_mva::Float64`: MV side rated power (MVA)
- `sn_lv_mva::Float64`: LV side rated power (MVA)
- `vn_hv_kv::Float64`: Rated voltage at HV side (kV)
- `vn_mv_kv::Float64`: Rated voltage at MV side (kV)
- `vn_lv_kv::Float64`: Rated voltage at LV side (kV)
- `vk_hv_percent::Float64`: HV side short-circuit impedance in percent
- `vk_mv_percent::Float64`: MV side short-circuit impedance in percent
- `vk_lv_percent::Float64`: LV side short-circuit impedance in percent
- `vkr_hv_percent::Float64`: HV side short-circuit resistance in percent
- `vkr_mv_percent::Float64`: MV side short-circuit resistance in percent
- `vkr_lv_percent::Float64`: LV side short-circuit resistance in percent
- `pfe_kw::Float64`: Iron losses (kW)
- `i0_percent::Float64`: No-load current in percent
- `shift_mv_degree::Float64`: MV side phase shift angle in degrees
- `shift_lv_degree::Float64`: LV side phase shift angle in degrees

# Tap changer parameters
- `tap_side::String`: Tap changer side (hv/mv/lv)
- `tap_neutral::Int`: Neutral tap position
- `tap_min::Int`: Minimum tap position
- `tap_max::Int`: Maximum tap position
- `tap_step_percent::Float64`: Tap step size in percent
- `tap_step_degree::Float64`: Tap step angle in degrees
- `tap_at_star_point::Bool`: Whether tap changer is at star point
- `tap_pos::Int`: Current tap position

# Other parameters
- `in_service::Bool`: Operating status flag

# Technical parameters
- `vector_group_hv_mv::String`: HV-MV vector group
- `vector_group_hv_lv::String`: HV-LV vector group
- `vector_group_mv_lv::String`: MV-LV vector group
- `hv_connection::String`: HV side connection type (Y/D)
- `mv_connection::String`: MV side connection type (Y/D)
- `lv_connection::String`: LV side connection type (Y/D)
- `thermal_capacity_mw::Float64`: Thermal capacity (MW)
- `cooling_type::String`: Cooling method
- `oil_volume_liters::Float64`: Oil volume (liters)
- `winding_material::String`: Winding material
"""
mutable struct Transformer3W <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    mv_bus::Int
    lv_bus::Int
    sn_hv_mva::Float64
    sn_mv_mva::Float64
    sn_lv_mva::Float64
    vn_hv_kv::Float64
    vn_mv_kv::Float64
    vn_lv_kv::Float64
    vk_hv_percent::Float64
    vk_mv_percent::Float64
    vk_lv_percent::Float64
    vkr_hv_percent::Float64
    vkr_mv_percent::Float64
    vkr_lv_percent::Float64
    pfe_kw::Float64
    i0_percent::Float64
    shift_mv_degree::Float64
    shift_lv_degree::Float64
    
    # Tap changer parameters
    tap_side::String
    tap_neutral::Int
    tap_min::Int
    tap_max::Int
    tap_step_percent::Float64
    tap_step_degree::Float64
    tap_at_star_point::Bool
    tap_pos::Int
    
    in_service::Bool
    
    # Technical parameters
    vector_group_hv_mv::String
    vector_group_hv_lv::String
    vector_group_mv_lv::String
    hv_connection::String
    mv_connection::String
    lv_connection::String
    thermal_capacity_mw::Float64
    cooling_type::String
    oil_volume_liters::Float64
    winding_material::String
    
    # Constructor
    """
        Transformer3W(index, name, std_type, hv_bus, mv_bus, lv_bus, 
                     sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                     vn_hv_kv, vn_mv_kv, vn_lv_kv,
                     vk_hv_percent, vk_mv_percent, vk_lv_percent,
                     vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                     pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; kwargs...)

    Create a new three-winding transformer instance.

    # Parameters
    - `index`: Unique identifier for the transformer
    - `name`: Transformer name
    - `std_type`: Standard type
    - `hv_bus`: High voltage bus number
    - `mv_bus`: Medium voltage bus number
    - `lv_bus`: Low voltage bus number
    - `sn_hv_mva`: HV side rated power (MVA)
    - `sn_mv_mva`: MV side rated power (MVA)
    - `sn_lv_mva`: LV side rated power (MVA)
    - `vn_hv_kv`: Rated voltage at HV side (kV)
    - `vn_mv_kv`: Rated voltage at MV side (kV)
    - `vn_lv_kv`: Rated voltage at LV side (kV)
    - `vk_hv_percent`: HV side short-circuit impedance in percent
    - `vk_mv_percent`: MV side short-circuit impedance in percent
    - `vk_lv_percent`: LV side short-circuit impedance in percent
    - `vkr_hv_percent`: HV side short-circuit resistance in percent
    - `vkr_mv_percent`: MV side short-circuit resistance in percent
    - `vkr_lv_percent`: LV side short-circuit resistance in percent
    - `pfe_kw`: Iron losses (kW)
    - `i0_percent`: No-load current in percent
    - `shift_mv_degree`: MV side phase shift angle in degrees
    - `shift_lv_degree`: LV side phase shift angle in degrees
    - `kwargs...`: Optional parameters

    # Optional parameters
    - `tap_side`: Tap changer side (hv/mv/lv)
    - `tap_neutral`: Neutral tap position
    - `tap_min`: Minimum tap position
    - `tap_max`: Maximum tap position
    - `tap_step_percent`: Tap step size in percent
    - `tap_step_degree`: Tap step angle in degrees
    - `tap_at_star_point`: Whether tap changer is at star point
    - `tap_pos`: Current tap position
    - `in_service`: Operating status flag
    - `vector_group_hv_mv`: HV-MV vector group
    - `vector_group_hv_lv`: HV-LV vector group
    - `vector_group_mv_lv`: MV-LV vector group
    - `hv_connection`: HV side connection type (Y/D)
    - `mv_connection`: MV side connection type (Y/D)
    - `lv_connection`: LV side connection type (Y/D)
    - `thermal_capacity_mw`: Thermal capacity (MW)
    - `cooling_type`: Cooling method
    - `oil_volume_liters`: Oil volume (liters)
    - `winding_material`: Winding material
    """
    function Transformer3W(index, name, std_type, hv_bus, mv_bus, lv_bus, 
                          sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                          vn_hv_kv, vn_mv_kv, vn_lv_kv,
                          vk_hv_percent, vk_mv_percent, vk_lv_percent,
                          vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                          pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; kwargs...)
        self = new(index, name, std_type, hv_bus, mv_bus, lv_bus, 
                  sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                  vn_hv_kv, vn_mv_kv, vn_lv_kv,
                  vk_hv_percent, vk_mv_percent, vk_lv_percent,
                  vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                  pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree)
        
        # Set other parameters
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end
