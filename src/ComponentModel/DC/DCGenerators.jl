"""
    StaticGeneratorDC <: AbstractComponent

Represents a static generator in a DC power flow model.

# Fields
- `index::Int`: Unique identifier for the generator
- `name::String`: Name of the generator
- `bus::Int`: Bus identifier where the generator is connected
- `p_mw::Float64`: Active power output in MW
- `scaling::Float64`: Scaling factor for the generator output
- `max_p_mw::Float64`: Maximum active power output in MW
- `min_p_mw::Float64`: Minimum active power output in MW
- `k::Float64`: Participation factor for power balancing
- `rx::Float64`: R/X ratio of the generator
- `in_service::Bool`: Operational status (true if in service)
- `type::String`: Generator type (WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.)
- `controllable::Bool`: Whether the generator output can be controlled
"""
mutable struct StaticGeneratorDC <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    k::Float64
    rx::Float64
    in_service::Bool
    type::String  # WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.
    controllable::Bool
    
    # Constructor
    """
        StaticGeneratorDC(index, name, bus, p_mw, scaling, max_p_mw, min_p_mw, k, rx, in_service, type, controllable)

    Create a new static generator for DC power flow analysis.

    # Parameters
    - `index`: Unique identifier for the generator
    - `name`: Name of the generator
    - `bus`: Bus identifier where the generator is connected
    - `p_mw`: Active power output in MW
    - `scaling`: Scaling factor for the generator output
    - `max_p_mw`: Maximum active power output in MW
    - `min_p_mw`: Minimum active power output in MW
    - `k`: Participation factor for power balancing
    - `rx`: R/X ratio of the generator
    - `in_service`: Operational status (true if in service)
    - `type`: Generator type (WP-wind power, PV-photovoltaic, CHP-combined heat and power, etc.)
    - `controllable`: Whether the generator output can be controlled
    """
    function StaticGeneratorDC(index, name, bus, p_mw, scaling, max_p_mw, min_p_mw, k, rx, in_service, type, controllable)
        return new(index, name, bus, p_mw, scaling, max_p_mw, min_p_mw, k, rx, in_service, type, controllable)
    end
end
