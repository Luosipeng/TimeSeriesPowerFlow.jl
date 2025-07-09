# AC Line Structure
"""
    Line <: AbstractComponent

Structure representing AC transmission or distribution lines in power systems.

# Fields
- `index::Int`: Unique identifier for the line
- `name::String`: Line name
- `from_bus::Int`: Starting bus number
- `to_bus::Int`: Ending bus number
- `length_km::Float64`: Line length (km)
- `r_ohm_per_km::Float64`: Positive sequence resistance per kilometer (Ω/km)
- `x_ohm_per_km::Float64`: Positive sequence reactance per kilometer (Ω/km)
- `c_nf_per_km::Float64`: Positive sequence capacitance per kilometer (nF/km)
- `r0_ohm_per_km::Float64`: Zero sequence resistance per kilometer (Ω/km)
- `x0_ohm_per_km::Float64`: Zero sequence reactance per kilometer (Ω/km)
- `c0_nf_per_km::Float64`: Zero sequence capacitance per kilometer (nF/km)
- `g_us_per_km::Float64`: Conductance per kilometer (μS/km)
- `max_i_ka::Float64`: Maximum allowable current (kA)
- `type::String`: Line type (cs-cable, ol-overhead line)
- `max_loading_percent::Float64`: Maximum loading percentage
- `parallel::Int`: Number of parallel lines
- `df::Float64`: Distribution factor
- `in_service::Bool`: Operating status flag

# Reliability Parameters
- `mtbf_hours::Float64`: Mean time between failures (hours)
- `mttr_hours::Float64`: Mean time to repair (hours)
- `sw_hours::Float64`: Switching operation time (hours)
- `rp_hours::Float64`: Repair preparation time (hours)
"""
mutable struct Line <: AbstractComponent
    index::Int
    name::String
    from_bus::Int
    to_bus::Int
    length_km::Float64
    r_ohm_per_km::Float64
    x_ohm_per_km::Float64
    c_nf_per_km::Float64
    r0_ohm_per_km::Float64
    x0_ohm_per_km::Float64
    c0_nf_per_km::Float64
    g_us_per_km::Float64
    max_i_ka::Float64
    type::String  # cs-cable, ol-overhead line
    max_loading_percent::Float64
    parallel::Int
    df::Float64
    in_service::Bool
    
    # Reliability parameters
    mtbf_hours::Float64
    mttr_hours::Float64
    sw_hours::Float64
    rp_hours::Float64
    
    # Constructor
    """
        Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
             c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km, 
             c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent, 
             parallel, df, in_service; reliability_params...)

    Create a new AC line instance.

    # Parameters
    - `index`: Unique identifier for the line
    - `name`: Line name
    - `from_bus`: Starting bus number
    - `to_bus`: Ending bus number
    - `length_km`: Line length (km)
    - `r_ohm_per_km`: Positive sequence resistance per kilometer (Ω/km)
    - `x_ohm_per_km`: Positive sequence reactance per kilometer (Ω/km)
    - `c_nf_per_km`: Positive sequence capacitance per kilometer (nF/km)
    - `r0_ohm_per_km`: Zero sequence resistance per kilometer (Ω/km)
    - `x0_ohm_per_km`: Zero sequence reactance per kilometer (Ω/km)
    - `c0_nf_per_km`: Zero sequence capacitance per kilometer (nF/km)
    - `g_us_per_km`: Conductance per kilometer (μS/km)
    - `max_i_ka`: Maximum allowable current (kA)
    - `type`: Line type (cs-cable, ol-overhead line)
    - `max_loading_percent`: Maximum loading percentage
    - `parallel`: Number of parallel lines
    - `df`: Distribution factor
    - `in_service`: Operating status flag
    - `reliability_params...`: Optional reliability parameters

    # Optional Parameters
    - `mtbf_hours`: Mean time between failures (hours)
    - `mttr_hours`: Mean time to repair (hours)
    - `sw_hours`: Switching operation time (hours)
    - `rp_hours`: Repair preparation time (hours)
    """
    function Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
                 c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km, 
                 c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent, 
                 parallel, df, in_service; reliability_params...)
        self = new(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                  c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km, 
                  c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent,
                  parallel, df, in_service)
        
        # Set reliability parameters
        for (key, value) in reliability_params
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end





# Switch Structure
"""
    Switch <: AbstractComponent

Structure representing switching devices in power systems.

# Fields
- `index::Int`: Unique identifier for the switch
- `name::String`: Switch name
- `bus_from::Int`: Starting bus number
- `bus_to::Int`: Ending bus number
- `element_type::String`: Connected element type (l-line, t-transformer, b-bus)
- `element_id::Int`: Connected element identifier
- `closed::Bool`: Switch status (closed/open)
- `type::String`: Switch type (CB-circuit breaker, LS-load switch, DS-disconnector switch)
- `z_ohm::Float64`: Switch impedance (Ω)
- `in_service::Bool`: Operating status flag
"""
mutable struct Switch <: AbstractComponent
    index::Int
    name::String
    bus_from::Int
    bus_to::Int
    element_type::String  # l-line, t-transformer, b-bus
    element_id::Int
    closed::Bool
    type::String  # CB-circuit breaker, LS-load switch, DS-disconnector switch
    z_ohm::Float64
    in_service::Bool
    
    # Constructor
    """
        Switch(index, name, bus_from, bus_to, element_type, element_id, 
               closed, type, z_ohm, in_service)

    Create a new switch instance.

    # Parameters
    - `index`: Unique identifier for the switch
    - `name`: Switch name
    - `bus_from`: Starting bus number
    - `bus_to`: Ending bus number
    - `element_type`: Connected element type (l-line, t-transformer, b-bus)
    - `element_id`: Connected element identifier
    - `closed`: Switch status (closed/open)
    - `type`: Switch type (CB-circuit breaker, LS-load switch, DS-disconnector switch)
    - `z_ohm`: Switch impedance (Ω)
    - `in_service`: Operating status flag
    """
    function Switch(index, name, bus_from, bus_to, element_type, element_id, 
                   closed, type, z_ohm, in_service)
        return new(index, name, bus_from, bus_to, element_type, element_id, 
                  closed, type, z_ohm, in_service)
    end
end

# High Voltage Circuit Breaker Structure
"""
    HighVoltageCircuitBreaker <: AbstractComponent

Structure representing high voltage circuit breakers in power systems.

# Fields
- `index::Int`: Unique identifier for the circuit breaker
- `name::String`: Circuit breaker name
- `bus_from::Int`: Starting bus number
- `bus_to::Int`: Ending bus number
- `element_type::String`: Connected element type (l-line, t-transformer, b-bus)
- `element_id::Int`: Connected element identifier
- `closed::Bool`: Circuit breaker status (closed/open)
- `type::String`: Circuit breaker type (CB-circuit breaker, LS-load switch, DS-disconnector switch)
- `z_ohm::Float64`: Circuit breaker impedance (Ω)
- `in_service::Bool`: Operating status flag
"""
mutable struct HighVoltageCircuitBreaker <: AbstractComponent
    index::Int
    name::String
    bus_from::Int
    bus_to::Int
    element_type::String  # l-line, t-transformer, b-bus
    element_id::Int
    closed::Bool
    type::String  # CB-circuit breaker, LS-load switch, DS-disconnector switch
    z_ohm::Float64
    in_service::Bool
    
    # Constructor
    """
        HighVoltageCircuitBreaker(index, name, bus_from, bus_to, element_type, element_id, 
                                 closed, type, z_ohm, in_service)

    Create a new high voltage circuit breaker instance.

    # Parameters
    - `index`: Unique identifier for the circuit breaker
    - `name`: Circuit breaker name
    - `bus_from`: Starting bus number
    - `bus_to`: Ending bus number
    - `element_type`: Connected element type (l-line, t-transformer, b-bus)
    - `element_id`: Connected element identifier
    - `closed`: Circuit breaker status (closed/open)
    - `type`: Circuit breaker type (CB-circuit breaker, LS-load switch, DS-disconnector switch)
    - `z_ohm`: Circuit breaker impedance (Ω)
    - `in_service`: Operating status flag
    """
    function HighVoltageCircuitBreaker(index, name, bus_from, bus_to, element_type, element_id, 
                   closed, type, z_ohm, in_service)
        return new(index, name, bus_from, bus_to, element_type, element_id, 
                  closed, type, z_ohm, in_service)
    end
end
