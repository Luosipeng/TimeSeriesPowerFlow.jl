"""
    LineDC <: AbstractComponent

Represents a DC transmission line in a power system.

# Fields
- `index::Int`: Unique identifier for the line
- `name::String`: Name of the line
- `from_bus::Int`: Bus identifier where the line starts
- `to_bus::Int`: Bus identifier where the line ends
- `length_km::Float64`: Length of the line in kilometers
- `r_ohm_per_km::Float64`: Resistance per kilometer in ohms
- `x_ohm_per_km::Float64`: Reactance per kilometer in ohms
- `g_us_per_km::Float64`: Conductance per kilometer in microsiemens
- `max_i_ka::Float64`: Maximum current capacity in kiloamperes
- `type::String`: Type of the line
- `max_loading_percent::Float64`: Maximum loading percentage allowed
- `parallel::Int`: Number of parallel lines
- `df::Float64`: Derating factor
- `in_service::Bool`: Operational status (true if in service)
- `mtbf_hours::Float64`: Mean time between failures in hours
- `mttr_hours::Float64`: Mean time to repair in hours
- `sw_hours::Float64`: Switching time in hours
- `rp_hours::Float64`: Replacement time in hours
"""
mutable struct LineDC <: AbstractComponent
    index::Int
    name::String
    from_bus::Int
    to_bus::Int
    length_km::Float64
    r_ohm_per_km::Float64
    x_ohm_per_km::Float64
    g_us_per_km::Float64
    max_i_ka::Float64
    type::String
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
        LineDC(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
              g_us_per_km, max_i_ka, type, max_loading_percent, 
              parallel, df, in_service)

    Create a new DC transmission line.

    # Parameters
    - `index`: Unique identifier for the line
    - `name`: Name of the line
    - `from_bus`: Bus identifier where the line starts
    - `to_bus`: Bus identifier where the line ends
    - `length_km`: Length of the line in kilometers
    - `r_ohm_per_km`: Resistance per kilometer in ohms
    - `x_ohm_per_km`: Reactance per kilometer in ohms
    - `g_us_per_km`: Conductance per kilometer in microsiemens
    - `max_i_ka`: Maximum current capacity in kiloamperes
    - `type`: Type of the line
    - `max_loading_percent`: Maximum loading percentage allowed
    - `parallel`: Number of parallel lines
    - `df`: Derating factor
    - `in_service`: Operational status (true if in service)

    Note: Reliability parameters (mtbf_hours, mttr_hours, sw_hours, rp_hours) need to be set separately after creation.
    """
    function LineDC(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
                   g_us_per_km, max_i_ka, type, max_loading_percent, 
                   parallel, df, in_service)
        return new(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                  g_us_per_km, max_i_ka, type, max_loading_percent,
                  parallel, df, in_service)
    end
end
