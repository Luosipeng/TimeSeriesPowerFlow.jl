"""
    Microgrid Modeling
    
    This module provides structures for modeling microgrids and their components
    in power system simulations.
"""

"""
    Microgrid Structure
    
    Represents a microgrid with generation, storage, and load resources that can
    operate connected to or isolated from the main grid.
"""
mutable struct Microgrid <: AbstractComponent
    index::Int
    name::String
    description::String
    control_area::String
    capacity_mw::Float64
    energy_mwh::Float64
    response_time_s::Float64
    ramp_rate_mw_per_min::Float64
    availability_percent::Float64
    operator::String
    in_service::Bool
    
    # Resource information
    resource_type::String
    resource_id::Int
    capacity_share_percent::Float64
    control_priority::Int
    resource_response_time_s::Float64
    max_duration_h::Float64
    
    # Load information
    timestamp::DateTime
    p_mw::Float64
    q_mvar::Float64
    flexibility_up_mw::Float64
    flexibility_down_mw::Float64
    flexibility_duration_h::Float64
    
    # Constructor
    function VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                              response_time_s, ramp_rate_mw_per_min, availability_percent,
                              operator, in_service; kwargs...)
        self = new(index, name, description, control_area, capacity_mw, energy_mwh,
                  response_time_s, ramp_rate_mw_per_min, availability_percent,
                  operator, in_service)
        
        # Set other parameters
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end
