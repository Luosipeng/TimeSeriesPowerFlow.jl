"""
    ChargingStation Structure
    
    Represents an electric vehicle charging station with multiple chargers and connection
    to the power system at a specific bus.
"""
mutable struct ChargingStation <: AbstractComponent
    index::Int
    name::String
    bus::Int
    location::String
    operator::String
    num_chargers::Int
    max_power_kw::Float64
    in_service::Bool
    
    # Constructor
    function ChargingStation(index, name, bus, location, operator, num_chargers, max_power_kw, in_service)
        return new(index, name, bus, location, operator, num_chargers, max_power_kw, in_service)
    end
end

"""
    Charger Structure
    
    Represents an individual charging unit within a charging station, with specific power
    capabilities and connector types.
"""
mutable struct Charger <: AbstractComponent
    index::Int
    name::String
    station_id::Int
    type::String  # AC_Level1, AC_Level2, DC_Fast, Ultra_Fast
    max_power_kw::Float64
    min_power_kw::Float64
    connector_type::String
    v2g_capable::Bool
    in_service::Bool
    availability::Float64
    
    # Constructor
    function Charger(index, name, station_id, type, max_power_kw, min_power_kw,
                    connector_type, v2g_capable, in_service, availability)
        return new(index, name, station_id, type, max_power_kw, min_power_kw,
                  connector_type, v2g_capable, in_service, availability)
    end
end

"""
    EVAggregator Structure
    
    Represents an entity that aggregates multiple electric vehicles to provide
    grid services and manage charging/discharging operations collectively.
"""
mutable struct EVAggregator <: AbstractComponent
    index::Int
    name::String
    num_evs::Int
    total_capacity_mwh::Float64
    max_power_mw::Float64
    control_strategy::String
    service_area::String
    operator::String
    
    # Constructor
    function EVAggregator(index, name, num_evs, total_capacity_mwh, max_power_mw,
                         control_strategy, service_area, operator)
        return new(index, name, num_evs, total_capacity_mwh, max_power_mw,
                  control_strategy, service_area, operator)
    end
end

"""
    V2GService Structure
    
    Represents a vehicle-to-grid service provided by an aggregator, including
    service type, time period, capacity, and pricing information.
"""
mutable struct V2GService <: AbstractComponent
    index::Int
    aggregator_id::Int
    service_type::String  # frequency_regulation, peak_shaving, voltage_support, emergency_backup
    start_time::DateTime
    end_time::DateTime
    capacity_mw::Float64
    energy_mwh::Float64
    price_per_mwh::Float64
    grid_area::String
    reliability_percent::Float64
    
    # Constructor
    function V2GService(index, aggregator_id, service_type, start_time, end_time,
                       capacity_mw, energy_mwh, price_per_mwh, grid_area, reliability_percent)
        return new(index, aggregator_id, service_type, start_time, end_time,
                  capacity_mw, energy_mwh, price_per_mwh, grid_area, reliability_percent)
    end
end
