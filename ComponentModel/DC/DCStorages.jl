"""
    Storage <: AbstractComponent

Represents a stationary energy storage system in a power network.

# Fields
- `index::Int`: Unique identifier for the storage system
- `name::String`: Name of the storage system
- `bus::Int`: Bus identifier where the storage system is connected
- `power_capacity_mw::Float64`: Maximum power capacity in MW
- `energy_capacity_mwh::Float64`: Energy storage capacity in MWh
- `soc_init::Float64`: Initial state of charge (0-1)
- `min_soc::Float64`: Minimum allowed state of charge (0-1)
- `max_soc::Float64`: Maximum allowed state of charge (0-1)
- `efficiency::Float64`: Round-trip efficiency of the storage system (0-1)
- `in_service::Bool`: Operational status (true if in service)
- `type::String`: Type of storage technology (e.g., "battery", "pumped_hydro")
- `controllable::Bool`: Whether the storage system can be controlled externally
"""
mutable struct Storage <: AbstractComponent
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
        Storage(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
               min_soc, max_soc, efficiency, in_service, type, controllable)

    Create a new stationary energy storage system.

    # Parameters
    - `index`: Unique identifier for the storage system
    - `name`: Name of the storage system
    - `bus`: Bus identifier where the storage system is connected
    - `power_capacity_mw`: Maximum power capacity in MW
    - `energy_capacity_mwh`: Energy storage capacity in MWh
    - `soc_init`: Initial state of charge (0-1)
    - `min_soc`: Minimum allowed state of charge (0-1)
    - `max_soc`: Maximum allowed state of charge (0-1)
    - `efficiency`: Round-trip efficiency of the storage system (0-1)
    - `in_service`: Operational status (true if in service)
    - `type`: Type of storage technology (e.g., "battery", "pumped_hydro")
    - `controllable`: Whether the storage system can be controlled externally
    """
    function Storage(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                     min_soc, max_soc, efficiency, in_service, type, controllable)
        return new(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                  min_soc, max_soc, efficiency, in_service, type, controllable)
    end
end

"""
    MobileStorage <: AbstractComponent

Represents a mobile energy storage system that can be relocated within a power network.

# Fields
- `index::Int`: Unique identifier for the mobile storage system
- `name::String`: Name of the mobile storage system
- `type::String`: Type of mobile storage (container, trailer, truck, ship)
- `capacity_mwh::Float64`: Energy storage capacity in MWh
- `power_mw::Float64`: Maximum power capacity in MW
- `soc_percent::Float64`: Current state of charge as a percentage (0-100)
- `charge_efficiency_percent::Float64`: Charging efficiency as a percentage (0-100)
- `discharge_efficiency_percent::Float64`: Discharging efficiency as a percentage (0-100)
- `owner::String`: Owner or operator of the mobile storage system
- `availability::Float64`: Availability factor (0-1)
- `current_location::String`: Current physical location of the mobile storage
- `status::String`: Current operational status (available, in_transit, connected, maintenance)
"""
mutable struct MobileStorage <: AbstractComponent
    index::Int
    name::String
    type::String  # container, trailer, truck, ship
    capacity_mwh::Float64
    power_mw::Float64
    soc_percent::Float64
    charge_efficiency_percent::Float64
    discharge_efficiency_percent::Float64
    owner::String
    availability::Float64
    current_location::String
    status::String  # available, in_transit, connected, maintenance
    
    # Constructor
    """
        MobileStorage(index, name, type, capacity_mwh, power_mw, soc_percent,
                     charge_efficiency_percent, discharge_efficiency_percent,
                     owner, availability, current_location, status)

    Create a new mobile energy storage system.

    # Parameters
    - `index`: Unique identifier for the mobile storage system
    - `name`: Name of the mobile storage system
    - `type`: Type of mobile storage (container, trailer, truck, ship)
    - `capacity_mwh`: Energy storage capacity in MWh
    - `power_mw`: Maximum power capacity in MW
    - `soc_percent`: Current state of charge as a percentage (0-100)
    - `charge_efficiency_percent`: Charging efficiency as a percentage (0-100)
    - `discharge_efficiency_percent`: Discharging efficiency as a percentage (0-100)
    - `owner`: Owner or operator of the mobile storage system
    - `availability`: Availability factor (0-1)
    - `current_location`: Current physical location of the mobile storage
    - `status`: Current operational status (available, in_transit, connected, maintenance)
    """
    function MobileStorage(index, name, type, capacity_mwh, power_mw, soc_percent,
                          charge_efficiency_percent, discharge_efficiency_percent,
                          owner, availability, current_location, status)
        return new(index, name, type, capacity_mwh, power_mw, soc_percent,
                  charge_efficiency_percent, discharge_efficiency_percent,
                  owner, availability, current_location, status)
    end
end

"""
    Storageetap <: AbstractComponent

Represents a detailed electrochemical energy storage model with electrical and thermal parameters.

# Fields
- `index::Int`: Unique identifier for the storage system
- `name::String`: Name of the storage system
- `bus::Int`: Bus identifier where the storage system is connected
- `ra::Float64`: Internal resistance in ohms
- `cell::Float64`: Number of cells in series per string
- `str::Float64`: Number of strings in parallel
- `package::Float64`: Number of packages in parallel
- `voc::Float64`: Open circuit voltage in volts
- `in_service::Bool`: Operational status (true if in service)
- `type::String`: Type of storage technology
- `controllable::Bool`: Whether the storage system can be controlled externally
"""
mutable struct Storageetap <: AbstractComponent
    index::Int
    name::String
    bus::Int
    ra::Float64
    cell::Float64
    str::Float64
    package::Float64
    voc::Float64
    in_service::Bool
    type::String
    controllable::Bool
    
    # Constructor
    """
        Storageetap(index, name, bus, ra, cell, str, package, voc, in_service, type, controllable)

    Create a new detailed electrochemical energy storage model.

    # Parameters
    - `index`: Unique identifier for the storage system
    - `name`: Name of the storage system
    - `bus`: Bus identifier where the storage system is connected
    - `ra`: Internal resistance in ohms
    - `cell`: Number of cells in series per string
    - `str`: Number of strings in parallel
    - `package`: Number of packages in parallel
    - `voc`: Open circuit voltage in volts
    - `in_service`: Operational status (true if in service)
    - `type`: Type of storage technology
    - `controllable`: Whether the storage system can be controlled externally
    """
    function Storageetap(index, name, bus, ra, cell, str, package, voc, in_service, type, controllable)
        return new(index, name, bus, ra, cell, str, package, voc, in_service, type, controllable)
    end
end
