"""
    BusDC <: AbstractComponent

Represents a DC bus (node) in a power system network.

# Fields
- `index::Int`: Unique identifier for the bus
- `name::String`: Name of the bus
- `vn_kv::Float64`: Nominal voltage in kilovolts
- `max_vm_pu::Float64`: Maximum allowed voltage magnitude in per unit
- `min_vm_pu::Float64`: Minimum allowed voltage magnitude in per unit
- `in_service::Bool`: Operational status (true if in service)
- `bus_id::Int`: External bus identifier
- `area_id::Int`: Area identifier the bus belongs to
- `zone_id::Int`: Zone identifier the bus belongs to
"""
mutable struct BusDC <: AbstractComponent
    index::Int
    name::String
    vn_kv::Float64
    max_vm_pu::Float64
    min_vm_pu::Float64
    in_service::Bool
    bus_id::Int
    area_id::Int
    zone_id::Int
    
    # Constructor
    """
        BusDC(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)

    Create a new DC bus (node) in a power system network.

    # Parameters
    - `index`: Unique identifier for the bus
    - `name`: Name of the bus
    - `vn_kv`: Nominal voltage in kilovolts
    - `max_vm_pu`: Maximum allowed voltage magnitude in per unit
    - `min_vm_pu`: Minimum allowed voltage magnitude in per unit
    - `in_service`: Operational status (true if in service)
    - `bus_id`: External bus identifier
    - `area_id`: Area identifier the bus belongs to
    - `zone_id`: Zone identifier the bus belongs to
    """
    function BusDC(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
        return new(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
    end
end
