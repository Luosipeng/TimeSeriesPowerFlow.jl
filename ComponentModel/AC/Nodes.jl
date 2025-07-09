"""
    Bus <: AbstractComponent

Structure representing a bus (node) in power systems.

# Fields
- `index::Int`: Unique identifier for the bus
- `name::String`: Bus name
- `vn_kv::Float64`: Nominal voltage level in kV
- `max_vm_pu::Float64`: Maximum voltage magnitude in per unit
- `min_vm_pu::Float64`: Minimum voltage magnitude in per unit
- `in_service::Bool`: Operating status flag
- `bus_id::Int`: Bus identifier
- `area_id::Int`: Area identifier the bus belongs to
- `zone_id::Int`: Zone identifier the bus belongs to
"""
mutable struct Bus <: AbstractComponent
    index::Int
    name::String
    vn_kv::Float64
    max_vm_pu::Float64
    min_vm_pu::Float64
    in_service::Bool
    bus_id::Int
    area_id::Int
    zone_id::Int
    
    """
        Bus(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)

    Create a new bus instance.

    # Parameters
    - `index`: Unique identifier for the bus
    - `name`: Bus name
    - `vn_kv`: Nominal voltage level in kV
    - `max_vm_pu`: Maximum voltage magnitude in per unit
    - `min_vm_pu`: Minimum voltage magnitude in per unit
    - `in_service`: Operating status flag
    - `bus_id`: Bus identifier
    - `area_id`: Area identifier the bus belongs to
    - `zone_id`: Zone identifier the bus belongs to
    """
    function Bus(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
        return new(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
    end
end
