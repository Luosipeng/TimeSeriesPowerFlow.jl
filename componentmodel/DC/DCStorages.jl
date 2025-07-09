# 储能结构
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
    
    # 构造函数
    function Storage(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                     min_soc, max_soc, efficiency, in_service, type, controllable)
        return new(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                  min_soc, max_soc, efficiency, in_service, type, controllable)
    end
end

# 移动储能结构
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
    
    # 构造函数
    function MobileStorage(index, name, type, capacity_mwh, power_mw, soc_percent,
                          charge_efficiency_percent, discharge_efficiency_percent,
                          owner, availability, current_location, status)
        return new(index, name, type, capacity_mwh, power_mw, soc_percent,
                  charge_efficiency_percent, discharge_efficiency_percent,
                  owner, availability, current_location, status)
    end
end

# 储能结构
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
    
    # 构造函数
    function Storageetap(index, name, bus, ra, cell, str, package, voc, in_service, type, controllable)
        return new(index, name, bus, ra, cell, str, package, voc, in_service, type, controllable)
    end
end