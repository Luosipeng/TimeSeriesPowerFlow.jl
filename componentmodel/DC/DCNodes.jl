# 节点结构
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
    
    # 构造函数
    function BusDC(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
        return new(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id)
    end
end