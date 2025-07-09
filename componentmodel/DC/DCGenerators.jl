# 静态发电机结构
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
    type::String  # WP-风电, PV-光伏, CHP-热电联产等
    controllable::Bool
    
    # 构造函数
    function StaticGeneratorDC(index, name, bus, p_mw, scaling, max_p_mw, min_p_mw, k, rx, in_service, type, controllable)
        return new(index, name, bus, p_mw, scaling, max_p_mw, min_p_mw, k, rx, in_service, type, controllable)
    end
end
