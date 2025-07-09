# 交流线路结构
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
    type::String  # cs-电缆, ol-架空线
    max_loading_percent::Float64
    parallel::Int
    df::Float64
    in_service::Bool
    
    # 可靠性参数
    mtbf_hours::Float64
    mttr_hours::Float64
    sw_hours::Float64
    rp_hours::Float64
    
    # 构造函数
    function Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
                 c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km, 
                 c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent, 
                 parallel, df, in_service; reliability_params...)
        self = new(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                  c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km, 
                  c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent,
                  parallel, df, in_service)
        
        # 设置可靠性参数
        for (key, value) in reliability_params
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

# # 交流线路结构
# mutable struct ThreePhaseBranch <: AbstractComponent
#     index::Int
#     name::String
#     from_bus::Int
#     to_bus::Int
#     length_km::Float64
#     r_ohm_per_km::Float64
#     x_ohm_per_km::Float64
#     c_nf_per_km::Float64
#     g_us_per_km::Float64
#     max_i_ka::Float64
#     type::String  # cs-电缆, ol-架空线
#     max_loading_percent::Float64
#     parallel::Int
#     df::Float64
#     in_service::Bool
    
#     # 可靠性参数
#     mtbf_hours::Float64
#     mttr_hours::Float64
#     sw_hours::Float64
#     rp_hours::Float64
    
#     # 构造函数
#     function Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
#                  c_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent, 
#                  parallel, df, in_service; reliability_params...)
#         self = new(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
#                   c_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent,
#                   parallel, df, in_service)
        
#         # 设置可靠性参数
#         for (key, value) in reliability_params
#             if hasfield(typeof(self), key)
#                 setfield!(self, key, value)
#             end
#         end
        
#         return self
#     end
# end

# 两卷变压器结构
mutable struct Transformer2W <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    lv_bus::Int
    sn_mva::Float64
    vn_hv_kv::Float64
    vn_lv_kv::Float64
    vk_percent::Float64
    vkr_percent::Float64
    pfe_kw::Float64
    i0_percent::Float64
    shift_degree::Float64
    
    # 分接头参数
    tap_side::String
    tap_neutral::Int
    tap_min::Int
    tap_max::Int
    tap_step_percent::Float64
    tap_step_degree::Float64
    tap_pos::Int
    tap_phase_shifter::Bool
    
    parallel::Int
    max_loading_percent::Float64
    df::Float64
    in_service::Bool
    oltc::Bool
    power_station_unit::Bool
    
    # 技术参数
    vector_group::String
    hv_connection::String
    lv_connection::String
    thermal_capacity_mw::Float64
    cooling_type::String
    oil_volume_liters::Float64
    winding_material::String
    
    # 构造函数
    function Transformer2W(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                          vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree; kwargs...)
        self = new(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
                  vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

mutable struct Transformer2Wetap <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    lv_bus::Int
    sn_mva::Float64
    vn_hv_kv::Float64
    vn_lv_kv::Float64
    z_percent::Float64
    x_r::Float64
    z0_percent::Float64
    x0_r0::Float64
   
    
    # 分接头参数
    tap_neutral::Int
    prim_tap::Float64
    prim_tap_min::Int
    prim_tap_max::Int
    sec_tap::Float64
    sec_tap_min::Int
    sec_tap_max::Int
    vectororwinding::String
    phaseshifthl::Float64
    phaseshiftps::Float64
    
    parallel::Int
    max_loading_percent::Float64
    df::Float64
    in_service::Bool
    oltc::Bool
    power_station_unit::Bool
    
    # 技术参数
    vector_group::String
    
    # 构造函数
    function Transformer2Wetap(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
        z_percent, x_r, z0_percent, x0_r0, tap_neutral; kwargs...)
        self = new(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
            z_percent, x_r, z0_percent, x0_r0, tap_neutral)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

# 三卷变压器结构
mutable struct Transformer3W <: AbstractComponent
    index::Int
    name::String
    std_type::String
    hv_bus::Int
    mv_bus::Int
    lv_bus::Int
    sn_hv_mva::Float64
    sn_mv_mva::Float64
    sn_lv_mva::Float64
    vn_hv_kv::Float64
    vn_mv_kv::Float64
    vn_lv_kv::Float64
    vk_hv_percent::Float64
    vk_mv_percent::Float64
    vk_lv_percent::Float64
    vkr_hv_percent::Float64
    vkr_mv_percent::Float64
    vkr_lv_percent::Float64
    pfe_kw::Float64
    i0_percent::Float64
    shift_mv_degree::Float64
    shift_lv_degree::Float64
    
    # 分接头参数
    tap_side::String
    tap_neutral::Int
    tap_min::Int
    tap_max::Int
    tap_step_percent::Float64
    tap_step_degree::Float64
    tap_at_star_point::Bool
    tap_pos::Int
    
    in_service::Bool
    
    # 技术参数
    vector_group_hv_mv::String
    vector_group_hv_lv::String
    vector_group_mv_lv::String
    hv_connection::String
    mv_connection::String
    lv_connection::String
    thermal_capacity_mw::Float64
    cooling_type::String
    oil_volume_liters::Float64
    winding_material::String
    
    # 构造函数
    function Transformer3W(index, name, std_type, hv_bus, mv_bus, lv_bus, 
                          sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                          vn_hv_kv, vn_mv_kv, vn_lv_kv,
                          vk_hv_percent, vk_mv_percent, vk_lv_percent,
                          vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                          pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; kwargs...)
        self = new(index, name, std_type, hv_bus, mv_bus, lv_bus, 
                  sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                  vn_hv_kv, vn_mv_kv, vn_lv_kv,
                  vk_hv_percent, vk_mv_percent, vk_lv_percent,
                  vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                  pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end


# 开关结构
mutable struct Switch <: AbstractComponent
    index::Int
    name::String
    bus_from::Int
    bus_to::Int
    element_type::String  # l-线路, t-变压器, b-节点
    element_id::Int
    closed::Bool
    type::String  # CB-断路器, LS-负荷开关, DS-隔离开关
    z_ohm::Float64
    in_service::Bool
    
    # 构造函数
    function Switch(index, name, bus_from, bus_to, element_type, element_id, 
                   closed, type, z_ohm, in_service)
        return new(index, name, bus_from, bus_to, element_type, element_id, 
                  closed, type, z_ohm, in_service)
    end
end

# 开关结构
mutable struct HighVoltageCircuitBreaker <: AbstractComponent
    index::Int
    name::String
    bus_from::Int
    bus_to::Int
    element_type::String  # l-线路, t-变压器, b-节点
    element_id::Int
    closed::Bool
    type::String  # CB-断路器, LS-负荷开关, DS-隔离开关
    z_ohm::Float64
    in_service::Bool
    
    # 构造函数
    function HighVoltageCircuitBreaker(index, name, bus_from, bus_to, element_type, element_id, 
                   closed, type, z_ohm, in_service)
        return new(index, name, bus_from, bus_to, element_type, element_id, 
                  closed, type, z_ohm, in_service)
    end
end