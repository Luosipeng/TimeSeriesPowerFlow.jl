# # 两卷变压器结构
# mutable struct Transformer2W <: AbstractComponent
#     index::Int
#     name::String
#     std_type::String
#     hv_bus::Int
#     lv_bus::Int
#     sn_mva::Float64
#     vn_hv_kv::Float64
#     vn_lv_kv::Float64
#     vk_percent::Float64
#     vkr_percent::Float64
#     pfe_kw::Float64
#     i0_percent::Float64
#     shift_degree::Float64
    
#     # 分接头参数
#     tap_side::String
#     tap_neutral::Int
#     tap_min::Int
#     tap_max::Int
#     tap_step_percent::Float64
#     tap_step_degree::Float64
#     tap_pos::Int
#     tap_phase_shifter::Bool
    
#     parallel::Int
#     max_loading_percent::Float64
#     df::Float64
#     in_service::Bool
#     oltc::Bool
#     power_station_unit::Bool
    
#     # 技术参数
#     vector_group::String
#     hv_connection::String
#     lv_connection::String
#     thermal_capacity_mw::Float64
#     cooling_type::String
#     oil_volume_liters::Float64
#     winding_material::String
    
#     # 构造函数
#     function Transformer2W(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
#                           vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree; kwargs...)
#         self = new(index, name, std_type, hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,
#                   vk_percent, vkr_percent, pfe_kw, i0_percent, shift_degree)
        
#         # 设置其他参数
#         for (key, value) in kwargs
#             if hasfield(typeof(self), key)
#                 setfield!(self, key, value)
#             end
#         end
        
#         return self
#     end
# end

# # 三卷变压器结构
# mutable struct Transformer3W <: AbstractComponent
#     index::Int
#     name::String
#     std_type::String
#     hv_bus::Int
#     mv_bus::Int
#     lv_bus::Int
#     sn_hv_mva::Float64
#     sn_mv_mva::Float64
#     sn_lv_mva::Float64
#     vn_hv_kv::Float64
#     vn_mv_kv::Float64
#     vn_lv_kv::Float64
#     vk_hv_percent::Float64
#     vk_mv_percent::Float64
#     vk_lv_percent::Float64
#     vkr_hv_percent::Float64
#     vkr_mv_percent::Float64
#     vkr_lv_percent::Float64
#     pfe_kw::Float64
#     i0_percent::Float64
#     shift_mv_degree::Float64
#     shift_lv_degree::Float64
    
#     # 分接头参数
#     tap_side::String
#     tap_neutral::Int
#     tap_min::Int
#     tap_max::Int
#     tap_step_percent::Float64
#     tap_step_degree::Float64
#     tap_at_star_point::Bool
#     tap_pos::Int
    
#     in_service::Bool
    
#     # 技术参数
#     vector_group_hv_mv::String
#     vector_group_hv_lv::String
#     vector_group_mv_lv::String
#     hv_connection::String
#     mv_connection::String
#     lv_connection::String
#     thermal_capacity_mw::Float64
#     cooling_type::String
#     oil_volume_liters::Float64
#     winding_material::String
    
#     # 构造函数
#     function Transformer3W(index, name, std_type, hv_bus, mv_bus, lv_bus, 
#                           sn_hv_mva, sn_mv_mva, sn_lv_mva, 
#                           vn_hv_kv, vn_mv_kv, vn_lv_kv,
#                           vk_hv_percent, vk_mv_percent, vk_lv_percent,
#                           vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
#                           pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; kwargs...)
#         self = new(index, name, std_type, hv_bus, mv_bus, lv_bus, 
#                   sn_hv_mva, sn_mv_mva, sn_lv_mva, 
#                   vn_hv_kv, vn_mv_kv, vn_lv_kv,
#                   vk_hv_percent, vk_mv_percent, vk_lv_percent,
#                   vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
#                   pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree)
        
#         # 设置其他参数
#         for (key, value) in kwargs
#             if hasfield(typeof(self), key)
#                 setfield!(self, key, value)
#             end
#         end
        
#         return self
#     end
# end