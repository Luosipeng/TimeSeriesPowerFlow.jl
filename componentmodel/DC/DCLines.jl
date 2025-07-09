# 直流线路结构
mutable struct LineDC <: AbstractComponent
    index::Int
    name::String
    from_bus::Int
    to_bus::Int
    length_km::Float64
    r_ohm_per_km::Float64
    x_ohm_per_km::Float64
    g_us_per_km::Float64
    max_i_ka::Float64
    type::String
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
    function LineDC(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km, 
                   g_us_per_km, max_i_ka, type, max_loading_percent, 
                   parallel, df, in_service)
        return new(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                  g_us_per_km, max_i_ka, type, max_loading_percent,
                  parallel, df, in_service)
    end
end