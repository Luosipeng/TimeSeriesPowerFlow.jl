# Symetric Load Model
mutable struct Load <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    const_z_percent::Float64
    const_i_percent::Float64
    const_p_percent::Float64
    scaling::Float64
    in_service::Bool
    type::String  # wye, delta等
    
    # 构造函数
    function Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service, type)
        return new(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service, type)
    end
end

# Symetric Load Model
mutable struct AsymmetricLoad <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    const_z_percent::Float64
    const_i_percent::Float64
    const_p_percent::Float64
    scaling::Float64
    in_service::Bool
    type::String  # wye, delta等
    
    # 构造函数
    function Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service, type)
        return new(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service, type)
    end
end

