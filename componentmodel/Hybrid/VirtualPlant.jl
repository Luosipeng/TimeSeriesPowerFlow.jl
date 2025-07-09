# 虚拟电厂结构
mutable struct VirtualPowerPlant <: AbstractComponent
    index::Int
    name::String
    description::String
    control_area::String
    capacity_mw::Float64
    energy_mwh::Float64
    response_time_s::Float64
    ramp_rate_mw_per_min::Float64
    availability_percent::Float64
    operator::String
    in_service::Bool
    
    # 资源信息
    resource_type::String
    resource_id::Int
    capacity_share_percent::Float64
    control_priority::Int
    resource_response_time_s::Float64
    max_duration_h::Float64
    
    # 负荷信息
    timestamp::DateTime
    p_mw::Float64
    q_mvar::Float64
    flexibility_up_mw::Float64
    flexibility_down_mw::Float64
    flexibility_duration_h::Float64
    
    # 构造函数
    function VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                              response_time_s, ramp_rate_mw_per_min, availability_percent,
                              operator, in_service; kwargs...)
        self = new(index, name, description, control_area, capacity_mw, energy_mwh,
                  response_time_s, ramp_rate_mw_per_min, availability_percent,
                  operator, in_service)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end


# 虚拟电厂结构
mutable struct FlexLoad <: AbstractComponent
    index::Int
    name::String
    description::String
    control_area::String
    capacity_mw::Float64
    energy_mwh::Float64
    response_time_s::Float64
    ramp_rate_mw_per_min::Float64
    availability_percent::Float64
    operator::String
    in_service::Bool
    
    # 资源信息
    resource_type::String
    resource_id::Int
    capacity_share_percent::Float64
    control_priority::Int
    resource_response_time_s::Float64
    max_duration_h::Float64
    
    # 负荷信息
    timestamp::DateTime
    p_mw::Float64
    q_mvar::Float64
    flexibility_up_mw::Float64
    flexibility_down_mw::Float64
    flexibility_duration_h::Float64
    
    # 构造函数
    function VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                              response_time_s, ramp_rate_mw_per_min, availability_percent,
                              operator, in_service; kwargs...)
        self = new(index, name, description, control_area, capacity_mw, energy_mwh,
                  response_time_s, ramp_rate_mw_per_min, availability_percent,
                  operator, in_service)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end