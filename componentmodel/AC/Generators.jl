# 静态发电机结构
mutable struct StaticGenerator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    k::Float64
    rx::Float64
    in_service::Bool
    type::String  # WP-风电, PV-光伏, CHP-热电联产等
    controllable::Bool
    
    # 构造函数
    function StaticGenerator(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                            max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable)
        return new(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                  max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable)
    end
end

# 常规发电机组结构
mutable struct Generator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    vm_pu::Float64
    sn_mva::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    vn_kv::Float64
    xdss_pu::Float64
    rdss_pu::Float64
    cos_phi::Float64
    controllable::Bool
    in_service::Bool
    type::String
    generator_type::String
    fuel_type::String
    
    # 启动和运行参数
    startup_time_cold_h::Float64
    startup_time_warm_h::Float64
    startup_time_hot_h::Float64
    min_up_time_h::Float64
    min_down_time_h::Float64
    ramp_up_rate_mw_per_min::Float64
    ramp_down_rate_mw_per_min::Float64
    
    # 效率和排放参数
    efficiency_percent::Float64
    heat_rate_mmbtu_per_mwh::Float64
    co2_emission_rate_kg_per_mwh::Float64
    
    # 构造函数
    function Generator(index, name, bus, p_mw, vm_pu, sn_mva, scaling, 
                      max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                      vn_kv, xdss_pu, rdss_pu, cos_phi, controllable, in_service,
                      type, generator_type, fuel_type; kwargs...)
        self = new(index, name, bus, p_mw, vm_pu, sn_mva, scaling, 
                  max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                  vn_kv, xdss_pu, rdss_pu, cos_phi, controllable, in_service,
                  type, generator_type, fuel_type)
        
        # 设置其他参数
        for (key, value) in kwargs
            if hasfield(typeof(self), key)
                setfield!(self, key, value)
            end
        end
        
        return self
    end
end

# 外部电网结构
mutable struct ExternalGrid <: AbstractComponent
    index::Int
    name::String
    bus::Int
    vm_pu::Float64
    va_degree::Float64
    in_service::Bool
    s_sc_max_mva::Float64
    s_sc_min_mva::Float64
    rx_max::Float64
    rx_min::Float64
    r0x0_max::Float64
    x0x_max::Float64
    controllable::Bool
    
    # 构造函数
    function ExternalGrid(index, name, bus, vm_pu, va_degree, in_service, 
                         s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable)
        return new(index, name, bus, vm_pu, va_degree, in_service, 
                  s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable)
    end
end

# 交流光伏系统
mutable struct ACPVSystem <: AbstractComponent
    index::Int
    name::String
    bus::Int
    # 逆变器参数
    p_mw::Float64
    q_mvar::Float64
    vm_ac_pu::Float64
    vm_dc_pu::Float64
    loss_percent::Float64
    loss_mw::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    # 光伏组件参数
    numpanelseries::Int
    numpanelparallel::Int
    vmpp::Float64
    impp::Float64
    voc::Float64
    isc::Float64
    pvanumcells::Int
    irradiance::Float64
    temperature::Float64
    α_isc::Float64
    β_voc::Float64
    control_mode::String
    controllable::Bool
    in_service::Bool
    
    # 构造函数
    function ACPVSystem(index, name, bus, p_mw, q_mvar, vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                        max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                        numpanelseries, numpanelparallel, vmpp, impp, voc, isc,
                        pvanumcells, irradiance, temperature, α_isc, β_voc,
                        control_mode, controllable, in_service)
        return new(index, name, bus, p_mw, q_mvar, vm_ac_pu, vm_dc_pu,
                  loss_percent, loss_mw, max_p_mw, min_p_mw,
                  max_q_mvar, min_q_mvar,
                  numpanelseries, numpanelparallel,
                  vmpp, impp, voc, isc,
                  pvanumcells,
                  irradiance, temperature,
                  α_isc, β_voc,
                  control_mode,
                  controllable,
                  in_service)
    end

end

# 储能结构
mutable struct StorageAC <: AbstractComponent
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
    function StorageAC(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                     min_soc, max_soc, efficiency, in_service, type, controllable)
        return new(index, name, bus, power_capacity_mw, energy_capacity_mwh, soc_init,
                  min_soc, max_soc, efficiency, in_service, type, controllable)
    end
end