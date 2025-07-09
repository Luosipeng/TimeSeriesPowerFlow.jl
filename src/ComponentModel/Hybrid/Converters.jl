"""
    Converter Structure
    
    Represents a power converter that connects AC and DC systems, with control capabilities
    and operational parameters.
"""
mutable struct Converter <: AbstractComponent
    index::Int
    name::String
    bus_ac::Int
    bus_dc::Int
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
    control_mode::String
    droop_kv::Float64
    in_service::Bool
    controllable::Bool
    
    # Constructor
    function Converter(index, name, bus_ac, bus_dc, p_mw, q_mvar, 
                      vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                      max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                      control_mode, droop_kv, in_service, controllable)
        return new(index, name, bus_ac, bus_dc, p_mw, q_mvar, 
                  vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                  max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                  control_mode, droop_kv, in_service, controllable)
    end
end
