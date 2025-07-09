# Symetric Load Model
mutable struct LoadDC <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    const_z_percent::Float64
    const_i_percent::Float64
    const_p_percent::Float64
    scaling::Float64
    in_service::Bool

    # 构造函数
    function LoadDC(index, name, bus, p_mw, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service)
        return new(index, name, bus, p_mw, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service)
    end
end

mutable struct PVArray <: AbstractComponent
    index::Int
    name::String
    bus::Int
    numpanelseries::Int
    numpanelparallel::Int
    vmpp::Float64
    impp::Float64
    voc::Float64
    isc::Float64
    pvanumcells::Int
    temperature::Float64
    irradiance::Float64
    α_isc::Float64
    β_voc::Float64
    in_service::Bool

    # 构造函数
    function PVArray(index, name, bus, numpanelseries, numpanelparallel, vmpp, impp,
                 voc, isc, pvanumcells, temperature, irradiance, α_isc,
                    β_voc, in_service)
        return new(index, name, bus, numpanelseries, numpanelparallel, vmpp, impp,
                  voc, isc, pvanumcells, temperature, irradiance, α_isc,
                  β_voc, in_service)
    end
end