# 设备碳排放结构
mutable struct EquipmentCarbon <: AbstractComponent
    index::Int
    element_type::String
    element_id::Int
    carbon_embodied_kgCO2e::Float64
    carbon_operational_kgCO2e_per_year::Float64
    lifetime_years::Int
    manufacturing_date::Date
    installation_date::Date
    recycling_rate_percent::Float64
    
    # 构造函数
    function EquipmentCarbon(index, element_type, element_id, carbon_embodied_kgCO2e,
                            carbon_operational_kgCO2e_per_year, lifetime_years,
                            manufacturing_date, installation_date, recycling_rate_percent)
        return new(index, element_type, element_id, carbon_embodied_kgCO2e,
                  carbon_operational_kgCO2e_per_year, lifetime_years,
                  manufacturing_date, installation_date, recycling_rate_percent)
    end
end