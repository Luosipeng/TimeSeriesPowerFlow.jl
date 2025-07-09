# 碳排放时间序列结构
mutable struct CarbonTimeSeries <: AbstractComponent
    index::Int
    timestamp::DateTime
    grid_carbon_intensity_kgCO2e_per_MWh::Float64
    renewable_generation_carbon_intensity_kgCO2e_per_MWh::Float64
    storage_carbon_intensity_kgCO2e_per_MWh::Float64
    
    # 构造函数
    function CarbonTimeSeries(index, timestamp, grid_carbon_intensity_kgCO2e_per_MWh,
                             renewable_generation_carbon_intensity_kgCO2e_per_MWh,
                             storage_carbon_intensity_kgCO2e_per_MWh)
        return new(index, timestamp, grid_carbon_intensity_kgCO2e_per_MWh,
                  renewable_generation_carbon_intensity_kgCO2e_per_MWh,
                  storage_carbon_intensity_kgCO2e_per_MWh)
    end
end

# 碳排放情景结构
mutable struct CarbonScenario <: AbstractComponent
    index::Int
    name::String
    description::String
    year::Int
    grid_carbon_intensity_kgCO2e_per_MWh::Float64
    renewable_penetration_percent::Float64
    ev_adoption_percent::Float64
    storage_adoption_percent::Float64
    
    # 构造函数
    function CarbonScenario(index, name, description, year, grid_carbon_intensity_kgCO2e_per_MWh,
                           renewable_penetration_percent, ev_adoption_percent, storage_adoption_percent)
        return new(index, name, description, year, grid_carbon_intensity_kgCO2e_per_MWh,
                  renewable_penetration_percent, ev_adoption_percent, storage_adoption_percent)
    end
end
