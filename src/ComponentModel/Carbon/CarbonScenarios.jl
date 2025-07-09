"""
    CarbonTimeSeries <: AbstractComponent

Represents a time series of carbon emission data in power systems.

# Fields
- `index::Int`: Unique identifier for the time series record
- `timestamp::DateTime`: Time point of the carbon emission data
- `grid_carbon_intensity_kgCO2e_per_MWh::Float64`: Grid carbon intensity in kg CO2 equivalent per MWh
- `renewable_generation_carbon_intensity_kgCO2e_per_MWh::Float64`: Carbon intensity of renewable generation in kg CO2 equivalent per MWh
- `storage_carbon_intensity_kgCO2e_per_MWh::Float64`: Carbon intensity of energy storage in kg CO2 equivalent per MWh
"""
mutable struct CarbonTimeSeries <: AbstractComponent
    index::Int
    timestamp::DateTime
    grid_carbon_intensity_kgCO2e_per_MWh::Float64
    renewable_generation_carbon_intensity_kgCO2e_per_MWh::Float64
    storage_carbon_intensity_kgCO2e_per_MWh::Float64
    
    # Constructor
    """
        CarbonTimeSeries(index, timestamp, grid_carbon_intensity_kgCO2e_per_MWh,
                         renewable_generation_carbon_intensity_kgCO2e_per_MWh,
                         storage_carbon_intensity_kgCO2e_per_MWh)

    Create a new carbon emission time series record.

    # Parameters
    - `index`: Unique identifier for the time series record
    - `timestamp`: Time point of the carbon emission data
    - `grid_carbon_intensity_kgCO2e_per_MWh`: Grid carbon intensity in kg CO2 equivalent per MWh
    - `renewable_generation_carbon_intensity_kgCO2e_per_MWh`: Carbon intensity of renewable generation in kg CO2 equivalent per MWh
    - `storage_carbon_intensity_kgCO2e_per_MWh`: Carbon intensity of energy storage in kg CO2 equivalent per MWh
    """
    function CarbonTimeSeries(index, timestamp, grid_carbon_intensity_kgCO2e_per_MWh,
                             renewable_generation_carbon_intensity_kgCO2e_per_MWh,
                             storage_carbon_intensity_kgCO2e_per_MWh)
        return new(index, timestamp, grid_carbon_intensity_kgCO2e_per_MWh,
                  renewable_generation_carbon_intensity_kgCO2e_per_MWh,
                  storage_carbon_intensity_kgCO2e_per_MWh)
    end
end

"""
    CarbonScenario <: AbstractComponent

Represents a carbon emission scenario for power system analysis.

# Fields
- `index::Int`: Unique identifier for the scenario
- `name::String`: Scenario name
- `description::String`: Detailed description of the scenario
- `year::Int`: Target year for the scenario
- `grid_carbon_intensity_kgCO2e_per_MWh::Float64`: Average grid carbon intensity in kg CO2 equivalent per MWh
- `renewable_penetration_percent::Float64`: Percentage of renewable energy penetration
- `ev_adoption_percent::Float64`: Percentage of electric vehicle adoption
- `storage_adoption_percent::Float64`: Percentage of energy storage adoption
"""
mutable struct CarbonScenario <: AbstractComponent
    index::Int
    name::String
    description::String
    year::Int
    grid_carbon_intensity_kgCO2e_per_MWh::Float64
    renewable_penetration_percent::Float64
    ev_adoption_percent::Float64
    storage_adoption_percent::Float64
    
    # Constructor
    """
        CarbonScenario(index, name, description, year, grid_carbon_intensity_kgCO2e_per_MWh,
                       renewable_penetration_percent, ev_adoption_percent, storage_adoption_percent)

    Create a new carbon emission scenario.

    # Parameters
    - `index`: Unique identifier for the scenario
    - `name`: Scenario name
    - `description`: Detailed description of the scenario
    - `year`: Target year for the scenario
    - `grid_carbon_intensity_kgCO2e_per_MWh`: Average grid carbon intensity in kg CO2 equivalent per MWh
    - `renewable_penetration_percent`: Percentage of renewable energy penetration
    - `ev_adoption_percent`: Percentage of electric vehicle adoption
    - `storage_adoption_percent`: Percentage of energy storage adoption
    """
    function CarbonScenario(index, name, description, year, grid_carbon_intensity_kgCO2e_per_MWh,
                           renewable_penetration_percent, ev_adoption_percent, storage_adoption_percent)
        return new(index, name, description, year, grid_carbon_intensity_kgCO2e_per_MWh,
                  renewable_penetration_percent, ev_adoption_percent, storage_adoption_percent)
    end
end
