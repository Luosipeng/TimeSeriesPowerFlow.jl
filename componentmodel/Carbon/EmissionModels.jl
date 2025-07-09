"""
    EquipmentCarbon <: AbstractComponent

Represents carbon emission data for power system equipment.

# Fields
- `index::Int`: Unique identifier for the equipment carbon record
- `element_type::String`: Type of equipment (e.g., "transformer", "line", "generator")
- `element_id::Int`: Identifier of the specific equipment
- `carbon_embodied_kgCO2e::Float64`: Embodied carbon in the equipment in kg CO2 equivalent
- `carbon_operational_kgCO2e_per_year::Float64`: Annual operational carbon emissions in kg CO2 equivalent per year
- `lifetime_years::Int`: Expected lifetime of the equipment in years
- `manufacturing_date::Date`: Date when the equipment was manufactured
- `installation_date::Date`: Date when the equipment was installed
- `recycling_rate_percent::Float64`: Percentage of materials that can be recycled at end of life
"""
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
    
    # Constructor
    """
        EquipmentCarbon(index, element_type, element_id, carbon_embodied_kgCO2e,
                        carbon_operational_kgCO2e_per_year, lifetime_years,
                        manufacturing_date, installation_date, recycling_rate_percent)

    Create a new equipment carbon emission record.

    # Parameters
    - `index`: Unique identifier for the equipment carbon record
    - `element_type`: Type of equipment (e.g., "transformer", "line", "generator")
    - `element_id`: Identifier of the specific equipment
    - `carbon_embodied_kgCO2e`: Embodied carbon in the equipment in kg CO2 equivalent
    - `carbon_operational_kgCO2e_per_year`: Annual operational carbon emissions in kg CO2 equivalent per year
    - `lifetime_years`: Expected lifetime of the equipment in years
    - `manufacturing_date`: Date when the equipment was manufactured
    - `installation_date`: Date when the equipment was installed
    - `recycling_rate_percent`: Percentage of materials that can be recycled at end of life
    """
    function EquipmentCarbon(index, element_type, element_id, carbon_embodied_kgCO2e,
                            carbon_operational_kgCO2e_per_year, lifetime_years,
                            manufacturing_date, installation_date, recycling_rate_percent)
        return new(index, element_type, element_id, carbon_embodied_kgCO2e,
                  carbon_operational_kgCO2e_per_year, lifetime_years,
                  manufacturing_date, installation_date, recycling_rate_percent)
    end
end
