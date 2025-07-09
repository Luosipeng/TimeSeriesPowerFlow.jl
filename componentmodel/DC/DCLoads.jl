"""
    LoadDC <: AbstractComponent

Represents a symmetric load model in a DC power flow analysis.

# Fields
- `index::Int`: Unique identifier for the load
- `name::String`: Name of the load
- `bus::Int`: Bus identifier where the load is connected
- `p_mw::Float64`: Active power consumption in MW
- `const_z_percent::Float64`: Percentage of load with constant impedance characteristic
- `const_i_percent::Float64`: Percentage of load with constant current characteristic
- `const_p_percent::Float64`: Percentage of load with constant power characteristic
- `scaling::Float64`: Scaling factor for the load
- `in_service::Bool`: Operational status (true if in service)
"""
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

    # Constructor
    """
        LoadDC(index, name, bus, p_mw, const_z_percent, const_i_percent, 
              const_p_percent, scaling, in_service)

    Create a new symmetric load model for DC power flow analysis.

    # Parameters
    - `index`: Unique identifier for the load
    - `name`: Name of the load
    - `bus`: Bus identifier where the load is connected
    - `p_mw`: Active power consumption in MW
    - `const_z_percent`: Percentage of load with constant impedance characteristic
    - `const_i_percent`: Percentage of load with constant current characteristic
    - `const_p_percent`: Percentage of load with constant power characteristic
    - `scaling`: Scaling factor for the load
    - `in_service`: Operational status (true if in service)
    """
    function LoadDC(index, name, bus, p_mw, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service)
        return new(index, name, bus, p_mw, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service)
    end
end

"""
    PVArray <: AbstractComponent

Represents a photovoltaic array in a power system.

# Fields
- `index::Int`: Unique identifier for the PV array
- `name::String`: Name of the PV array
- `bus::Int`: Bus identifier where the PV array is connected
- `numpanelseries::Int`: Number of panels connected in series
- `numpanelparallel::Int`: Number of panel strings connected in parallel
- `vmpp::Float64`: Voltage at maximum power point in volts
- `impp::Float64`: Current at maximum power point in amperes
- `voc::Float64`: Open circuit voltage in volts
- `isc::Float64`: Short circuit current in amperes
- `pvanumcells::Int`: Number of cells in each PV panel
- `temperature::Float64`: Operating temperature in degrees Celsius
- `irradiance::Float64`: Solar irradiance in W/m²
- `α_isc::Float64`: Temperature coefficient of short circuit current in %/°C
- `β_voc::Float64`: Temperature coefficient of open circuit voltage in %/°C
- `in_service::Bool`: Operational status (true if in service)
"""
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

    # Constructor
    """
        PVArray(index, name, bus, numpanelseries, numpanelparallel, vmpp, impp,
               voc, isc, pvanumcells, temperature, irradiance, α_isc,
               β_voc, in_service)

    Create a new photovoltaic array model.

    # Parameters
    - `index`: Unique identifier for the PV array
    - `name`: Name of the PV array
    - `bus`: Bus identifier where the PV array is connected
    - `numpanelseries`: Number of panels connected in series
    - `numpanelparallel`: Number of panel strings connected in parallel
    - `vmpp`: Voltage at maximum power point in volts
    - `impp`: Current at maximum power point in amperes
    - `voc`: Open circuit voltage in volts
    - `isc`: Short circuit current in amperes
    - `pvanumcells`: Number of cells in each PV panel
    - `temperature`: Operating temperature in degrees Celsius
    - `irradiance`: Solar irradiance in W/m²
    - `α_isc`: Temperature coefficient of short circuit current in %/°C
    - `β_voc`: Temperature coefficient of open circuit voltage in %/°C
    - `in_service`: Operational status (true if in service)
    """
    function PVArray(index, name, bus, numpanelseries, numpanelparallel, vmpp, impp,
                 voc, isc, pvanumcells, temperature, irradiance, α_isc,
                    β_voc, in_service)
        return new(index, name, bus, numpanelseries, numpanelparallel, vmpp, impp,
                  voc, isc, pvanumcells, temperature, irradiance, α_isc,
                  β_voc, in_service)
    end
end
