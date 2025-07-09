"""
    Load <: AbstractComponent

Structure representing a symmetric load in power systems.

# Fields
- `index::Int`: Unique identifier for the load
- `name::String`: Load name
- `bus::Int`: Bus number where the load is connected
- `p_mw::Float64`: Active power demand (MW)
- `q_mvar::Float64`: Reactive power demand (MVAr)
- `const_z_percent::Float64`: Percentage of constant impedance load model
- `const_i_percent::Float64`: Percentage of constant current load model
- `const_p_percent::Float64`: Percentage of constant power load model
- `scaling::Float64`: Scaling factor for the load
- `in_service::Bool`: Operating status flag
- `type::String`: Load connection type (wye, delta, etc.)
"""
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
    type::String  # wye, delta, etc.
    
    """
        Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
             const_p_percent, scaling, in_service, type)

    Create a new symmetric load instance.

    # Parameters
    - `index`: Unique identifier for the load
    - `name`: Load name
    - `bus`: Bus number where the load is connected
    - `p_mw`: Active power demand (MW)
    - `q_mvar`: Reactive power demand (MVAr)
    - `const_z_percent`: Percentage of constant impedance load model
    - `const_i_percent`: Percentage of constant current load model
    - `const_p_percent`: Percentage of constant power load model
    - `scaling`: Scaling factor for the load
    - `in_service`: Operating status flag
    - `type`: Load connection type (wye, delta, etc.)
    """
    function Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service, type)
        return new(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service, type)
    end
end

"""
    AsymmetricLoad <: AbstractComponent

Structure representing an asymmetric load in power systems.

# Fields
- `index::Int`: Unique identifier for the load
- `name::String`: Load name
- `bus::Int`: Bus number where the load is connected
- `p_mw::Float64`: Active power demand (MW)
- `q_mvar::Float64`: Reactive power demand (MVAr)
- `const_z_percent::Float64`: Percentage of constant impedance load model
- `const_i_percent::Float64`: Percentage of constant current load model
- `const_p_percent::Float64`: Percentage of constant power load model
- `scaling::Float64`: Scaling factor for the load
- `in_service::Bool`: Operating status flag
- `type::String`: Load connection type (wye, delta, etc.)
"""
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
    type::String  # wye, delta, etc.
    
    """
        Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
             const_p_percent, scaling, in_service, type)

    Create a new asymmetric load instance.

    # Parameters
    - `index`: Unique identifier for the load
    - `name`: Load name
    - `bus`: Bus number where the load is connected
    - `p_mw`: Active power demand (MW)
    - `q_mvar`: Reactive power demand (MVAr)
    - `const_z_percent`: Percentage of constant impedance load model
    - `const_i_percent`: Percentage of constant current load model
    - `const_p_percent`: Percentage of constant power load model
    - `scaling`: Scaling factor for the load
    - `in_service`: Operating status flag
    - `type`: Load connection type (wye, delta, etc.)
    """
    function Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                 const_p_percent, scaling, in_service, type)
        return new(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent, 
                  const_p_percent, scaling, in_service, type)
    end
end
