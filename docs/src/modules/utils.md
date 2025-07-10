# Utils Module

This module provides utility functions for power system analysis, including data conversion, network topology analysis, and various helper functions for working with power system models.

## Data Conversion Functions

### External to Internal Indexing

```julia
ext2int(bus::Matrix{Float64}, gen::Matrix{Float64}, branch::Matrix{Float64})
```

Converts external bus numbers to consecutive internal bus numbers.

- **Arguments**:
  - `bus`: Bus data matrix
  - `gen`: Generator data matrix
  - `branch`: Branch data matrix
- **Returns**: Converted matrices and mapping between external and internal indices

### Internal to External Indexing

```julia
int2ext(i2e::Vector{Int}, bus::Matrix{Float64}, gen::Matrix{Float64}, branch::Matrix{Float64})
```

Converts internal consecutive bus numbers back to original external bus numbers.

- **Arguments**:
  - `i2e`: Mapping from internal to external bus numbers
  - `bus`: Bus data matrix with internal indexing
  - `gen`: Generator data matrix with internal indexing
  - `branch`: Branch data matrix with internal indexing
- **Returns**: Matrices with external indexing restored

## Network Topology Analysis

### Island Detection

```julia
find_islands(jpc::JPC)
```

Identifies electrically isolated islands in an AC power system.

- **Arguments**:
  - `jpc`: Power system case data structure
- **Returns**:
  - `islands`: List of islands, each containing connected bus indices
  - `isolated_buses`: List of isolated buses

```julia
extract_islands(jpc::JPC)
```

Extracts electrically isolated islands from a power system case.

- **Arguments**:
  - `jpc`: Power system case data structure
- **Returns**:
  - `Vector{JPC}`: List of JPC objects, each representing a separate island

### Node Mapping

```julia
create_node_mapping(case::JuliaPowerCase)
```

Creates numbering mapping for nodes and detects duplicate nodes.

- **Arguments**:
  - `case`: JuliaPowerCase object
- **Returns**: A dictionary mapping node names to IDs

```julia
resolve_node_mapping(node_id, node_merge_map)
```

Resolves the final mapping of a node by traversing through a node merge map.

- **Arguments**:
  - `node_id`: The initial node ID to resolve
  - `node_merge_map`: A dictionary mapping source nodes to destination nodes
- **Returns**: The final resolved node ID after following all mappings

### Virtual Node Handling

```julia
merge_virtual_nodes(case::JuliaPowerCase)
```

Merges virtual nodes in a power system case and updates all connected elements.

- **Arguments**:
  - `case`: JuliaPowerCase object
- **Returns**: Updated case with virtual nodes merged

## Data Extraction

```julia
extract_data(sheet_name::String, xf)
```

Extracts data from an Excel worksheet into a matrix.

- **Arguments**:
  - `sheet_name`: Name of the worksheet
  - `xf`: Excel file object
- **Returns**: Matrix containing the extracted data

## Data Format Conversion

### MATLAB to Julia Conversion

```julia
parse_matlab_case_file(filepath)
```

Parses a MATPOWER case file in MATLAB format and converts it to Julia data structures.

- **Arguments**:
  - `filepath`: Path to the MATLAB case file
- **Returns**: Parsed power system data in Julia format

### JuliaPowerCase to JPC Conversion

```julia
JuliaPowerCase2Jpc(case::JuliaPowerCase)
```

Converts a JuliaPowerCase object to a JPC (Julia Power Case) object.

- **Arguments**:
  - `case`: JuliaPowerCase object
- **Returns**: JPC object

```julia
JuliaPowerCase2Jpc_3ph(case::JuliaPowerCase)
```

Converts a JuliaPowerCase object to a three-phase JPC_3ph object.

- **Arguments**:
  - `case`: JuliaPowerCase object
- **Returns**:
  - `jpc_3ph`: The converted three-phase power system model
  - Additional parameters for three-phase analysis

## ETAP Data Import

```julia
load_julia_power_data(file_path::String)
```

Loads power system data from an Excel file and converts it to a JuliaPowerCase structure.

- **Arguments**:
  - `file_path`: Path to the Excel file
- **Returns**: JuliaPowerCase structure containing all power system components

## Indexing Constants

The `idx.jl` file provides constants for indexing into power system data matrices:

```julia
function idx_bus()
    # define bus types
    PQ      = 1;
    PV      = 2;
    REF     = 3;
    NONE    = 4;
    
    # define the indices
    BUS_I       = 1;    # bus number (1 to 29997)
    BUS_TYPE    = 2;    # bus type (1-PQ, 2-PV, 3-ref, 4-isolated)
    PD          = 3;    # Pd, real power demand (MW)
    QD          = 4;    # Qd, reactive power demand (MVAr)
    GS          = 5;    # Gs, shunt conductance (MW at V = 1.0 p.u.)
    BS          = 6;    # Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
    BUS_AREA    = 7;    # area number, 1-100
    VM          = 8;    # Vm, voltage magnitude (p.u.)
    VA          = 9;    # Va, voltage angle (degrees)
    BASE_KV     = 10;   # baseKV, base voltage (kV)
    ZONE        = 11;   # zone, loss zone (1-999)
    VMAX        = 12;   # maxVm, maximum voltage magnitude (p.u.)
    VMIN        = 13;   # minVm, minimum voltage magnitude (p.u.)
    
    return PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN
end
```

Similar indexing functions are provided for generators, branches, and other power system components.

## Types

The `Types.jl` file defines data structures for power system analysis, including:

```julia
mutable struct MicrogridPlanningProblem
    # Constants
    bigM::Float64
    
    # Planning parameters
    nPV::Int
    nBESS::Int
    DeltaPV::Float64
    DeltaBESS::Float64
    
    # Technical parameters
    soc0::Float64
    soc_max::Float64
    soc_min::Float64
    eta_ug::Float64
    eff_ch::Float64
    eff_dc::Float64
    
    # Economic parameters
    pv_unit_cost::Float64
    ess_unit_cost_e::Float64
    ess_unit_cost_p::Float64
    ug_unit_cost::Float64
    carbon_permit_price::Float64
    carbon_emission_coeff::Float64
    npv::Float64

    # Reliability parameters
    ug_prob::Float64
    mg_prob::Float64
    power_requirements::Float64
    energy_requirements::Float64
    pv_power_factor::Float64
    pv_energy_factor::Float64
    
    # Carbon emission limit
    carbon_limit::Float64
    # Uncertainty parameters
    ns::Int
    ws::Vector{Float64}
    T::Int
    ppv_un::Int
    pl::Int
    fl::Int
    ug_status::Int
    ele_price::Int
    carbon_price::Int
    carbon_emission_ug::Int
    nu::Int
    NU::Int
    
    # Indices
    IPV::Int
    IESS::Int
    IPESS::Int
    PUG::Int
    CARBON_PERMIT::Int
    CARBON::Int
    NX::Int
    
    # First stage constraints
    lx::Vector{Float64}
    ux::Vector{Float64}
    vtypex::Vector{Char}
    Aeq::SparseMatrixCSC{Float64, Int}
    beq::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    b::Vector{Float64}
    
    # Objective function
    c::Vector{Float64}

    # Second stage indices
    pug::Int
    pug_sold::Int
    ppv::Int
    pess_ch::Int
    pess_dc::Int
    eess::Int
    iess_ch::Int
    iess_dc::Int
    wpv::Int
    wpess_dc::Int
    wpess_ch::Int
    plc::Int
    pf::Int
    ny::Int
    NY::Int
end
```

```julia
mutable struct JuliaPowerCase
    version::String
    baseMVA::Float32
    basef::Float32
    
    # AC Network Components
    busesAC::Vector{Bus}                # Matches busAC in JPC
    branchesAC::Vector{Line}            # Matches branchAC in JPC
    loadsAC::Vector{Load}               # Matches loadAC in JPC
    loadsAC_flex::Vector{FlexLoad}      # Matches loadAC_flex in JPC
    loadsAC_asymm::Vector{AsymmetricLoad}  # Matches loadAC_asymm in JPC
    # branches3ph::Vector{ThreePhaseBranch}  # Matches branch3ph in JPC
    gensAC::Vector{Generator}            # Matches genAC in JPC
    
    # DC Network Components
    busesDC::Vector{BusDC}              # Matches busDC in JPC
    branchesDC::Vector{LineDC}          # Matches branchDC in JPC
    loadsDC::Vector{LoadDC}          # Matches loadDC in JPC
    
    # Distributed Energy Resources
    sgensAC::Vector{StaticGenerator}     # Matches sgenAC in JPC
    storages::Vector{Storage}            # Matches storage in JPC
    storageetap::Vector{Storageetap}      # Matches storageetap in JPC
    sgensDC::Vector{StaticGeneratorDC}   # Matches sgenDC in JPC
    pvarray::Vector{PVArray}          # Matches pv_array in JPC
    ACPVSystems::Vector{ACPVSystem}  # Matches AC PV systems in JPC
    
    # Special Components
    converters::Vector{Converter}        # Matches converter in JPC
    ext_grids::Vector{ExternalGrid}      # Matches ext_grid in JPC
    hvcbs::Vector{HighVoltageCircuitBreaker}  # Matches hvcb in JPC
    microgrids::Vector{Microgrid}        # Matches microgrid in JPC
    
    # Original components without direct JPC equivalents (preserved)
    transformers_2w::Vector{Transformer2W}
    transformers_3w::Vector{Transformer3W}
    transformers_2w_etap::Vector{Transformer2Wetap}
    charging_stations::Vector{ChargingStation}
    chargers::Vector{Charger}
    ev_aggregators::Vector{EVAggregator}
    virtual_power_plants::Vector{VirtualPowerPlant}
    carbon_time_series::Vector{CarbonTimeSeries}
    equipment_carbon::Vector{EquipmentCarbon}
    
    # Lookup dictionaries
    bus_name_to_id::Dict{String, Int}
    busdc_name_to_id::Dict{String, Int}
    zone_to_id::Dict{String, Int}
    area_to_id::Dict{String, Int}
end
```

```julia
mutable struct JPC
    version::String
    baseMVA::Float64
    success::Bool
    iterationsAC::Int
    iterationsDC::Int
    
    # AC Network Components
    busAC::Array{Float64,2}         # Bus data 
    genAC::Array{Float64,2}         # Generator data 
    branchAC::Array{Float64,2}      # Branch data
    loadAC::Array{Float64,2}        # Load data
    loadAC_flex::Array{Float64,2}   # Flexible load data 
    loadAC_asymm::Array{Float64,2}  # Asymmetric load data
    branch3ph::Array{Float64,2}     # Three-phase branch data
    
    # DC Network Components
    busDC::Array{Float64,2}         # DC bus data 
    branchDC::Array{Float64,2}      # DC branch data 
    genDC::Array{Float64,2}         # DC generator data 
    loadDC::Array{Float64,2}        # DC load data
    
    # Distributed Energy Resources
    sgenAC::Array{Float64,2}        # Solar/PV generation data 
    storageetap::Array{Float64,2}  # Energy storage system data
    storage::Array{Float64,2}       # Energy storage system data 
    sgenDC::Array{Float64,2}        # DC-connected PV data 
    pv::Array{Float64,2}          # PV array data 
    pv_acsystem::Array{Float64,2}  # AC PV systems data 
    
    # Special Components
    converter::Array{Float64,2}     # AC/DC converter data 
    ext_grid::Array{Float64,2}      # External grid data 
    hvcb::Array{Float64,2}          # High voltage circuit breaker data 
    microgrid::Array{Float64,2}     # Microgrid data 
end
```

```julia
mutable struct JPC_3ph
    version::String
    baseMVA::Float32
    basef::Float32
    mode::String
    success::Bool  # Flag indicating if three-phase power flow converged
    iterations::Int  # Number of iterations for three-phase power flow
    
    # AC Network Components - Matrix representations of JuliaPowerCase components
    busAC_0::Array{Float64,2}  # Phase A bus data
    busAC_1::Array{Float64,2}  # Phase B bus data
    busAC_2::Array{Float64,2}  # Phase C bus data

    branchAC_0::Array{Float64,2}  # Phase A branch data
    branchAC_1::Array{Float64,2}  # Phase B branch data
    branchAC_2::Array{Float64,2}  # Phase C branch data

    loadAC_0::Array{Float64,2}  # Phase A load data
    loadAC_1::Array{Float64,2}  # Phase B load data
    loadAC_2::Array{Float64,2}  # Phase C load data

    genAC_0::Array{Float64,2}  # Phase A generator data
    genAC_1::Array{Float64,2}  # Phase B generator data
    genAC_2::Array{Float64,2}  # Phase C generator data

    storageAC::Array{Float64,2}  # Storage data
    
    # DC Network Components
    busDC::Array{Float64,2}  # DC bus data
    branchDC::Array{Float64,2}  # DC branch data
    loadDC::Array{Float64,2}  # DC load data
    genDC::Array{Float64,2}  # DC generator data
    storageDC::Array{Float64,2}  # DC storage data
    
    # Special Components
    ext_grid::Array{Float64,2}  # External grid data
    switche::Array{Float64,2}  # Switch data

    # Three-phase power flow results
    res_bus_3ph::Array{Float64,2}  # Bus results from three-phase power flow
    res_loadsAC_3ph::Array{Float64,2}  # Load results from three-phase power flow
    res_ext_grid_3ph::Array{Float64,2}  # External grid results from three-phase power flow
    
    # Lookup dictionaries
    bus_name_to_id::Dict{String, Int}
    zone_to_id::Dict{String, Int}
    area_to_id::Dict{String, Int}
end
```

This structure contains all parameters and variables for microgrid planning optimization.