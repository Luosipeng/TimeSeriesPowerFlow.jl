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
    # ... (additional fields)
end
```

This structure contains all parameters and variables for microgrid planning optimization.