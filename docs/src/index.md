# TSPF.jl

*Time Series Power Flow Analysis Tool*

TSPF.jl is a Julia package tool for distribution system simulation design. It provides AC/DC hybrid power flow calculation and time series power flow calculation functions for power systems that support the access of photovoltaic, energy storage and VSC equipment.

## Features

- **Complete Components Model**: AC system models, DC system models, Hybrid models and carbon emission models, etc.
- **Complete case structure**: JuliaPowerCase, JPC
- **Independent functional Module**: ComponentModel, Utils, PowerFlow, TimeSeriesPowerFlow
- **DOPF Fusion**: Using relaxed OPF to bridge time varying loads and generation profiles on the energy status of ESSs
- **Multi Control Mode Aggregation**: Support 7 VSC control modes
- **Aggregation of externally input and internally generated studies**: Support external input cases and internal generated cases

## Project Structure

TSPF.jl contains four main modules:

- **ComponentModel**: Power system component modeling
- **Utils**: General utility functions
- **PowerFlow**: AC/DC hybrid power flow calculation
- **TimeSeriesPowerFlow**: Time series power flow analysis

## Quick Start

```julia
# Add project path
push!(LOAD_PATH, "/path/to/TSPF")

# Import modules
using TSPF.PowerFlow
using TSPF.TimeSeriesPowerFlow

# Load test case
case = load_julia_power_data("data/test_case.xlsx")

# Run power flow calculation
results = runpf(case)
