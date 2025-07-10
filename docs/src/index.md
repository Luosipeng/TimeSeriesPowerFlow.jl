# TSPF.jl

*Time Series Power Flow Analysis Tool*

TSPF.jl is a Julia package developed by the HR-PES team of Xi'an Jiaotong University, which provides a dynamic power flow simulation framework for distribution systems incorporating various renewable resources.

## 1. Features
## Features

- **Integration of PF and OPF**: Using relaxed OPF to bridge time-varying loads and generation profiles with ESS status, integrating VSC power allocation into Power Flow calculations
- **Customizable Simulation Environment**: Comprehensive parameter customization including irradiance data, electricity prices, load profiles, and flexible network topology configuration
- **Comprehensive Renewable Energy Models**: Validated mathematical models for energy storage, photovoltaic generation, and VSCs for reliable distribution network analysis
- **Advanced VSC Control Framework**: Implementation of seven control modes (voltage, reactive power, droop control, etc.) with corresponding solution algorithms
- **Flexible Case Management**: Support for both externally imported and internally generated case studies with complete data structures

## Project Structure

TSPF.jl contains four main modules:

- **TSPF**: Power system component modeling
- **Utils**: General utility functions
- **PowerFlow**: AC/DC hybrid power flow calculation
- **TimeSeriesPowerFlow**: Time series power flow analysis

## 2. Installation

## Installation

You can install this package through Julia's package manager:

```julia
using Pkg
Pkg.add("TSPF")
```

Or, if you want to use the latest development version:

```julia
using Pkg
Pkg.add(url="https://github.com/Luosipeng/TSPF.jl.git")
```

## 3. Quick Start Example

## Quick Start
Here's a simple example showing how to use TSPF.jl to run the time series power flow:

```julia
# Add project path
push!(LOAD_PATH, "/path/to/TSPF")

# Import modules
using TSPF

file_path = joinpath(pwd(), "data", "test_case.xlsx")

case = load_julia_power_data(file_path)

# Topology processing
results, new_case = topology_analysis(case, output_file="topology_results.xlsx")

# View results
println("Found ", nrow(results["cycles"]), " cycles")
println("Network is divided into ", length(unique(results["nodes"].Partition)), " partitions")

jpc = JuliaPowerCase2Jpc(new_case)

opt = options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;

jpc_list, isolated = extract_islands_acdc(jpc)
n_islands = length(jpc_list)
println("Extracted $(n_islands) islands in total")

# Create results array
results_array = Vector{Any}(undef, n_islands)

println("Starting multi-threaded calculation...")
t_start = time()

# Use multi-threading to calculate power flow for each island
@threads for i in 1:n_islands
    results_array[i] = runhpf(jpc_list[i], opt)
end

t_end = time()
elapsed = t_end - t_start

# Construct result similar to @timed return
results = (value=results_array, time=elapsed)

# # Get voltage results for all nodes
voltage_results = get_bus_voltage_results_acdc(results, new_case)
```

## 4. Documentation Structure

## Documentation Structure

This documentation is divided into the following sections:

- **[Module](modules/componentmodel.md)**: Detailed usage and examples of various modules
  - [ComponentModel](modules/componentmodel.md): Detailed introduction to ComponentModel.jl
  - [Utils](modules/utils.md): Detailed introduction to Utils.jl
  - [PowerFlow](modules/powerflow.md): Detailed introduction to PowerFlow.jl
  - [TimeSeriesPowerFlow](modules/timeseriespowerflow.md): Detailed introduction to TimeSeriesPowerFlow.jl

- **[API](api/componentmodel.md)**: API document of various modules
  - [ComponentModel API](api/componentmodel.md): API document of ComponentModel.jl
  - [Utils API](api/utils.md): API document of Utils.jl
  - [PowerFlow API](api/powerflow.md): API document of PowerFlow.jl
  - [TimeSeriesPowerFlow API](api/timeseriespowerflow.md): API document of TimeSeriesPowerFlow.jl

- **[References](references.md)**: Citation information for TSPF.jl and acknowledgments of referenced works and dependencies

## 5. Contribution and License Information

## Contribution

Contributions to TSPF.jl are welcome! Please refer to the [Contribution Guidelines](https://github.com/Luosipeng/TSPF.jl/blob/master/CONTRIBUTING.md) for more information.

## License

TSPF.jl is licensed under the [MIT License](https://github.com/Luosipeng/TSPF.jl/blob/master/LICENSE).

## 6. Citation Information

## Citation

If you use TSPF.jl in your research, please cite:

```bibtex
@misc{TSPF.jl,
  author = {Sipeng Luo,Tianyang Zhao,Chaohong Bie},
  title = {TSPF.jl: a Julia package for distribution system dynamic power flow},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/Luosipeng/TSPF.jl}
}
```