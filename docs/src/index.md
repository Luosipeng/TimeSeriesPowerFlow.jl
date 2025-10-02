# HyDistFlow.jl

*Hybrid AC/DC Distribution System Analysis Tool*

HyDistFlow.jl is a Julia package developed by the HR-PES team of Xi'an Jiaotong University, which provides a static/dynamic power flow simulation framework for hybrid AC/DC distribution systems incorporating various renewable resources.

## 1. Features
### Features

- **Integration of PF and OPF**: Using relaxed OPF to bridge time-varying loads and generation profiles with ESS status, integrating VSC power allocation into Power Flow calculations
- **Customizable Simulation Environment**: Comprehensive parameter customization including irradiance data, electricity prices, load profiles, and flexible network topology configuration
- **Comprehensive Renewable Energy Models**: Validated mathematical models for energy storage, photovoltaic generation, and VSCs for reliable distribution network analysis
- **Advanced VSC Control Framework**: Implementation of seven control modes (voltage, reactive power, droop control, etc.) with corresponding solution algorithms
- **Flexible Case Management**: Support for both externally imported and internally generated case studies with complete data structures

### Project Structure

HyDistFlow.jl contains four main modules:

- **ComponentModel**: Power system component modeling
- **Utils**: General utility functions
- **PowerFlow**: AC/DC hybrid power flow calculation
- **TimeDomainPowerFlow**: Time series power flow analysis

## 2. Installation

### Installation

You can install this package through Julia's package manager:

```julia
using Pkg
Pkg.add("HyDistFlow")
```

Or, if you want to use the latest development version:

```julia
using Pkg
Pkg.add(url="https://github.com/Luosipeng/HyDistFlow.jl.git")
```

## 3. Quick Start Example

### Quick Start
Here's a simple example showing how to use HyDistFlow.jl to run the time series power flow:

```julia
# Add project path
push!(LOAD_PATH, "/path/to/HyDistFlow")

# Import modules
using Dates
using XLSX
using DataFrames
using Base.Threads
using HyDistFlow

# Input Data
file_path = joinpath(pwd(), "data", "test_case.xlsx")
load_path = joinpath(pwd(), "data", "load.xlsx")  
price_path = joinpath(pwd(), "data", "price.xlsx")  
irradiance_path = joinpath(pwd(), "data", "irradiance.xlsx")  

# Process Data
case = load_julia_power_data(file_path)
time_column, time_str_column, load_names, data = read_load_data(load_path) 
time_column, time_str_column, price_profiles = read_price_data(price_path)  
time_column, time_str_column, irradiance_profiles = read_irradiance_data(irradiance_path) 

# Topology processing
results, new_case = topology_analysis(case, output_file="topology_results.xlsx")

# Clear existing storage data and add a new battery storage system
empty!(new_case.storageetap)
push!(new_case.storages, Storage(1, "Battery_ESS_1", 3, 0.75, 1.5, 0.3, 0.05, 0.95, 0.9, true, "lithium_ion", true))

# Set control mode for converters to Droop_Udc_Us (voltage droop control)
new_case.converters[3].control_mode = "Droop_Udc_Us"
new_case.converters[2].control_mode = "Droop_Udc_Us"
new_case.converters[1].control_mode = "Droop_Udc_Us"

opt = options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;

# Run time-series power flow calculation and measure execution time
@time results = runtdpf(new_case, data, load_names, price_profiles, irradiance_profiles, opt)

# # Get voltage results for all nodes
plot_result = plot_voltage_time_series(results, "Bus_21", new_case, 366, "AC"; save_path="voltage_plot")
```

## 4. Documentation Structure

### Documentation Structure

This documentation is divided into the following sections:

- **[Module](modules/componentmodel.md)**: Detailed usage and examples of various modules
  - [ComponentModel](modules/componentmodel.md): Detailed introduction to ComponentModel.jl
  - [Utils](modules/utils.md): Detailed introduction to Utils.jl
  - [PowerFlow](modules/powerflow.md): Detailed introduction to PowerFlow.jl
  - [TimeDomainPowerFlow](modules/timedomainpowerflow.md): Detailed introduction to TimeDomainPowerFlow.jl

- **[API](api/componentmodel.md)**: API document of various modules
  - [ComponentModel API](api/componentmodel.md): API document of ComponentModel.jl
  - [Utils API](api/utils.md): API document of Utils.jl
  - [PowerFlow API](api/powerflow.md): API document of PowerFlow.jl
  - [TimeDomainPowerFlow API](api/timedomainpowerflow.md): API document of TimeDomainPowerFlow.jl

- **[References](references.md)**: Citation information for HyDistFlow.jl and acknowledgments of referenced works and dependencies

## 5. Contribution and License Information

### Contribution

Contributions to HyDistFlow.jl are welcome! Please refer to the [Contribution Guidelines](https://github.com/Luosipeng/HyDistFlow.jl/blob/master/CONTRIBUTING.md) for more information.

### License

HyDistFlow.jl is licensed under the [MIT License](https://github.com/Luosipeng/HyDistFlow.jl/blob/master/LICENSE).

## 6. Citation Information

### Citation

If you use HyDistFlow.jl in your research, please cite:

```bibtex
@software{HyDistFlow2025,
  author = {Sipeng Luo,Tianyang Zhao,Zhaohong Bie},
  title = {HyDistFlow.jl: a Julia package for time series power flow analysis},
  year = {2025},
  url = {https://github.com/Luosipeng/HyDistFlow.jl}
}
```