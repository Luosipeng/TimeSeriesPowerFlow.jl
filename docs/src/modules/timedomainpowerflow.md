# TimeDomainPowerFlow Module

The TimeDomainPowerFlow module provides tools and functions for time-series power flow analysis in power systems, including dynamic dispatch, voltage profile analysis, and renewable integration.

## Module Overview

```julia
module TimeDomainPowerFlow

    using Printf
    using SparseArrays
    using Plots
    using Statistics
    using DataFrames
    using XLSX
    using Dates

    # Module exports and includes
    # ...
end
```

## Flow Chart

![TimeDomainPowerFlow flow chart](../assets/TimeDomainPowerFlow.png)

## Time Series Data Processing

### Reading Load Data

```julia
read_load_data(file_path)
```

Read load time series data from an Excel file.

**Arguments**:
- `file_path`: Path to the Excel file containing load data

**Returns**:
- `time_column`: Vector containing time values from the first column
- `time_str_column`: Vector containing time string representations from the second column
- `load_names`: Vector of load names extracted from column headers
- `data`: DataFrame containing the entire load dataset

This function reads load time series data from an Excel file into a DataFrame. It extracts the time values, time string representations, and the names of the loads from the column headers.

## Simulation Functions

### Run Single Day Simulation

```julia
run_single_day(old_jpc, opt, day_load_matrix, day_price_line, day_irradiance_line)
```

Run a day-ahead simulation for a hybrid AC-DC power system with renewable generation and energy storage.

**Arguments**:
- `old_jpc`: Original power system data structure
- `opt`: Power flow options
- `day_load_matrix`: Matrix containing hourly load data (hour, bus_id, active_power, reactive_power)
- `day_price_line`: Vector containing hourly electricity prices
- `day_irradiance_line`: Matrix containing hourly solar irradiance data

**Returns**:
- Simulation results including bus voltages, branch flows, and economic data

### Run Dynamic Dispatch

```julia
run_dynamic_dispatch(new_jpc, Cld_ac, Cld_dc, loadAC_PD, loadAC_QD, loadDC_PD, genAC_PG,
                    Cgen_ac, Cconv_ac, Cconv_dc, η_rec, η_inv, Cpv_ac, Cpv_dc,
                    pv_ac_p_mw_ratio, pv_ac_p_mw, pv_max_p_mw, pv_max_p_mw_ratio,
                    Cstorage_ac, ess_initial_soc, ess_max_soc, ess_min_soc,
                    ess_charge_efficiency, ess_discharge_efficiency, ess_p_mw,
                    ess_e_mwh, ess_max_p_mw, ess_min_p_mw, time_step_h, time_periods,
                    price_line, irradiance_line, opt)
```

Performs dynamic dispatch optimization for a hybrid AC-DC power system with renewable generation and energy storage.

**Arguments**:
- Power system parameters (jpc, loads, generators)
- Converter parameters (efficiency, costs)
- PV parameters (capacity, ratios)
- Storage parameters (SOC limits, efficiency, capacity)
- Time parameters and price/irradiance data
- Power flow options

**Returns**:
- Optimized dispatch results and power flow solutions

## Visualization Functions

### Plot Voltage Time Series

```julia
plot_voltage_time_series(results, bus_name, case, time_day, bus_type = "AC"; 
                         save_path = nothing, save_format = "pdf")
```

Plot voltage time series for a specified bus in a power system.

**Arguments**:
- `results`: Simulation results containing bus voltage data
- `bus_name`: Name of the bus to plot voltage for
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `bus_type`: Type of bus to analyze ("AC" or "DC")
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot

### Plot Power Flow Violations

```julia
plot_flow_violations(results, case, time_day, flow_limit = 3.0, plot_type = "summary", 
                    flow_direction = "max"; save_path = nothing, save_format = "pdf")
```

Plot power flow violations in power system branches.

**Arguments**:
- `results`: Simulation results containing branch flow data
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `flow_limit`: Power flow limit in MW (default: 3.0)
- `plot_type`: Type of plot to generate:
  - "summary": Overall statistics of violations
  - "worst": Shows the worst branches with violations
  - "all": Shows all branches with violations
- `flow_direction`: How to evaluate flow violations:
  - "max": Maximum absolute value of flow in either direction
  - "both": Same as "max"
  - "forward": Only check flow from from-bus to to-bus
  - "reverse": Only check flow from to-bus to from-bus
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot

### Plot Losses Time Series

```julia
plot_losses_time_series(results, case, time_day, plot_type = "total", loss_type = "active"; 
                        save_path = nothing, save_format = "pdf")
```

Plot power system losses time series for AC and DC branches.

**Arguments**:
- `results`: Simulation results containing branch flow data
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `plot_type`: Type of plot to generate:
  - "total": Overall system losses
  - "branch": Individual branch losses for major branches
- `loss_type`: Type of losses to display:
  - "active": Only active power losses (MW)
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot

### Plot Active Load Time Series

```julia
plot_PD_time_series(results, bus_name, case, time_day; 
                   save_path = nothing, save_format = "pdf")
```

Plot the time series of active load.

**Arguments**:
- `results`: Result dataset
- `bus_name`: Bus name
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot

### Record Voltage Violations

```julia
record_voltage_violation(results, bus_name, case, time_day, bus_type = "AC"; 
                        save_path = nothing, save_format = "pdf")
```

Analyze and visualize voltage violations for a specified bus in a power system.

**Arguments**:
- `results`: Simulation results containing bus voltage data
- `bus_name`: Name of the bus to analyze voltage violations for
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `bus_type`: Type of bus to analyze:
  - "AC": AC bus
  - "DC": DC bus
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot

## Utility Functions

### Renumber Hybrid System

```julia
renumber_hybrid_system(jpc)
```

Renumber buses and branches in a hybrid AC-DC power system to create a more organized numbering scheme.

**Arguments**:
- `jpc`: A structure containing the hybrid power system data with fields:
  - `busAC`: Matrix of AC bus data
  - `busDC`: Matrix of DC bus data
  - `branchAC`: Matrix of AC branch data
  - `branchDC`: Matrix of DC branch data
  - `convAC`: Matrix of AC converter data
  - `convDC`: Matrix of DC converter data

**Returns**:
- Updated `jpc` with renumbered components

## Testing

The module includes a test script (`test_time_domain_pf.jl`) that demonstrates the use of the TimeDomainPowerFlow module for analyzing a power system over a time period. The script:

1. Loads required packages
2. Reads power system data
3. Configures simulation parameters
4. Runs time-domain power flow analysis
5. Generates visualization of results

```julia
using Dates
using XLSX
using DataFrames
using Base.Filesystem
using HyDistFlow.TimeDomainPowerFlowlow

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