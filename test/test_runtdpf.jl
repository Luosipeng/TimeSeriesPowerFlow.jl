using  Dates
using  XLSX
using  DataFrames
using  Base.Threads
using  HyDistFlow

file_path = joinpath(pwd(), "data", "test_case.xlsx")

load_path = joinpath(pwd(), "data", "load.xlsx")  # Load time series data for loads
price_path = joinpath(pwd(), "data", "price.xlsx")  # Load
irradiance_path = joinpath(pwd(), "data", "irradiance.xlsx")  # Load solar irradiance time series data

storage_profile = joinpath(pwd(), "data", "storage_profile.xlsx")  # Load storage profile data

# Load power system data and time series data
case = load_julia_power_data(file_path)  # Load power system case data
time_column, time_str_column, load_names, data = read_load_data(load_path)  # Read load time series data
time_column, time_str_column, price_profiles = read_price_data(price_path)  # Read price time series data
time_column, time_str_column, irradiance_profiles = read_irradiance_data(irradiance_path)  # Read solar irradiance time series data
time_column, time_str_column, storage_profiles = read_storage_profile_data(storage_profile)  # Read storage profile data

# Topology processing
# Analyze network topology and save results to Excel file
results, new_case = topology_analysis(case, output_file="topology_results.xlsx")

# Clear existing storage data and add a new battery storage system
empty!(new_case.storageetap)
push!(new_case.storages, Storage(1, "Battery_ESS_1", 3, 0.75, 1.5, 0.3, 0.05, 0.95, 0.9, true, "lithium_ion", true))
# Parameters: id, name, bus, power_rating, energy_capacity, initial_SOC, self_discharge_rate, max_SOC, min_SOC, is_available, type, can_discharge

# Set control mode for converters to Droop_Udc_Us (voltage droop control)
new_case.converters[3].control_mode = "Droop_Udc_Us"
new_case.converters[2].control_mode = "Droop_Udc_Us"
new_case.converters[1].control_mode = "Droop_Udc_Us"

# Configure power flow calculation options
opt = options()  # Initialize with default settings
opt["PF"]["NR_ALG"] = "bicgstab";    # Set Newton-Raphson algorithm to BiCGSTAB (Biconjugate Gradient Stabilized method)
opt["PF"]["ENFORCE_Q_LIMS"] = 0;     # Disable enforcement of reactive power limits
opt["PF"]["DC_PREPROCESS"] = 0;      # Enable DC power flow preprocessing

# Run time-series power flow calculation and measure execution time
@time results = runtdpf(new_case, data, load_names, price_profiles, irradiance_profiles, opt)

# get the voltage results for all nodes
plot_result, voltage_magnitude, voltage_angle = plot_voltage_time_series(results, "Bus_21", new_case, 366, "AC"; save_path="voltage_plot")
# plot_PD_time_series(results, "Bus_21", case, 30, "DC")
plot_result, violation_stats = record_voltage_violation(results, "Bus_21", new_case, 366, "AC", save_path="voltage_violation_analysis")

# # Plot and analyze system losses (active power losses in both AC and DC networks)
plot_result, total_ac_active_losses, total_ac_reactive_losses, total_dc_active_losses = plot_losses_time_series(results, new_case, 366, "total", "active", save_path="system_losses")

# # Analyze power flow violations in branches for time step 1
plot_result, violation_count, max_violation_percent, total_violation_severity, violation_details, branch_violation_stats = plot_flow_violations(results, case, 366, 1,"summary", "max", save_path="flow_violations")
