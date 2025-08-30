"""
    read_load_data(file_path)

Read load time series data from an Excel file.

# Arguments
- `file_path`: Path to the Excel file containing load data

# Returns
- `time_column`: Vector containing time values from the first column
- `time_str_column`: Vector containing time string representations from the second column
- `load_names`: Vector of load names extracted from column headers
- `data`: DataFrame containing the entire load dataset

# Description
This function reads load time series data from an Excel file into a DataFrame.
It extracts the time values (first column), time string representations (second column),
and the names of the loads from the column headers. The function returns these
components separately along with the complete DataFrame for further processing.
"""
function read_load_data(file_path)
    # Read Excel file into DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # Extract time column and time string column
    time_column = data[!, 1]  # First column is time
    time_str_column = data[!, 2]  # Second column is time string
    
    # Get load names (column names)
    load_names = names(data)[3:end]
    
    return time_column, time_str_column, load_names, data
end

"""
    read_price_data(file_path)

Read electricity price time series data from an Excel file.

# Arguments
- `file_path`: Path to the Excel file containing price data

# Returns
- `time_column`: Vector containing time values from the first column
- `time_str_column`: Vector containing time string representations from the second column
- `data`: DataFrame containing the entire price dataset

# Description
This function reads electricity price time series data from an Excel file into a DataFrame.
It extracts the time values (first column) and time string representations (second column).
The function returns these components separately along with the complete DataFrame for further processing.
"""
function read_price_data(file_path)
    # Read Excel file into DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # Extract time column and time string column
    time_column = data[!, 1]  # First column is time
    time_str_column = data[!, 2]  # Second column is time string
    
    
    return time_column, time_str_column, data
end

"""
    read_irradiance_data(file_path)

Read solar irradiance time series data from an Excel file.

# Arguments
- `file_path`: Path to the Excel file containing irradiance data

# Returns
- `hour_column`: Vector containing hour values from the first column
- `time_str_column`: Vector containing time string representations from the second column
- `data`: DataFrame containing the entire irradiance dataset

# Description
This function reads solar irradiance time series data from an Excel file into a DataFrame.
It extracts the hour values (first column) and time string representations (second column).
The function returns these components separately along with the complete DataFrame for further processing.
"""
function read_irradiance_data(file_path)
    # Read Excel file into DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # Extract hour column
    hour_column = data[!, 1]  # First column is hour
    # Extract time string column
    time_str_column = data[!, 2]  # Second column is time string
    
    return hour_column, time_str_column, data
end


"""
    read_storage_profile_data(file_path)

Read storage profile time series data from an Excel file.

# Arguments
- `file_path`: Path to the Excel file containing storage profile data

# Returns
- `hour_column`: Vector containing hour values from the first column
- `time_str_column`: Vector containing time string representations from the second column
- `data`: DataFrame containing the entire storage profile dataset

# Description
This function reads storage profile time series data from an Excel file into a DataFrame.
It extracts the hour values (first column) and time string representations (second column).
The function returns these components separately along with the complete DataFrame for further processing.
"""
function read_storage_profile_data(file_path)
    # Read Excel file into DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # Extract hour column
    hour_column = data[!, 1]  # First column is hour
    # Extract time string column
    time_str_column = data[!, 2]  # Second column is time string
    
    return hour_column, time_str_column, data
end

"""
    create_time_series_loads(case::Utils.JuliaPowerCase, data, load_names, num_days)

Create time series load data for power system simulation based on input data.

# Arguments
- `case`: Power system case data structure
- `data`: DataFrame containing load time series data
- `load_names`: Vector of load names
- `num_days`: Number of days to process

# Returns
- Dictionary mapping day number to load matrix for that day

# Description
This function processes time series load data and creates daily load matrices for power system simulation.
For each day, it creates a matrix where each row represents a specific bus at a specific hour, with columns:
[hour, bus_id, active_power, reactive_power, const_z_percent, const_i_percent, const_p_percent]

The function performs the following steps:
1. Validates that the number of time points is divisible by 24 (hours per day)
2. Filters in-service loads from the power system case
3. Maps load names to column indices in the input data
4. For each day and hour:
   - Accumulates loads connected to the same bus
   - Calculates actual active and reactive power based on apparent power and power factor
   - Computes weighted load characteristics (ZIP model parameters)
   - Ensures load characteristic percentages sum to 1
5. Returns a dictionary where keys are day numbers and values are the corresponding load matrices

The function handles cases where load names in the data match or don't match those in the power system case,
maintaining power factor and appropriately scaling loads based on the input data.
"""
function create_time_series_loads(case::Utils.JuliaPowerCase, data, load_names, num_days)
    # Get number of time points
    num_timepoints = size(data, 1)
    hours_per_day = 24
    
    # Check if time points can be divided by 24
    if num_timepoints % hours_per_day != 0
        @warn "Number of time points ($num_timepoints) is not a multiple of 24, will only process complete days"
        num_days = div(num_timepoints, hours_per_day)  # Adjust to complete days
        num_timepoints = num_days * hours_per_day
    else
        num_days = div(num_timepoints, hours_per_day)
    end
    
    # Filter in-service loads
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # If no in-service loads, return empty result
    if isempty(in_service_loads)
        @warn "No in-service loads"
        return Dict{Int, Array{Float64, 2}}()
    end
    
    # Create a dictionary to store column indices for each load name in the data table
    load_column_indices = Dict{String, Int}()
    
    # Get column names of the data table as strings
    data_column_names = string.(names(data))
    
    # Find column indices for each load name in the data table
    for name in load_names
        if name in data_column_names
            col_idx = findfirst(x -> string(x) == name, names(data))
            if col_idx !== nothing
                load_column_indices[name] = col_idx
            end
        end
    end
    
    # Create a dictionary to store load time series data by day
    # Format: day_loads[day] = [timepoint, bus_id, p_mw, q_mvar, const_z, const_i, const_p]
    day_loads = Dict{Int, Array{Float64, 2}}()
    
    # Get all unique bus IDs in the system
    unique_bus_ids = unique([load.bus for load in in_service_loads])
    num_buses = length(unique_bus_ids)
    
    # Create mapping from bus ID to index
    bus_id_to_index = Dict(bus_id => i for (i, bus_id) in enumerate(unique_bus_ids))
    
    # Process by day
    for day in 1:num_days
        @info "Processing day $day load data, total $num_days days"
        
        # Calculate time point range for the day
        start_t = (day - 1) * hours_per_day + 1
        end_t = day * hours_per_day
        
        # Create load matrix for the day
        # One row per time point per bus, columns are [time index, bus ID, P, Q, const_z, const_i, const_p]
        day_load_matrix = zeros(num_buses * hours_per_day, 7)
        
        # Process each time point of the day
        for hour in 1:hours_per_day
            t = start_t + hour - 1  # Global time point index
            
            # Create a dictionary to accumulate loads connected to the same bus
            bus_load_sum = Dict{Int, Vector{Float64}}()
            
            # Initialize bus load characteristic accumulator
            # [p_mw, q_mvar, p_mw*const_z, q_mvar*const_z, p_mw*const_i, q_mvar*const_i, p_mw*const_p, q_mvar*const_p]
            for bus_id in unique_bus_ids
                bus_load_sum[bus_id] = zeros(8)
            end
            
            # Create a dictionary for each load to accumulate matching load values
            load_power_dict = Dict{String, Float64}()
            
            # Get load values for the specified time point from the data table
            for (name, col_idx) in load_column_indices
                load_power_dict[name] = data[t, col_idx]
            end
            
            # Process each load
            for load in in_service_loads
                load_name = load.name
                bus_id = load.bus
                
                # Find all loads matching the current load name
                matching_loads = [name for name in keys(load_power_dict) if startswith(name, load_name)]
                
                # Calculate actual active and reactive loads
                if !isempty(matching_loads)
                    # Found matching loads, accumulate all matching load apparent power values
                    total_apparent_power = sum(load_power_dict[name] for name in matching_loads)
                    
                    # Get original power values
                    old_p = load.p_mw
                    old_q = load.q_mvar
                    
                    # Calculate original power factor (keep unchanged)
                    if old_p != 0 || old_q != 0
                        # Calculate original apparent power and power factor
                        old_s = sqrt(old_p^2 + old_q^2)
                        if old_s > 0
                            pf = old_p / old_s  # Power factor = P/S
                            # Keep power factor unchanged, calculate active and reactive power based on new apparent power
                            new_s = total_apparent_power / 1000.0  # Assume data unit is kVA, convert to MVA
                            new_p = new_s * pf  # P = S * cos(φ)
                            new_q = new_s * sqrt(1 - pf^2)  # Q = S * sin(φ)
                            
                            # Keep sign of reactive power (inductive or capacitive)
                            if old_q < 0
                                new_q = -new_q
                            end
                        else
                            # If original apparent power is 0, assume power factor is 0.9
                            new_s = total_apparent_power / 1000.0
                            new_p = new_s * 0.9
                            new_q = new_s * sqrt(1 - 0.9^2)
                        end
                    else
                        # If original powers are both 0, assume power factor is 0.9
                        new_s = total_apparent_power / 1000.0
                        new_p = new_s * 0.9
                        new_q = new_s * sqrt(1 - 0.9^2)
                    end
                else
                    # For loads not in the data table, keep original values
                    new_p = load.p_mw * load.scaling
                    new_q = load.q_mvar * load.scaling
                end
                
                # Get load characteristics
                const_z = load.const_z_percent / 100.0
                const_i = load.const_i_percent / 100.0
                const_p = load.const_p_percent / 100.0
                
                # Accumulate to the corresponding bus load
                bus_load_sum[bus_id][1] += new_p  # Total active power
                bus_load_sum[bus_id][2] += new_q  # Total reactive power
                bus_load_sum[bus_id][3] += new_p * const_z  # Constant impedance active
                bus_load_sum[bus_id][4] += new_q * const_z  # Constant impedance reactive
                bus_load_sum[bus_id][5] += new_p * const_i  # Constant current active
                bus_load_sum[bus_id][6] += new_q * const_i  # Constant current reactive
                bus_load_sum[bus_id][7] += new_p * const_p  # Constant power active
                bus_load_sum[bus_id][8] += new_q * const_p  # Constant power reactive
            end
            
            # Fill load matrix, merge loads on the same bus
            for bus_id in unique_bus_ids
                idx = bus_id_to_index[bus_id]
                row_idx = (hour - 1) * num_buses + idx
                
                # Get accumulated load data
                load_data = bus_load_sum[bus_id]
                total_p = load_data[1]
                total_q = load_data[2]
                
                # Calculate merged load characteristic percentages
                const_z_percent = total_p > 0 ? load_data[3] / total_p : 0.0
                const_i_percent = total_p > 0 ? load_data[5] / total_p : 0.0
                const_p_percent = total_p > 0 ? load_data[7] / total_p : 0.0
                
                # Ensure percentages sum to 1
                sum_percent = const_z_percent + const_i_percent + const_p_percent
                if sum_percent > 0
                    const_z_percent /= sum_percent
                    const_i_percent /= sum_percent
                    const_p_percent /= sum_percent
                else
                    # Default to constant power load
                    const_p_percent = 1.0
                    const_z_percent = 0.0
                    const_i_percent = 0.0
                end
                
                # Fill the corresponding row of the load matrix
                day_load_matrix[row_idx, :] = [
                    hour,           # Time point index (hour)
                    bus_id,         # Bus ID
                    total_p,        # Active load (MW)
                    total_q,        # Reactive load (MVAr)
                    const_z_percent,  # Constant impedance load percentage
                    const_i_percent,  # Constant current load percentage
                    const_p_percent   # Constant power load percentage
                ]
            end
        end
        
        # Store the day's load data in the dictionary
        day_loads[day] = day_load_matrix
    end
    
    return day_loads
end

"""
    create_time_series_prices(price_profiles, num_days=nothing)

Create time series price data for power system economic analysis.

# Arguments
- `price_profiles`: Matrix containing electricity price data, where rows represent days and columns represent hours
- `num_days`: Optional parameter specifying the number of days to process (defaults to all available days)

# Returns
- Dictionary mapping day number to price matrix for that day

# Description
This function processes electricity price profiles and organizes them into daily price matrices.
For each day, it creates a matrix where each row represents a specific hour, with columns:
[hour, price]

The function performs the following steps:
1. Determines the total number of days available in the price profiles
2. Limits processing to the specified number of days if provided
3. For each day and hour:
   - Extracts the price value from the price profiles
   - Creates a matrix with hour and price information
4. Returns a dictionary where keys are day numbers and values are the corresponding price matrices

The price profiles are expected to have a structure where the first column contains date information
and subsequent columns (2-25) contain hourly price data.
"""
function create_time_series_prices(price_profiles, num_days=nothing)
    # Get number of days
    total_days = size(price_profiles, 1)
    hours_per_day = 24
    
    # If days are specified, ensure not exceeding available days
    if num_days === nothing
        num_days = total_days
    else
        num_days = min(num_days, total_days)
    end
    
    @info "Processing price data, total $num_days days"
    
    # Create a dictionary to store price time series data by day
    # Format: day_prices[day] = [hour, price]
    day_prices = Dict{Int, Array{Float64, 2}}()
    
    # Process by day
    for day in 1:num_days
        @info "Processing day $day price data, total $num_days days"
        
        # Create price matrix for the day [hour, price]
        day_price_matrix = zeros(hours_per_day, 2)
        
        # Process each hour of the day
        for hour in 1:hours_per_day
            # Column index (first column is date, so start from second column)
            col_idx = hour + 1
            
            # Get price
            price = price_profiles[day, col_idx]
            
            # Fill price matrix
            day_price_matrix[hour, :] = [hour, price]
        end
        
        # Store the day's price data in the dictionary
        day_prices[day] = day_price_matrix
    end
    
    return day_prices
end

"""
    create_time_series_irradiance(irradiance_profiles, num_days=nothing)

Create time series solar irradiance data for renewable energy simulation.

# Arguments
- `irradiance_profiles`: Matrix containing solar irradiance data, where rows represent days and columns represent hours
- `num_days`: Optional parameter specifying the number of days to process (defaults to all available days)

# Returns
- Dictionary mapping day number to irradiance matrix for that day

# Description
This function processes solar irradiance profiles and organizes them into daily irradiance matrices.
For each day, it creates a matrix where each row represents a specific hour, with columns:
[hour, irradiance]

The function performs the following steps:
1. Determines the total number of days available in the irradiance profiles
2. Limits processing to the specified number of days if provided
3. For each day and hour:
   - Extracts the irradiance value from the irradiance profiles
   - Creates a matrix with hour and irradiance information
4. Returns a dictionary where keys are day numbers and values are the corresponding irradiance matrices

The irradiance profiles are expected to have a structure where the first column contains date information
and subsequent columns (2-25) contain hourly irradiance data.
"""
function create_time_series_irradiance(irradiance_profiles, num_days=nothing)
    # Get number of days
    total_days = size(irradiance_profiles, 1)
    hours_per_day = 24
    
    # If days are specified, ensure not exceeding available days
    if num_days === nothing
        num_days = total_days
    else
        num_days = min(num_days, total_days)
    end
    
    @info "Processing irradiance data, total $num_days days"
    
    # Create a dictionary to store irradiance time series data by day
    # Format: day_irradiance[day] = [hour, irradiance]
    day_irradiance = Dict{Int, Array{Float64, 2}}()
    
    # Process by day
    for day in 1:num_days
        @info "Processing day $day irradiance data, total $num_days days"
        
        # Create irradiance matrix for the day [hour, irradiance]
        day_irradiance_matrix = zeros(hours_per_day, 2)
        
        # Process each hour of the day
        for hour in 1:hours_per_day
            # Column index (first column is date, so start from second column)
            col_idx = hour + 1
            
            # Get irradiance
            irradiance = irradiance_profiles[day, col_idx]
            
            # Fill irradiance matrix
            day_irradiance_matrix[hour, :] = [hour, irradiance]
        end
        
        # Store the day's irradiance data in the dictionary
        day_irradiance[day] = day_irradiance_matrix
    end
    
    return day_irradiance
end

"""
    create_time_series_storage_profile(storage_profiles, num_days=nothing)

Create time series storage profile data for energy system simulation.

# Arguments
- `storage_profiles`: Matrix containing storage profile data, where rows represent days and columns represent hours
- `num_days`: Optional parameter specifying the number of days to process (defaults to all available days)

# Returns
- Dictionary mapping day number to storage profile matrix for that day

# Description
This function processes storage profiles and organizes them into daily storage matrices.
For each day, it creates a matrix where each row represents a specific hour, with columns:
[hour, storage_value]

The function performs the following steps:
1. Determines the total number of days available in the storage profiles
2. Limits processing to the specified number of days if provided
3. For each day and hour:
   - Extracts the storage value from the storage profiles
   - Creates a matrix with hour and storage value information
4. Returns a dictionary where keys are day numbers and values are the corresponding storage profile matrices

The storage profiles are expected to have a structure where the first column contains date information
and subsequent columns (2-25) contain hourly storage data.
"""
function create_time_series_storage_profile(storage_profiles, num_days=nothing)
    # Get number of days
    total_days = size(storage_profiles, 1)
    hours_per_day = 24
    
    # If days are specified, ensure not exceeding available days
    if num_days === nothing
        num_days = total_days
    else
        num_days = min(num_days, total_days)
    end
    
    @info "Processing storage profile data, total $num_days days"
    
    # Create a dictionary to store storage profile time series data by day
    # Format: day_storage[day] = [hour, storage_value]
    day_storage = Dict{Int, Array{Float64, 2}}()
    
    # Process by day
    for day in 1:num_days
        @info "Processing day $day storage profile data, total $num_days days"
        
        # Create storage profile matrix for the day [hour, storage_value]
        day_storage_matrix = zeros(hours_per_day, 2)
        
        # Process each hour of the day
        for hour in 1:hours_per_day
            # Column index (first column is date, so start from second column)
            col_idx = hour + 1
            
            # Get storage value
            storage_value = storage_profiles[day, col_idx]
            
            # Fill storage profile matrix
            day_storage_matrix[hour, :] = [hour, storage_value]
        end
        
        # Store the day's storage profile data in the dictionary
        day_storage[day] = day_storage_matrix
    end
    
    return day_storage
end



"""
    runtdpf(case, data, load_names, price_profiles, irradiance_profiles, opt)

Run time-domain power flow analysis for multiple days with varying loads, prices, and solar irradiance.

# Arguments
- `case`: Power system case data structure
- `data`: DataFrame containing load time series data
- `load_names`: Vector of load names
- `price_profiles`: Matrix containing electricity price data
- `irradiance_profiles`: Matrix containing solar irradiance data
- `opt`: Options for power flow calculation

# Returns
- 3D array of results where dimensions represent [island, day, hour]

# Description
This function performs time-domain power flow analysis across multiple days, considering varying loads,
electricity prices, and solar irradiance. It processes the data day by day and handles multiple grid islands.

The function performs the following steps:
1. Determines the number of days from the time series data
2. Validates that the number of time points is divisible by 24 (hours per day)
3. Generates daily load, price, and irradiance data using helper functions
4. Converts the power system case to JPC format
5. Extracts grid islands from the power system
6. For each day and each island (in parallel):
   - Extracts the relevant load data for the island
   - Performs power flow calculations for each hour using the `run_single_day` function
7. Returns a 3D array of results organized by island, day, and hour

The parallel processing is implemented using Julia's threading capabilities (@threads),
with the outer loop over days being parallelized for efficiency.
"""
function runtdpf(case, data, load_names, price_profiles, irradiance_profiles, opt; storage_profiles=nothing)
    # Get number of time points and days
    num_timepoints = size(data, 1)
    hours_per_day = 24
    num_days = div(num_timepoints, hours_per_day)
    
    # Check if time points can be divided by 24
    if num_timepoints % hours_per_day != 0
        @warn "Number of time points ($num_timepoints) is not a multiple of 24, will only process complete days"
        num_timepoints = num_days * hours_per_day  # Adjust to complete days
    end
    
    # Generate load data organized by day (with loads on the same bus merged)
    day_loads = TimeDomainPowerFlow.create_time_series_loads(case, data, load_names, num_days)
    day_price = TimeDomainPowerFlow.create_time_series_prices(price_profiles, num_days)
    day_irradiance = TimeDomainPowerFlow.create_time_series_irradiance(irradiance_profiles, num_days)

    if storage_profiles !== nothing
        day_storage = TimeDomainPowerFlow.create_time_series_storage_profile(storage_profiles, num_days)
    end
    # Convert powerflow case to JPC format
    jpc = TimeDomainPowerFlow.JuliaPowerCase2Jpc(case)

    # Extract grid islands
    jpc_list, isolated = TimeDomainPowerFlow.extract_islands_acdc(jpc)
    n_islands = length(jpc_list)

    # Create results array
    results = Array{Any}(undef, n_islands, num_days, 24)
    
    # Process each day in parallel
    @threads for d in 1:num_days
        @info "Processing day $d, total $num_days days"
        
        # Get load data for the day
        if !haskey(day_loads, d)
            @error "No load data for day $d"
            continue
        end
        
        day_load_matrix = day_loads[d]
        day_price_line = day_price[d]
        day_irradiance_line = day_irradiance[d]
        if storage_profiles !== nothing
            day_storage_line = day_storage[d]
        else
            day_storage_line = nothing
        end
        # Extract load matrix for the day
        load_matrix_list, isolated_load_matrix = TimeDomainPowerFlow.extract_load_matrix_by_islands(day_load_matrix, jpc_list)
        
        # Perform power flow calculation for each island
        # Note: Not using @threads here because outer loop is already parallel
        for i in 1:n_islands
            results[i,d,:] = run_single_day(jpc_list[i], opt, load_matrix_list[i], day_price_line, day_irradiance_line, day_storage_line)
        end

        @info "Completed processing day $d"
    end
    
    return results
end

"""
    extract_load_matrix_by_islands(day_load_matrix, jpc_list)

Extract and organize load data by power system islands.

# Arguments
- `day_load_matrix`: Matrix containing load data for a day, where rows represent loads and columns include time, bus ID, and load values
- `jpc_list`: List of power system cases representing different islands in the network

# Returns
- `load_matrix_list`: List of load matrices, one for each island
- `isolated_load_matrix`: Matrix containing load data for buses not belonging to any island

# Description
This function distributes load data among different power system islands based on bus IDs.
It separates the daily load data into distinct matrices for each island in the power system.

The function performs the following steps:
1. Extracts bus IDs from the load data matrix (assumed to be in the second column)
2. For each island in the power system:
   - Identifies the AC buses belonging to the island
   - Finds load data rows associated with those buses
   - Creates a load matrix specific to that island
   - If no loads are found for an island, creates an empty matrix with the same column structure
3. Identifies loads on buses that don't belong to any valid island (isolated buses)
4. Returns both the list of island-specific load matrices and the matrix for isolated loads

This function is essential for distributed power flow analysis where each island 
needs to be processed separately with its corresponding load data.
"""
function extract_load_matrix_by_islands(day_load_matrix, jpc_list)
    # Initialize result list to store load matrix for each island
    load_matrix_list = []
    
    # Assume first column of day_load_matrix is time, second column is bus number
    # Extract bus IDs from day_load_matrix
    bus_ids_in_load = day_load_matrix[:, 2]
    
    # Process each island
    for i in eachindex(jpc_list)
        jpck = jpc_list[i]
        
        # Get AC bus IDs in the current island
        island_ac_buses = jpck.busAC[:, BUS_I]
        
        # Find rows in day_load_matrix belonging to the current island
        load_indices = findall(bus_id -> bus_id in island_ac_buses, bus_ids_in_load)
        
        if !isempty(load_indices)
            # Extract corresponding load matrix
            island_load_matrix = day_load_matrix[load_indices, :]
            push!(load_matrix_list, island_load_matrix)
        else
            # If current island has no corresponding loads, add an empty matrix
            push!(load_matrix_list, similar(day_load_matrix, 0, size(day_load_matrix, 2)))
        end
    end
    
    # Process isolated buses (buses not in any valid island)
    isolated_load_indices = findall(bus_id -> !any(bus_id in jpck.busAC[:, BUS_I] for jpck in jpc_list), bus_ids_in_load)
    isolated_load_matrix = isempty(isolated_load_indices) ? similar(day_load_matrix, 0, size(day_load_matrix, 2)) : day_load_matrix[isolated_load_indices, :]
    
    return load_matrix_list, isolated_load_matrix
end
