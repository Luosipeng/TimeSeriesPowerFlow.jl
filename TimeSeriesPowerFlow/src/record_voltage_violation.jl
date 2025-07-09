"""
    record_voltage_violation(results, bus_name, case, time_day, bus_type = "AC"; save_path = nothing, save_format = "pdf")

Analyze and visualize voltage violations for a specified bus in a power system.

# Arguments
- `results`: Simulation results containing bus voltage data
- `bus_name`: Name of the bus to analyze voltage violations for
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `bus_type`: Type of bus to analyze:
  - "AC": AC bus
  - "DC": DC bus
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot (default: "pdf")

# Returns
- `p`: The generated plot showing voltage violations
- `stats`: A named tuple containing detailed violation statistics:
  - `bus_name`: Name of the analyzed bus
  - `total_violations`: Total number of voltage violations
  - `under_limit_count`: Number of violations below the lower limit
  - `over_limit_count`: Number of violations above the upper limit
  - `worst_under_limit`: Lowest voltage value during violations
  - `worst_over_limit`: Highest voltage value during violations
  - `under_limit_indices`: Time indices of under-voltage violations
  - `over_limit_indices`: Time indices of over-voltage violations
  - `max_under_continuous`: Maximum duration of continuous under-voltage violations
  - `max_over_continuous`: Maximum duration of continuous over-voltage violations
  - `max_any_continuous`: Maximum duration of any continuous violations
  - `under_continuous_periods`: Details of continuous under-voltage violation periods
  - `over_continuous_periods`: Details of continuous over-voltage violation periods
  - `any_continuous_periods`: Details of all continuous violation periods

# Description
This function analyzes voltage violations for a specified bus by comparing voltage magnitude values against
defined upper and lower limits. It identifies individual violations, continuous violation periods, and
calculates statistics such as violation counts, maximum deviations, and longest violation durations.

The function creates a comprehensive visualization showing:
1. Voltage magnitude time series
2. Upper and lower voltage limits
3. Highlighted violation points
4. Marked most severe violations
5. Highlighted longest continuous violation periods
6. Detailed violation statistics

Additionally, the function prints detailed tables of individual violations and continuous violation periods,
including time information, violation types, voltage values, and deviation percentages.
"""

function record_voltage_violation(results, bus_name, case, time_day, bus_type = "AC"; save_path = nothing, save_format = "pdf")
    default(fontfamily="Microsoft YaHei")
    # Create a time series voltage magnitude matrix, first column represents the node number
    num_bus_AC = length(case.busesAC)
    voltage_magnitude_time_series_AC = zeros(num_bus_AC, time_day*24 + 1)  # First column for node number, subsequent columns for voltage magnitude at each time point
    voltage_magnitude_time_series_AC[:, 1] = 1:num_bus_AC  # First column for node number
    voltage_max_limit_AC = zeros(num_bus_AC, 2)  # Voltage maximum limit
    voltage_min_limit_AC = zeros(num_bus_AC, 2)
    voltage_max_limit_AC[:, 1] = 1:num_bus_AC  # First column for node number
    voltage_min_limit_AC[:, 1] = 1:num_bus_AC  # First column for node number

    for i in 1:size(results,1)
        for d in 1:time_day
            for hour in 1:24
                # Get voltage magnitude at current time point
                voltage_magnitude_time_series_AC[Int.(results[i, d, hour].busAC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busAC[:, VM]
            end
        end
    end

    for i in 1:size(results,1)
        # Get voltage maximum and minimum limits at current time point
        voltage_max_limit_AC[Int.(results[i, 1, 1].busAC[:, BUS_I]),2] .= results[i, 1, 1].busAC[:, VMAX]
        voltage_min_limit_AC[Int.(results[i, 1, 1].busAC[:, BUS_I]),2] .= results[i, 1, 1].busAC[:, VMIN]
    end

    num_bus_DC = length(case.busesDC)
    voltage_magnitude_time_series_DC = zeros(num_bus_DC, time_day*24 + 1)  # First column for node number, subsequent columns for voltage magnitude at each time point
    voltage_magnitude_time_series_DC[:, 1] = 1:num_bus_DC  # First column for node number
    voltage_max_limit_DC = zeros(num_bus_DC, 2)  # Voltage maximum limit
    voltage_min_limit_DC = zeros(num_bus_DC, 2)
    voltage_max_limit_DC[:, 1] = 1:num_bus_DC  # First column for node number
    voltage_min_limit_DC[:, 1] = 1:num_bus_DC  # First column for node number
    
    for i in 1:size(results,1)
        for d in 1:time_day
            for hour in 1:24
                # Get voltage magnitude at current time point
                voltage_magnitude_time_series_DC[Int.(results[i, d, hour].busDC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busDC[:, VM]
            end
        end
    end

    for i in 1:size(results,1)
        # Get voltage maximum and minimum limits at current time point
        voltage_max_limit_DC[Int.(results[i, 1, 1].busDC[:, BUS_I]),2] .= results[i, 1, 1].busDC[:, VMAX]
        voltage_min_limit_DC[Int.(results[i, 1, 1].busDC[:, BUS_I]),2] .= results[i, 1, 1].busDC[:, VMIN]
    end

    # Improved time axis label handling
    time_labels = String[]
    
    # Choose appropriate label format based on total days
    if time_day <= 3  # 3 days or less
        # One label every 4 hours
        for d in 1:time_day
            for h in 1:24
                if h % 4 == 0 || h == 1
                    push!(time_labels, "D$(d)-H$(h)")
                else
                    push!(time_labels, "")
                end
            end
        end
    elseif time_day <= 7  # One week or less
        # One label every 6 hours
        for d in 1:time_day
            for h in 1:24
                if h % 6 == 0 || h == 1
                    push!(time_labels, "D$(d)-H$(h)")
                else
                    push!(time_labels, "")
                end
            end
        end
    elseif time_day <= 31  # One month or less
        # Show only one label per day (first hour)
        for d in 1:time_day
            for h in 1:24
                if h == 1
                    push!(time_labels, "D$(d)")
                else
                    push!(time_labels, "")
                end
            end
        end
    else  # More than one month
        # Show labels every few days
        label_interval = max(1, round(Int, time_day / 20))  # Dynamically calculate label interval to ensure no more than about 20 labels
        for d in 1:time_day
            for h in 1:24
                if h == 1 && (d % label_interval == 0 || d == 1 || d == time_day)
                    push!(time_labels, "D$(d)")
                else
                    push!(time_labels, "")
                end
            end
        end
    end
    
    # Create time point index array to specify which points show labels
    time_points = 1:time_day*24
    xtick_indices = findall(x -> x != "", time_labels)
    xtick_labels = time_labels[xtick_indices]
    
    # Select appropriate data based on bus type
    if bus_type == "AC"
        bus_name_to_index = case.bus_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        voltage_magnitude_series = voltage_magnitude_time_series_AC[voltage_row, 2:end]
        voltage_min_limit = voltage_min_limit_AC[voltage_row, 2]
        voltage_max_limit = voltage_max_limit_AC[voltage_row, 2]
        bus_title = "AC Bus"
    else  # DC bus
        bus_name_to_index = case.busdc_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        voltage_magnitude_series = voltage_magnitude_time_series_DC[voltage_row, 2:end]
        voltage_min_limit = voltage_min_limit_DC[voltage_row, 2]
        voltage_max_limit = voltage_max_limit_DC[voltage_row, 2]
        bus_title = "DC Bus"
    end

    # Find time points with violations
    under_limit_indices = findall(v -> v < voltage_min_limit, voltage_magnitude_series)
    over_limit_indices = findall(v -> v > voltage_max_limit, voltage_magnitude_series)
    
    # Calculate violation counts
    under_limit_count = length(under_limit_indices)
    over_limit_count = length(over_limit_indices)
    total_violations = under_limit_count + over_limit_count
    
    # Find the most severe violations
    worst_under_limit = under_limit_count > 0 ? minimum(voltage_magnitude_series[under_limit_indices]) : NaN
    worst_over_limit = over_limit_count > 0 ? maximum(voltage_magnitude_series[over_limit_indices]) : NaN
    
    # Find the time points of the most severe violations
    worst_under_limit_idx = under_limit_count > 0 ? under_limit_indices[argmin(voltage_magnitude_series[under_limit_indices])] : 0
    worst_over_limit_idx = over_limit_count > 0 ? over_limit_indices[argmax(voltage_magnitude_series[over_limit_indices])] : 0
    
    # Calculate continuous violation duration
    # Create violation status arrays
    under_limit_status = zeros(Bool, length(voltage_magnitude_series))
    over_limit_status = zeros(Bool, length(voltage_magnitude_series))
    
    under_limit_status[under_limit_indices] .= true
    over_limit_status[over_limit_indices] .= true
    
    # Calculate continuous violation duration
    # Statistics for continuous violations below the lower limit
    current_under_continuous = 0
    max_under_continuous = 0
    under_continuous_periods = []  # Store all continuous violation intervals [(start, end, length), ...]
    
    for i in 1:length(under_limit_status)
        if under_limit_status[i]
            current_under_continuous += 1
        else
            if current_under_continuous > 0
                # Record continuous violation interval
                start_idx = i - current_under_continuous
                end_idx = i - 1
                push!(under_continuous_periods, (start_idx, end_idx, current_under_continuous))
                
                # Update maximum continuous violation duration
                max_under_continuous = max(max_under_continuous, current_under_continuous)
                current_under_continuous = 0
            end
        end
    end
    
    # Handle continuous violations at the end of the sequence
    if current_under_continuous > 0
        start_idx = length(under_limit_status) - current_under_continuous + 1
        end_idx = length(under_limit_status)
        push!(under_continuous_periods, (start_idx, end_idx, current_under_continuous))
        max_under_continuous = max(max_under_continuous, current_under_continuous)
    end
    
    # Statistics for continuous violations above the upper limit
    current_over_continuous = 0
    max_over_continuous = 0
    over_continuous_periods = []  # Store all continuous violation intervals
    
    for i in 1:length(over_limit_status)
        if over_limit_status[i]
            current_over_continuous += 1
        else
            if current_over_continuous > 0
                # Record continuous violation interval
                start_idx = i - current_over_continuous
                end_idx = i - 1
                push!(over_continuous_periods, (start_idx, end_idx, current_over_continuous))
                
                # Update maximum continuous violation duration
                max_over_continuous = max(max_over_continuous, current_over_continuous)
                current_over_continuous = 0
            end
        end
    end
    
    # Handle continuous violations at the end of the sequence
    if current_over_continuous > 0
        start_idx = length(over_limit_status) - current_over_continuous + 1
        end_idx = length(over_limit_status)
        push!(over_continuous_periods, (start_idx, end_idx, current_over_continuous))
        max_over_continuous = max(max_over_continuous, current_over_continuous)
    end
    
    # Calculate total continuous violation duration (below lower limit or above upper limit)
    any_limit_status = under_limit_status .| over_limit_status
    current_any_continuous = 0
    max_any_continuous = 0
    any_continuous_periods = []
    
    for i in 1:length(any_limit_status)
        if any_limit_status[i]
            current_any_continuous += 1
        else
            if current_any_continuous > 0
                # Record continuous violation interval
                start_idx = i - current_any_continuous
                end_idx = i - 1
                push!(any_continuous_periods, (start_idx, end_idx, current_any_continuous))
                
                # Update maximum continuous violation duration
                max_any_continuous = max(max_any_continuous, current_any_continuous)
                current_any_continuous = 0
            end
        end
    end
    
    # Handle continuous violations at the end of the sequence
    if current_any_continuous > 0
        start_idx = length(any_limit_status) - current_any_continuous + 1
        end_idx = length(any_limit_status)
        push!(any_continuous_periods, (start_idx, end_idx, current_any_continuous))
        max_any_continuous = max(max_any_continuous, current_any_continuous)
    end
    
    # Find the longest continuous violation interval
    longest_under_period = isempty(under_continuous_periods) ? nothing : under_continuous_periods[argmax([p[3] for p in under_continuous_periods])]
    longest_over_period = isempty(over_continuous_periods) ? nothing : over_continuous_periods[argmax([p[3] for p in over_continuous_periods])]
    longest_any_period = isempty(any_continuous_periods) ? nothing : any_continuous_periods[argmax([p[3] for p in any_continuous_periods])]
    
    # Plot voltage magnitude time series
    p = plot(time_points, voltage_magnitude_series, 
             label="Voltage Magnitude", 
             linewidth=2, 
             color=:blue,
             title="$(bus_title) \"$(bus_name)\" Voltage Violation Analysis",
             ylabel="Voltage Magnitude (p.u.)",
             grid=true,
             legend=:topright,
             size=(900, 600))  # Increase chart width to provide more space for labels
    
    # Add upper and lower limit lines
    hline!([voltage_min_limit], label="Lower Limit ($(round(voltage_min_limit, digits=3)))", color=:red, linestyle=:dash, linewidth=1.5)
    hline!([voltage_max_limit], label="Upper Limit ($(round(voltage_max_limit, digits=3)))", color=:red, linestyle=:dash, linewidth=1.5)
    
    # Mark violation points
    if under_limit_count > 0
        scatter!(under_limit_indices, voltage_magnitude_series[under_limit_indices], 
                 label="Below Lower Limit ($(under_limit_count) times)", 
                 color=:orange, 
                 markersize=4)
    end
    
    if over_limit_count > 0
        scatter!(over_limit_indices, voltage_magnitude_series[over_limit_indices], 
                 label="Above Upper Limit ($(over_limit_count) times)", 
                 color=:purple, 
                 markersize=4)
    end
    
    # Mark the most severe violation points
    if under_limit_count > 0
        scatter!([worst_under_limit_idx], [worst_under_limit], 
                 label="Lowest Point: $(round(worst_under_limit, digits=4))", 
                 color=:red, 
                 markershape=:star5, 
                 markersize=8)
    end
    
    if over_limit_count > 0
        scatter!([worst_over_limit_idx], [worst_over_limit], 
                 label="Highest Point: $(round(worst_over_limit, digits=4))", 
                 color=:darkred, 
                 markershape=:star5, 
                 markersize=8)
    end
    
    # Highlight the longest continuous violation interval
    if !isnothing(longest_under_period)
        start_idx, end_idx, length_period = longest_under_period
        highlight_x = start_idx:end_idx
        highlight_y = voltage_magnitude_series[highlight_x]
        plot!(highlight_x, highlight_y, 
              linewidth=4, 
              color=:orange, 
              alpha=0.5, 
              label="Longest Below Lower Limit Interval ($(length_period) hours)")
    end
    
    if !isnothing(longest_over_period)
        start_idx, end_idx, length_period = longest_over_period
        highlight_x = start_idx:end_idx
        highlight_y = voltage_magnitude_series[highlight_x]
        plot!(highlight_x, highlight_y, 
              linewidth=4, 
              color=:purple, 
              alpha=0.5, 
              label="Longest Above Upper Limit Interval ($(length_period) hours)")
    end
    
    # Set optimized x-axis tick labels
    p = plot!(p, xticks=(xtick_indices, xtick_labels), xrotation=45)
    
    # Add auxiliary grid lines to make time points easier to match
    p = plot!(p, minorgrid=true, minorgridalpha=0.1)
    
    # Add date range annotation
    if time_day > 7
        day_range_text = "Time span: $(time_day) days"
        annotate!(0.5*length(time_points), voltage_min_limit - 0.02, 
                 text(day_range_text, :center, 8))
    end
    
    # Add annotation with violation statistics
    violation_text = "Total Violations: $(total_violations)\n" *
                     "Below Lower Limit: $(under_limit_count)\n" *
                     "Above Upper Limit: $(over_limit_count)\n" *
                     "Max Continuous Violation: $(max_any_continuous) hours\n" *
                     "Max Continuous Below Limit: $(max_under_continuous) hours\n" *
                     "Max Continuous Above Limit: $(max_over_continuous) hours\n" *
                     "Lowest Voltage: $(round(isnan(worst_under_limit) ? minimum(voltage_magnitude_series) : worst_under_limit, digits=4))\n" *
                     "Highest Voltage: $(round(isnan(worst_over_limit) ? maximum(voltage_magnitude_series) : worst_over_limit, digits=4))"
    
    annotate!(0.8 * length(voltage_magnitude_series), 
              voltage_min_limit + 0.6 * (voltage_max_limit - voltage_min_limit), 
              text(violation_text, :left, 8))
    
    # Save the plot if a path is provided
    if save_path !== nothing
        # If save_path doesn't include file extension, add it based on save_format
        if !contains(save_path, ".")
            save_path = "$(save_path).$(save_format)"
        end
        savefig(p, save_path)
        println("Plot saved to: $(save_path)")
    end
    
    # Create violation details table
    if total_violations > 0
        # Prepare violation details data
        violation_details = []
        
        # Add details for violations below the lower limit
        for i in under_limit_indices
            # Get actual time label, even if it's an empty string in time_labels
            actual_time_label = "D$(ceil(Int, i/24))-H$((i-1)%24+1)"
            
            push!(violation_details, (
                Time = actual_time_label,
                Type = "Below Lower Limit",
                Voltage = round(voltage_magnitude_series[i], digits=4),
                Deviation = round(voltage_magnitude_series[i] - voltage_min_limit, digits=4),
                DeviationPercent = round((voltage_magnitude_series[i] - voltage_min_limit) / voltage_min_limit * 100, digits=2)
            ))
        end
        
        # Add details for violations above the upper limit
        for i in over_limit_indices
            # Get actual time label, even if it's an empty string in time_labels
            actual_time_label = "D$(ceil(Int, i/24))-H$((i-1)%24+1)"
            
            push!(violation_details, (
                Time = actual_time_label,
                Type = "Above Upper Limit",
                Voltage = round(voltage_magnitude_series[i], digits=4),
                Deviation = round(voltage_magnitude_series[i] - voltage_max_limit, digits=4),
                DeviationPercent = round((voltage_magnitude_series[i] - voltage_max_limit) / voltage_max_limit * 100, digits=2)
            ))
        end
        
        # Sort by time
        sort!(violation_details, by = x -> x.Time)
        
        # Create violation details table
        violation_table = DataFrame(
            Time = [v.Time for v in violation_details],
            Type = [v.Type for v in violation_details],
            Voltage = [v.Voltage for v in violation_details],
            Deviation = [v.Deviation for v in violation_details],
            DeviationPercent = [v.DeviationPercent for v in violation_details]
        )
        
        # Output violation details table
        println("\nNode $(bus_name) Voltage Violation Details (Sorted by Time):")
        println(violation_table)
    else
        println("\nNode $(bus_name) has no voltage violations.")
    end
    
    # Create continuous violation interval table
    if !isempty(any_continuous_periods)
        # Prepare continuous violation interval data
        continuous_periods_data = []
        
        for (i, period) in enumerate(any_continuous_periods)
            start_idx, end_idx, length_period = period
            
            # Get actual start and end time labels
            start_time_label = "D$(ceil(Int, start_idx/24))-H$((start_idx-1)%24+1)"
            end_time_label = "D$(ceil(Int, end_idx/24))-H$((end_idx-1)%24+1)"
            
            # Determine violation type
            if all(under_limit_status[start_idx:end_idx])
                period_type = "Below Lower Limit"
            elseif all(over_limit_status[start_idx:end_idx])
                period_type = "Above Upper Limit"
            else
                period_type = "Mixed Violations"
            end
            
            # Calculate maximum deviation in the interval
            if period_type == "Below Lower Limit" || period_type == "Mixed Violations"
                under_values = voltage_magnitude_series[start_idx:end_idx][under_limit_status[start_idx:end_idx]]
                if !isempty(under_values)
                    worst_value = minimum(under_values)
                    worst_deviation = worst_value - voltage_min_limit
                    worst_deviation_pct = worst_deviation / voltage_min_limit * 100
                else
                    worst_value = NaN
                    worst_deviation = NaN
                    worst_deviation_pct = NaN
                end
            elseif period_type == "Above Upper Limit"
                over_values = voltage_magnitude_series[start_idx:end_idx][over_limit_status[start_idx:end_idx]]
                worst_value = maximum(over_values)
                worst_deviation = worst_value - voltage_max_limit
                worst_deviation_pct = worst_deviation / voltage_max_limit * 100
            end
            
            push!(continuous_periods_data, (
                StartTime = start_time_label,
                EndTime = end_time_label,
                Duration = length_period,
                ViolationType = period_type,
                WorstVoltage = round(worst_value, digits=4),
                MaxDeviation = round(worst_deviation, digits=4),
                MaxDeviationPercent = round(worst_deviation_pct, digits=2)
            ))
        end
        
        # Sort by duration
        sort!(continuous_periods_data, by = x -> x.Duration, rev = true)
        
        # Create continuous violation interval table
        periods_table = DataFrame(
            StartTime = [p.StartTime for p in continuous_periods_data],
            EndTime = [p.EndTime for p in continuous_periods_data],
            Duration = [p.Duration for p in continuous_periods_data],
            ViolationType = [p.ViolationType for p in continuous_periods_data],
            WorstVoltage = [p.WorstVoltage for p in continuous_periods_data],
            MaxDeviation = [p.MaxDeviation for p in continuous_periods_data],
            MaxDeviationPercent = [p.MaxDeviationPercent for p in continuous_periods_data]
        )
        
        # Output continuous violation interval table
        println("\nNode $(bus_name) Continuous Voltage Violation Intervals (Sorted by Duration):")
        println(periods_table)
    end
    
    # Return plot result and statistics
    stats = (
        bus_name = bus_name,
        total_violations = total_violations,
        under_limit_count = under_limit_count,
        over_limit_count = over_limit_count,
        worst_under_limit = worst_under_limit,
        worst_over_limit = worst_over_limit,
        under_limit_indices = under_limit_indices,
        over_limit_indices = over_limit_indices,
        max_under_continuous = max_under_continuous,
        max_over_continuous = max_over_continuous,
        max_any_continuous = max_any_continuous,
        under_continuous_periods = under_continuous_periods,
        over_continuous_periods = over_continuous_periods,
        any_continuous_periods = any_continuous_periods
    )
    
    return p, stats
end
