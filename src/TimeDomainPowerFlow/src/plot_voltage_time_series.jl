"""
    plot_voltage_time_series(results, bus_name, case, time_day, bus_type = "AC"; save_path = nothing, save_format = "pdf")

Plot voltage time series for a specified bus in a power system.

# Arguments
- `results`: Simulation results containing bus voltage data
- `bus_name`: Name of the bus to plot voltage for
- `case`: Power system case data
- `time_day`: Number of days in the simulation
- `bus_type`: Type of bus to plot:
  - "AC": AC bus (plots both voltage magnitude and angle)
  - "DC": DC bus (plots only voltage magnitude)
- `save_path`: Optional path to save the plot
- `save_format`: Format to save the plot (default: "pdf")

# Returns
- `plot_result`: The generated plot showing voltage time series
- `voltage_magnitude`: Vector of voltage magnitude values for the specified bus
- `voltage_angle`: Vector of voltage angle values for the specified bus (only for AC buses, empty for DC buses)

# Description
This function extracts and visualizes voltage data for a specified bus over time. For AC buses,
it creates two subplots showing both voltage magnitude (in per unit) and voltage angle (in degrees).
For DC buses, it creates a single plot showing only voltage magnitude.

The function automatically adjusts the time axis labels based on the simulation duration to ensure
readability. For short time periods, it shows detailed hour labels, while for longer periods,
it shows day labels at appropriate intervals.

The plot includes grid lines and time span information for easy reference, and can be saved
to a file if a save path is provided.

In addition to the plot, the function also returns the actual voltage magnitude and angle values
for further analysis or processing.
"""

function plot_voltage_time_series(results, bus_name, case, time_day, bus_type = "AC"; save_path = nothing, save_format = "pdf")
    default(fontfamily="Microsoft YaHei")
    # Create a time series voltage matrix, the first column represents the node number
    num_bus_AC = length(case.busesAC)
    voltage_magnitude_time_series_AC = zeros(num_bus_AC, time_day*24 + 1)  # First column for node number, subsequent columns for voltage magnitude at each time point
    voltage_angle_time_series_AC = zeros(num_bus_AC, time_day*24 + 1)  # First column for node number, subsequent columns for voltage angle at each time point
    voltage_magnitude_time_series_AC[:, 1] = 1:num_bus_AC  # First column for node number
    voltage_angle_time_series_AC[:, 1] = 1:num_bus_AC  # First column for node number
    for i in eachindex(results[:,1,1])
        for d in 1:time_day
            for hour in 1:24
                # Get voltage magnitude at current time point
                voltage_magnitude_time_series_AC[Int.(results[i, d, hour].busAC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busAC[:, VM]  # Assuming VM is the column for voltage magnitude
                # Get voltage angle at current time point
                voltage_angle_time_series_AC[Int.(results[i, d, hour].busAC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busAC[:, VA]  # Assuming VA is the column for voltage angle
            end
        end
    end

    num_bus_DC = length(case.busesDC)
    voltage_magnitude_time_series_DC = zeros(num_bus_DC, time_day*24 + 1)  # First column for node number, subsequent columns for voltage magnitude at each time point
    voltage_magnitude_time_series_DC[:, 1] = 1:num_bus_DC  # First column for node number
    for i in eachindex(results[:,1,1])
        for d in 1:time_day
            for hour in 1:24
                # Get voltage magnitude at current time point
                voltage_magnitude_time_series_DC[Int.(results[i, d, hour].busDC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busDC[:, VM]  # Assuming VM is the column for voltage magnitude
            end
        end
    end

    # Create array of time point indices
    time_points = 1:time_day*24
    
    # Dynamically adjust label density based on time span
    # Calculate appropriate label interval to avoid label overlap
    if time_day <= 3  # 3 days or less
        # One label every 8 hours
        label_interval_hours = 8
        label_format = "D%d-H%d"  # Show day and hour
    elseif time_day <= 7  # One week or less
        # One label every 12 hours
        label_interval_hours = 12
        label_format = "D%d-H%d"  # Show day and hour
    elseif time_day <= 31  # One month or less
        # One label per day
        label_interval_hours = 24
        label_format = "D%d"  # Show only day
    else  # More than one month
        # One label every few days
        days_interval = max(1, round(Int, time_day / 15))  # Dynamically calculate day interval to ensure about 15 labels total
        label_interval_hours = days_interval * 24
        label_format = "D%d"  # Show only day
    end
    
    # Create labels and corresponding indices
    xtick_indices = []
    xtick_labels = []
    
    for t in 1:length(time_points)
        day = ceil(Int, t/24)
        hour = (t-1)%24 + 1
        
        # Decide whether to add a label based on the interval
        if (t-1) % label_interval_hours == 0 || t == 1 || t == length(time_points)
            push!(xtick_indices, t)
            
            if label_format == "D%d-H%d"
                push!(xtick_labels, @sprintf("D%d-H%d", day, hour))
            else
                push!(xtick_labels, @sprintf("D%d", day))
            end
        end
    end
    
    # Set plot size, adjusting based on time span
    plot_width = min(1000, max(800, time_day * 20))  # Dynamically adjust width based on number of days
    plot_size = (plot_width, 600)
    
    # Initialize return values
    voltage_magnitude = Float64[]
    voltage_angle = Float64[]
    
    if bus_type == "AC"
        bus_name_to_index = case.bus_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        voltage_magnitude_series = voltage_magnitude_time_series_AC[voltage_row, 2:end]  # Get voltage magnitude for the corresponding bus
        voltage_angle_series = voltage_angle_time_series_AC[voltage_row, 2:end]
        
        # Store values for return
        voltage_magnitude = copy(voltage_magnitude_series)
        voltage_angle = copy(voltage_angle_series)
        
        # Create two subplots
        p1 = plot(time_points, voltage_magnitude_series, 
                 label="", 
                 linewidth=2, 
                 color=:blue,
                 title="AC Bus \"$(bus_name)\" Voltage Magnitude Time Series",
                 ylabel="Voltage Magnitude (p.u.)",
                 grid=true)
                 
        p2 = plot(time_points, voltage_angle_series, 
                 label="", 
                 linewidth=2, 
                 color=:red,
                 title="AC Bus \"$(bus_name)\" Voltage Angle Time Series",
                 xlabel="Time",
                 ylabel="Voltage Angle (degree)",
                 grid=true)
        
        # Set x-axis tick labels
        plot!(p1, xticks=(xtick_indices, xtick_labels), xrotation=45)
        plot!(p2, xticks=(xtick_indices, xtick_labels), xrotation=45)
        
        # Add auxiliary grid lines to make time points easier to correspond
        plot!(p1, minorgrid=true, minorgridalpha=0.1)
        plot!(p2, minorgrid=true, minorgridalpha=0.1)
        
        # Combine the two subplots
        plot_result = plot(p1, p2, layout=(2,1), size=plot_size)
        
        # Add time span information
        if time_day > 7
            annotate!(subplot=1, 0.5*length(time_points), minimum(voltage_magnitude_series) + 0.9*(maximum(voltage_magnitude_series)-minimum(voltage_magnitude_series)), 
                     text("Time span: $(time_day) days", :right, 8))
        end
        
    else  # DC bus
        bus_name_to_index = case.busdc_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        voltage_magnitude_series = voltage_magnitude_time_series_DC[voltage_row, 2:end]  # Get voltage magnitude for the corresponding bus

        # Store values for return
        voltage_magnitude = copy(voltage_magnitude_series)
        # For DC buses, voltage angle is not applicable, so return empty array
        voltage_angle = Float64[]
        
        # Plot voltage magnitude time series
        plot_result = plot(time_points, voltage_magnitude_series, 
                          label="", 
                          linewidth=2, 
                          color=:green,
                          title="DC Bus \"$(bus_name)\" Voltage Magnitude Time Series",
                          xlabel="Time",
                          ylabel="Voltage Magnitude (p.u.)",
                          grid=true,
                          size=plot_size)
        
        # Set x-axis tick labels
        plot!(plot_result, xticks=(xtick_indices, xtick_labels), xrotation=45)
        
        # Add auxiliary grid lines
        plot!(plot_result, minorgrid=true, minorgridalpha=0.1)
        
        # Add time span information
        if time_day > 7
            annotate!(0.5*length(time_points), minimum(voltage_magnitude_series) + 0.9*(maximum(voltage_magnitude_series)-minimum(voltage_magnitude_series)), 
                     text("Time span: $(time_day) days", :right, 8))
        end
    end
    
    # Save the plot if a path is provided
    if save_path !== nothing
        # If save_path doesn't include file extension, add it based on save_format
        if !contains(save_path, ".")
            save_path = "$(save_path).$(save_format)"
        end
        savefig(plot_result, save_path)
        println("Plot saved to: $(save_path)")
    end
    
    return plot_result, voltage_magnitude, voltage_angle
end
