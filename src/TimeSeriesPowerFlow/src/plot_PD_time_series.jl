"""
Plot the time series of active load
Parameters:
results - Result dataset
bus_name - Bus name
case - System case
time_day - Number of days
bus_type - Bus type (default is "AC")
"""
function plot_PD_time_series(results, bus_name, case, time_day, bus_type = "AC")
    default(fontfamily="Microsoft YaHei")
    # Create a time series active load matrix, first column represents node number
    num_bus_AC = length(case.busesAC)
    PD_time_series_AC = zeros(num_bus_AC, time_day*24 + 1)  # First column for node number, subsequent columns for active load at each time point
    PD_time_series_AC[:, 1] = 1:num_bus_AC  # First column is node number
    
    for i in eachindex(results[:,1,1])
        for d in 1:time_day
            for hour in 1:24
                # Get active load at current time point
                PD_time_series_AC[Int.(results[i, d, hour].busAC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busAC[:, PD]  # Read PD data
            end
        end
    end

    # If DC bus active load is needed, process similarly
    num_bus_DC = length(case.busesDC)
    PD_time_series_DC = zeros(num_bus_DC, time_day*24 + 1)  # First column for node number, subsequent columns for active load at each time point
    PD_time_series_DC[:, 1] = 1:num_bus_DC  # First column is node number
    
    for i in eachindex(results[:,1,1])
        for d in 1:time_day
            for hour in 1:24
                # Get active load at current time point
                if size(results[i, d, hour].busDC, 2) >= PD  # Ensure DC bus has PD column
                    PD_time_series_DC[Int.(results[i, d, hour].busDC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busDC[:, PD]
                end
            end
        end
    end

    # Create time axis labels
    time_labels = String[]
    for d in 1:time_day
        for h in 1:24
            push!(time_labels, "D$(d)-H$(h)")
        end
    end
    
    if bus_type == "AC"
        bus_name_to_index = case.bus_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        PD_series = PD_time_series_AC[voltage_row, 2:end]  # Get active load for corresponding bus

        # Create time points array
        time_points = 1:length(PD_series)
        
        # Create chart
        plot_result = plot(time_points, PD_series, 
                         label="", 
                         linewidth=2, 
                         color=:blue,
                         title="AC Bus \"$(bus_name)\" Active Load Time Series",
                         xlabel="Time",
                         ylabel="Active Load (MW)",
                         grid=true)
        
        # Set x-axis tick labels
        if time_day * 24 <= 48  # If total time points don't exceed 48, show all
            xtick_indices = 1:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels))
        else  # Otherwise show one label every 6 hours
            xtick_indices = 1:6:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels[xtick_indices]))
        end
        
    else  # DC bus
        bus_name_to_index = case.busdc_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        PD_series = PD_time_series_DC[voltage_row, 2:end]  # Get active load for corresponding bus

        # Create time points array
        time_points = 1:length(PD_series)
        
        # Plot active load time series
        plot_result = plot(time_points, PD_series, 
                          label="", 
                          linewidth=2, 
                          color=:green,
                          title="DC Bus \"$(bus_name)\" Active Load Time Series",
                          xlabel="Time",
                          ylabel="Active Load (MW)",
                          grid=true)
        
        # Set x-axis tick labels
        if time_day * 24 <= 48
            xtick_indices = 1:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels))
        else
            xtick_indices = 1:6:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels[xtick_indices]))
        end
    end
    
    return plot_result
end
