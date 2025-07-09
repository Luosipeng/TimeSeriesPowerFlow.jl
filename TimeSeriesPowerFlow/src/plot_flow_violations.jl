"""
    plot_flow_violations(results, case, time_day, flow_limit = 3.0, plot_type = "summary", flow_direction = "max"; save_path = nothing, save_format = "pdf")

Plot power flow violations in power system branches.

# Arguments
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
- `save_format`: Format to save the plot (default: "pdf")

# Returns
- `plot_result`: The generated plot
- `violation_count`: Number of branches with violations at each time point
- `max_violation_percent`: Maximum violation percentage at each time point
- `total_violation_severity`: Sum of violation severity at each time point
- `violation_details`: Detailed information about each violation
- `branch_violation_stats`: Statistics for each branch with violations
"""
function plot_flow_violations(results, case, time_day, flow_limit = 3.0, plot_type = "summary", flow_direction = "max"; save_path = nothing, save_format = "pdf")
    # Ensure flow_limit is a numerical type
    flow_limit = float(flow_limit)  # Add this line to ensure flow_limit is a float
    
    default(fontfamily="Microsoft YaHei")
    # Get number of islands
    num_islands = size(results, 1)
    
    # Create arrays to store violation statistics
    violation_count = zeros(Int, time_day*24)  # Number of branches with violations at each time point
    max_violation_percent = zeros(time_day*24)  # Maximum violation percentage at each time point
    total_violation_severity = zeros(time_day*24)  # Sum of violation severity at each time point
    
    # Create dictionaries to store power flows for all branches
    branch_flows = Dict()
    branch_names = Dict()
    violation_status = Dict()  # Record whether each branch has a violation at each time point
    
    # Create reverse mapping from ID to name (since JuliaPowerCase only has bus_name_to_id)
    id_to_bus_name = Dict{Int, String}()
    for (name, id) in case.bus_name_to_id
        id_to_bus_name[id] = name
    end
    
    # Create data structure to store violation details
    violation_details = []
    
    # Create dictionary to track violation statistics for each branch
    branch_violation_stats = Dict()
    
    # Process each island separately
    for island in 1:num_islands
        # Get number of AC branches in this island
        num_branch_AC = size(results[island, 1, 1].branchAC, 1)
        
        # Get number of DC branches in this island
        num_branch_DC = size(results[island, 1, 1].branchDC, 1)

        # Initialize arrays to store power flows for each branch
        for b in 1:num_branch_AC
            branch_key = "AC_$(island)_$(b)"
            branch_flows[branch_key] = zeros(time_day*24)
            violation_status[branch_key] = zeros(Bool, time_day*24)
            
            # Get branch from and to buses
            from_bus = Int(results[island, 1, 1].branchAC[b, F_BUS])
            to_bus = Int(results[island, 1, 1].branchAC[b, T_BUS])
            
            # Try to get bus names
            from_name = haskey(id_to_bus_name, from_bus) ? id_to_bus_name[from_bus] : "Bus$from_bus"
            to_name = haskey(id_to_bus_name, to_bus) ? id_to_bus_name[to_bus] : "Bus$to_bus"
            
            branch_names[branch_key] = "$from_name - $to_name"
            
            # Initialize branch violation statistics
            branch_violation_stats[branch_key] = (
                total_violations = 0,           # Total number of violations
                max_continuous_violations = 0,  # Maximum continuous violation duration
                current_continuous = 0,         # Current continuous violation count
                max_flow = 0.0,                 # Maximum power flow
                max_flow_time = "",             # Time of maximum power flow
                max_violation_percent = 0.0,    # Maximum violation percentage
                max_violation_time = ""         # Time of maximum violation percentage
            )
        end
        
        for b in 1:num_branch_DC
            branch_key = "DC_$(island)_$(b)"
            branch_flows[branch_key] = zeros(time_day*24)
            violation_status[branch_key] = zeros(Bool, time_day*24)
            
            # Get DC branch from and to buses
            from_bus = Int(results[island, 1, 1].branchDC[b, F_BUS])
            to_bus = Int(results[island, 1, 1].branchDC[b, T_BUS])
            
            # Try to get bus names
            from_name = haskey(id_to_bus_name, from_bus) ? id_to_bus_name[from_bus] : "Bus$from_bus"
            to_name = haskey(id_to_bus_name, to_bus) ? id_to_bus_name[to_bus] : "Bus$to_bus"
            
            branch_names[branch_key] = "$from_name - $to_name (DC)"
            
            # Initialize branch violation statistics
            branch_violation_stats[branch_key] = (
                total_violations = 0,           # Total number of violations
                max_continuous_violations = 0,  # Maximum continuous violation duration
                current_continuous = 0,         # Current continuous violation count
                max_flow = 0.0,                 # Maximum power flow
                max_flow_time = "",             # Time of maximum power flow
                max_violation_percent = 0.0,    # Maximum violation percentage
                max_violation_time = ""         # Time of maximum violation percentage
            )
        end

        # Extract branch power flows from results and check for violations at each time point
        for d in 1:time_day
            for hour in 1:24
                time_idx = (d-1)*24 + hour
                time_label = "D$(d)-H$(hour)"
                
                # Process AC branches
                if num_branch_AC > 0
                    branch_data = results[island, d, hour].branchAC
                    
                    # Ensure branch data contains necessary columns
                    if size(branch_data, 2) >= PT
                        # Check power flow for each branch
                        for b in 1:num_branch_AC
                            branch_key = "AC_$(island)_$(b)"
                            
                            # Determine how to check power flow based on flow_direction parameter
                            if flow_direction == "max"
                                # Take the maximum absolute value of power flow
                                flow_value = max(abs(branch_data[b, PF]), abs(branch_data[b, PT]))
                            elseif flow_direction == "both"
                                # Check both forward and reverse power flows, take maximum
                                flow_value = max(abs(branch_data[b, PF]), abs(branch_data[b, PT]))
                            elseif flow_direction == "forward"
                                # Only check forward power flow (from from to to)
                                flow_value = abs(branch_data[b, PF])
                            elseif flow_direction == "reverse"
                                # Only check reverse power flow (from to to from)
                                flow_value = abs(branch_data[b, PT])
                            end
                            
                            branch_flows[branch_key][time_idx] = flow_value
                            
                            # Update maximum power flow record
                            if flow_value > branch_violation_stats[branch_key].max_flow
                                branch_violation_stats[branch_key] = merge(
                                    branch_violation_stats[branch_key],
                                    (max_flow = flow_value, max_flow_time = time_label)
                                )
                            end
                            
                            # Check for violation
                            if flow_value > flow_limit
                                violation_status[branch_key][time_idx] = true
                                violation_count[time_idx] += 1
                                
                                # Calculate violation percentage
                                violation_percent = (flow_value - flow_limit) / flow_limit * 100
                                
                                # Update maximum violation percentage
                                if violation_percent > max_violation_percent[time_idx]
                                    max_violation_percent[time_idx] = violation_percent
                                end
                                
                                # Add to total violation severity
                                total_violation_severity[time_idx] += violation_percent
                                
                                # Record violation details
                                push!(violation_details, (
                                    branch_name = branch_names[branch_key],
                                    branch_key = branch_key,
                                    time_label = time_label,
                                    time_idx = time_idx,
                                    flow_value = flow_value,
                                    violation_percent = violation_percent
                                ))
                                
                                # Update branch violation statistics
                                stats = branch_violation_stats[branch_key]
                                
                                # Increment total violations
                                total_violations = stats.total_violations + 1
                                
                                # Increment current continuous violations
                                current_continuous = stats.current_continuous + 1
                                
                                # Update maximum continuous violations
                                max_continuous = max(stats.max_continuous_violations, current_continuous)
                                
                                # Update maximum violation percentage
                                max_viol_percent = stats.max_violation_percent
                                max_viol_time = stats.max_violation_time
                                
                                if violation_percent > max_viol_percent
                                    max_viol_percent = violation_percent
                                    max_viol_time = time_label
                                end
                                
                                # Update statistics
                                branch_violation_stats[branch_key] = (
                                    total_violations = total_violations,
                                    max_continuous_violations = max_continuous,
                                    current_continuous = current_continuous,
                                    max_flow = stats.max_flow,
                                    max_flow_time = stats.max_flow_time,
                                    max_violation_percent = max_viol_percent,
                                    max_violation_time = max_viol_time
                                )
                            else
                                # If no violation at current time point, reset continuous violation count
                                if branch_violation_stats[branch_key].current_continuous > 0
                                    branch_violation_stats[branch_key] = merge(
                                        branch_violation_stats[branch_key],
                                        (current_continuous = 0,)
                                    )
                                end
                            end
                        end
                    end
                end
                
                # Process DC branches
                if num_branch_DC > 0 
                    dc_branch_data = results[island, d, hour].branchDC
                    
                    # Ensure DC branch data contains necessary columns
                    if size(dc_branch_data, 2) >= PT
                        # Check power flow for each DC branch
                        for b in 1:num_branch_DC
                            branch_key = "DC_$(island)_$(b)"
                            
                            # Determine how to check power flow based on flow_direction parameter
                            if flow_direction == "max" || flow_direction == "both"
                                # Take the maximum absolute value of power flow
                                flow_value = max(abs(dc_branch_data[b, PF]), abs(dc_branch_data[b, PT]))
                            elseif flow_direction == "forward"
                                # Only check forward power flow
                                flow_value = abs(dc_branch_data[b, PF])
                            elseif flow_direction == "reverse"
                                # Only check reverse power flow
                                flow_value = abs(dc_branch_data[b, PT])
                            end
                            
                            branch_flows[branch_key][time_idx] = flow_value
                            
                            # Update maximum power flow record
                            if flow_value > branch_violation_stats[branch_key].max_flow
                                branch_violation_stats[branch_key] = merge(
                                    branch_violation_stats[branch_key],
                                    (max_flow = flow_value, max_flow_time = time_label)
                                )
                            end
                            
                            # Check for violation
                            if flow_value > flow_limit
                                violation_status[branch_key][time_idx] = true
                                violation_count[time_idx] += 1
                                
                                # Calculate violation percentage
                                violation_percent = (flow_value - flow_limit) / flow_limit * 100
                                
                                # Update maximum violation percentage
                                if violation_percent > max_violation_percent[time_idx]
                                    max_violation_percent[time_idx] = violation_percent
                                end
                                
                                # Add to total violation severity
                                total_violation_severity[time_idx] += violation_percent
                                
                                # Record violation details
                                push!(violation_details, (
                                    branch_name = branch_names[branch_key],
                                    branch_key = branch_key,
                                    time_label = time_label,
                                    time_idx = time_idx,
                                    flow_value = flow_value,
                                    violation_percent = violation_percent
                                ))
                                
                                # Update branch violation statistics
                                stats = branch_violation_stats[branch_key]
                                
                                # Increment total violations
                                total_violations = stats.total_violations + 1
                                
                                # Increment current continuous violations
                                current_continuous = stats.current_continuous + 1
                                
                                # Update maximum continuous violations
                                max_continuous = max(stats.max_continuous_violations, current_continuous)
                                
                                # Update maximum violation percentage
                                max_viol_percent = stats.max_violation_percent
                                max_viol_time = stats.max_violation_time
                                
                                if violation_percent > max_viol_percent
                                    max_viol_percent = violation_percent
                                    max_viol_time = time_label
                                end
                                
                                # Update statistics
                                branch_violation_stats[branch_key] = (
                                    total_violations = total_violations,
                                    max_continuous_violations = max_continuous,
                                    current_continuous = current_continuous,
                                    max_flow = stats.max_flow,
                                    max_flow_time = stats.max_flow_time,
                                    max_violation_percent = max_viol_percent,
                                    max_violation_time = max_viol_time
                                )
                            else
                                # If no violation at current time point, reset continuous violation count
                                if branch_violation_stats[branch_key].current_continuous > 0
                                    branch_violation_stats[branch_key] = merge(
                                        branch_violation_stats[branch_key],
                                        (current_continuous = 0,)
                                    )
                                end
                            end
                        end
                    end
                end
            end
        end
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
    
    # Set chart size based on time span
    plot_size = time_day > 14 ? (900, 600) : (800, 500)
    
    # Create charts based on plot type
    if plot_type == "summary"
        # Create violation statistics chart
        plot_result = plot(layout=(2,1), size=plot_size)
        
        # Plot number of branches with violations
        plot!(subplot=1, time_points, violation_count, 
             title="Number of Branches with Violations",
             xlabel="",
             ylabel="Number of Violations",
             grid=true,
             linewidth=2,
             color=:red,
             label="")
        
        # Plot maximum violation percentage and total violation severity
        plot!(subplot=2, time_points, max_violation_percent, 
             title="Violation Severity",
             xlabel="Time",
             ylabel="Violation Percentage (%)",
             grid=true,
             linewidth=2,
             color=:blue,
             label="Maximum Violation Percentage")
        
        plot!(subplot=2, time_points, total_violation_severity, 
             linewidth=2,
             color=:purple,
             linestyle=:dash,
             label="Total Violation Severity")
        
        # Calculate statistics
        total_violations = sum(violation_count)
        max_violations = maximum(violation_count)
        max_violations_time_idx = argmax(violation_count)
        max_violation_pct = maximum(max_violation_percent)
        max_violation_pct_time_idx = argmax(max_violation_percent)
        
        # Get actual time label for maximum violation time
        max_viol_day = ceil(Int, max_violations_time_idx/24)
        max_viol_hour = (max_violations_time_idx-1)%24+1
        max_viol_time_label = "D$(max_viol_day)-H$(max_viol_hour)"
        
        # Get actual time label for maximum violation percentage time
        max_pct_day = ceil(Int, max_violation_pct_time_idx/24)
        max_pct_hour = (max_violation_pct_time_idx-1)%24+1
        max_pct_time_label = "D$(max_pct_day)-H$(max_pct_hour)"
        
        # Add statistics annotation
        stats_text = "Total Violations: $total_violations\n" *
                    "Max Branches with Violations: $max_violations ($(max_viol_time_label))\n" *
                    "Max Violation Percentage: $(round(max_violation_pct, digits=2))% ($(max_pct_time_label))"
        
        annotate!(subplot=1, 0.7 * length(time_points), 0.8 * max_violations, 
                 text(stats_text, :left, 8))
        
        # Set optimized x-axis tick labels
        plot!(subplot=1, xticks=([], []))  # No x-axis labels for top plot
        plot!(subplot=2, xticks=(xtick_indices, xtick_labels), xrotation=45)
        
        # Add auxiliary grid lines to make time points easier to match
        plot!(subplot=1, minorgrid=true, minorgridalpha=0.1)
        plot!(subplot=2, minorgrid=true, minorgridalpha=0.1)
        
    elseif plot_type == "worst"
        # Find branches with most severe violations
        branch_violation_severity = Dict()
        for branch_key in keys(branch_flows)
            flows = branch_flows[branch_key]
            violations = flows .> flow_limit
            if any(violations)
                # Calculate total violation severity for this branch
                violation_severity = sum((flows[violations] .- flow_limit) ./ flow_limit .* 100)
                branch_violation_severity[branch_key] = violation_severity
            end
        end
        
        # If no violations, display message
        if isempty(branch_violation_severity)
            plot_result = plot(title="No Branch Violations",
                              xlabel="Time",
                              ylabel="Power Flow (MW)",
                              grid=true,
                              size=plot_size)
            
            # Save the plot if a path is provided
            if save_path !== nothing
                # If save_path doesn't include file extension, add it based on save_format
                if !contains(save_path, ".")
                    save_path = "$(save_path).$(save_format)"
                end
                savefig(plot_result, save_path)
                println("Plot saved to: $(save_path)")
            end
            
            return plot_result, violation_count, max_violation_percent, total_violation_severity, violation_details, branch_violation_stats
        end
        
        # Sort by violation severity, select top 5 worst branches
        sorted_branches = sort(collect(branch_violation_severity), by=x->x[2], rev=true)
        top_branches = [branch for (branch, _) in sorted_branches[1:min(5, length(sorted_branches))]]
        
        # Plot power flows for worst violation branches
        plot_result = plot(title="Power Flows for Worst Violation Branches",
                          xlabel="Time",
                          ylabel="Power Flow (MW)",
                          grid=true,
                          size=plot_size,
                          legend=:topright)
        
        # Add power flow limit line
        hline!([flow_limit], linestyle=:dash, color=:black, label="Power Flow Limit ($flow_limit MW)")
        
        # Plot power flow curves for each major violation branch
        colors = [:red, :blue, :green, :purple, :orange]
        for (i, branch_key) in enumerate(top_branches)
            flows = branch_flows[branch_key]
            branch_name = branch_names[branch_key]
            
            # Calculate statistics
            max_flow = maximum(flows)
            max_flow_time_idx = argmax(flows)
            violation_times = count(flows .> flow_limit)
            avg_violation_pct = mean((flows[flows .> flow_limit] .- flow_limit) ./ flow_limit .* 100)
            
            # Get actual time label for maximum flow time
            max_flow_day = ceil(Int, max_flow_time_idx/24)
            max_flow_hour = (max_flow_time_idx-1)%24+1
            max_flow_time_label = "D$(max_flow_day)-H$(max_flow_hour)"
            
            # Simplify label to avoid legend being too long
            short_branch_name = length(branch_name) > 20 ? branch_name[1:18] * "..." : branch_name
            label_text = "$short_branch_name (Max: $(round(max_flow, digits=2)) MW)"
            
            plot!(time_points, flows, 
                 label=label_text, 
                 linewidth=2, 
                 color=colors[i])
        end
        
        # Set optimized x-axis tick labels
        plot!(xticks=(xtick_indices, xtick_labels), xrotation=45)
        
        # Add auxiliary grid lines to make time points easier to match
        plot!(minorgrid=true, minorgridalpha=0.1)
        
        # Add violation statistics information
        total_violations = sum(violation_count)
        violation_info = "Total violations: $total_violations times across $(length(branch_violation_severity)) branches"
        title!(plot_result, "Power Flows for Worst Violation Branches\n$violation_info")
        
    elseif plot_type == "all"
        # Find all branches with violations
        violated_branches = []
        for branch_key in keys(branch_flows)
            if any(violation_status[branch_key])
                push!(violated_branches, branch_key)
            end
        end
        
        # If no violations, display message
        if isempty(violated_branches)
            plot_result = plot(title="No Branch Violations",
                              xlabel="Time",
                              ylabel="Power Flow (MW)",
                              grid=true,
                              size=plot_size)
            
            # Save the plot if a path is provided
            if save_path !== nothing
                # If save_path doesn't include file extension, add it based on save_format
                if !contains(save_path, ".")
                    save_path = "$(save_path).$(save_format)"
                end
                savefig(plot_result, save_path)
                println("Plot saved to: $(save_path)")
            end
            
            return plot_result, violation_count, max_violation_percent, total_violation_severity, violation_details, branch_violation_stats
        end
        
        # Decide how to plot based on number of violation branches
        if length(violated_branches) <= 5
            # If not too many violation branches, show all in one plot
            plot_result = plot(title="Power Flows for All Violation Branches",
                              xlabel="Time",
                              ylabel="Power Flow (MW)",
                              grid=true,
                              size=plot_size,
                              legend=:topright)
            
            # Add power flow limit line
            hline!([flow_limit], linestyle=:dash, color=:black, label="Power Flow Limit ($flow_limit MW)")
            
            # Plot power flow curves for each violation branch
            colors = [:red, :blue, :green, :purple, :orange, :cyan, :magenta, :yellow]
            for (i, branch_key) in enumerate(violated_branches)
                flows = branch_flows[branch_key]
                branch_name = branch_names[branch_key]
                
                # Calculate statistics
                violation_times = count(flows .> flow_limit)
                
                # Simplify label to avoid legend being too long
                short_branch_name = length(branch_name) > 20 ? branch_name[1:18] * "..." : branch_name
                
                color_idx = ((i-1) % length(colors)) + 1
                plot!(time_points, flows, 
                     label="$short_branch_name (Viol: $violation_times)", 
                     linewidth=2, 
                     color=colors[color_idx])
            end
            
            # Set optimized x-axis tick labels
            plot!(xticks=(xtick_indices, xtick_labels), xrotation=45)
            
            # Add auxiliary grid lines to make time points easier to match
            plot!(minorgrid=true, minorgridalpha=0.1)
            
        else
            # Too many violation branches, create multiple subplots
            num_plots = min(9, length(violated_branches))
            plot_layout = num_plots <= 3 ? (num_plots, 1) : 
                         num_plots <= 6 ? (2, 3) : (3, 3)
            
            plot_result = plot(layout=plot_layout, size=(900, 700),
                              title="Branch Violation Power Flows")
            
            # Choose appropriate x-axis label interval for each subplot
            if time_day <= 7
                # For short time series, use the same label interval for each subplot
                subplot_xtick_indices = xtick_indices
                subplot_xtick_labels = xtick_labels
            else
                # For long time series, use fewer labels
                day_interval = max(1, round(Int, time_day / 5))
                subplot_xtick_indices = []
                subplot_xtick_labels = []
                
                for d in 1:day_interval:time_day
                    idx = (d-1)*24 + 1
                    push!(subplot_xtick_indices, idx)
                    push!(subplot_xtick_labels, "D$d")
                end
            end
            
            for i in 1:num_plots
                branch_key = violated_branches[i]
                flows = branch_flows[branch_key]
                branch_name = branch_names[branch_key]
                
                # Calculate statistics
                max_flow = maximum(flows)
                violation_times = count(flows .> flow_limit)
                
                # Simplify label to fit subplot
                short_branch_name = length(branch_name) > 25 ? branch_name[1:22] * "..." : branch_name
                
                # Plot power flow curve in subplot
                plot!(subplot=i, time_points, flows, 
                     title=short_branch_name,
                     xlabel=i > (plot_layout[1]-1)*plot_layout[2] ? "Time" : "",
                     ylabel="Power Flow (MW)",
                     grid=true,
                     linewidth=2,
                     color=:red,
                     label="")
                
                # Add power flow limit line
                hline!(subplot=i, [flow_limit], 
                      linestyle=:dash, 
                      color=:black, 
                      label="")
                
                # Add statistics
                annotate!(subplot=i, 0.7 * length(time_points), 0.8 * max_flow, 
                         text("Max: $(round(max_flow, digits=2)) MW\nViol: $violation_times", :left, 6))
                
                # Set x-axis labels for subplot
                if i > (plot_layout[1]-1)*plot_layout[2]  # Only show x-axis labels on bottom subplots
                    plot!(subplot=i, xticks=(subplot_xtick_indices, subplot_xtick_labels), xrotation=45)
                else
                    plot!(subplot=i, xticks=(subplot_xtick_indices, []))
                end
                
                # Add auxiliary grid lines to make time points easier to match
                plot!(subplot=i, minorgrid=true, minorgridalpha=0.1)
            end
            
            if length(violated_branches) > num_plots
                # Add message indicating more violation branches not shown
                remaining = length(violated_branches) - num_plots
                title!(plot_result, "Branch Violation Power Flows (Showing $num_plots of $(length(violated_branches)) branches)")
            end
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
    
    return plot_result, violation_count, max_violation_percent, total_violation_severity, violation_details, branch_violation_stats
end
