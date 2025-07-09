function plot_losses_time_series(results, case, time_day, plot_type = "total", loss_type = "active"; save_path = nothing, save_format = "pdf")
    default(fontfamily="Microsoft YaHei")
    # Get number of islands
    num_islands = size(results, 1)
    
    # Create arrays to store total losses
    total_ac_active_losses = zeros(time_day*24)
    total_ac_reactive_losses = zeros(time_day*24)
    total_dc_active_losses = zeros(time_day*24)
    
    # Process each island separately
    for island in 1:num_islands
        # Get number of AC branches in this island
        num_branch_AC = 0
        num_branch_AC = size(results[island, 1, 1].branchAC, 1)
        
        # Get number of DC branches in this island
        num_branch_DC = 0
        num_branch_DC = size(results[island, 1, 1].branchDC, 1)

        
        # Create time series loss matrices - AC branches
        ac_active_losses = zeros(num_branch_AC, time_day*24 + 1)
        ac_reactive_losses = zeros(num_branch_AC, time_day*24 + 1)
        ac_active_losses[:, 1] = 1:num_branch_AC
        ac_reactive_losses[:, 1] = 1:num_branch_AC
        
        # Create time series loss matrices - DC branches
        dc_active_losses = zeros(num_branch_DC, time_day*24 + 1)
        dc_active_losses[:, 1] = 1:num_branch_DC
        
        # Extract branch power flows from results and calculate losses for each time point
        for d in 1:time_day
            for hour in 1:24
                time_idx = (d-1)*24 + hour
                
                # Process AC branches
                if num_branch_AC > 0 
                    branch_data = results[island, d, hour].branchAC
                    
                    # Ensure branch data contains necessary columns
                    if size(branch_data, 2) >= QT
                        # Calculate losses for each branch
                        for b in 1:num_branch_AC
                            # Calculate active power loss, considering power flow direction
                            if branch_data[b, PF] * branch_data[b, PT] <= 0  # Power flow direction consistent
                                # Determine power flow direction
                                if abs(branch_data[b, PF]) >= abs(branch_data[b, PT])
                                    # From start node to end node
                                    p_loss = abs(branch_data[b, PF]) - abs(branch_data[b, PT])
                                else
                                    # From end node to start node
                                    p_loss = abs(branch_data[b, PT]) - abs(branch_data[b, PF])
                                end
                            else  # Power flow direction inconsistent (possibly a ring network)
                                p_loss = abs(branch_data[b, PF] + branch_data[b, PT])
                            end
                            
                            ac_active_losses[b, time_idx+1] = p_loss
                            total_ac_active_losses[time_idx] += p_loss
                            
                            # Calculate reactive power loss, considering power flow direction
                            if branch_data[b, QF] * branch_data[b, QT] <= 0  # Power flow direction consistent
                                # Determine power flow direction
                                if abs(branch_data[b, QF]) >= abs(branch_data[b, QT])
                                    # From start node to end node
                                    q_loss = abs(branch_data[b, QF]) - abs(branch_data[b, QT])
                                else
                                    # From end node to start node
                                    q_loss = abs(branch_data[b, QT]) - abs(branch_data[b, QF])
                                end
                            else  # Power flow direction inconsistent (possibly a ring network)
                                q_loss = abs(branch_data[b, QF] + branch_data[b, QT])
                            end
                            
                            ac_reactive_losses[b, time_idx+1] = q_loss
                            total_ac_reactive_losses[time_idx] += q_loss
                        end
                    end
                end
                
                # Process DC branches
                if num_branch_DC > 0 
                    dc_branch_data = results[island, d, hour].branchDC
                    
                    # Ensure DC branch data contains necessary columns (assuming power flow is in similar positions)
                    if size(dc_branch_data, 2) >= PT
                        # Calculate losses for each DC branch
                        for b in 1:num_branch_DC
                            # For DC branches, active power loss calculation needs to consider power flow direction
                            if dc_branch_data[b, PF] * dc_branch_data[b, PT] <= 0  # Power flow direction consistent
                                # Determine power flow direction
                                if abs(dc_branch_data[b, PF]) >= abs(dc_branch_data[b, PT])
                                    # From start node to end node
                                    p_loss = abs(dc_branch_data[b, PF]) - abs(dc_branch_data[b, PT])
                                else
                                    # From end node to start node
                                    p_loss = abs(dc_branch_data[b, PT]) - abs(dc_branch_data[b, PF])
                                end
                            else  # Power flow direction inconsistent (possibly a ring network)
                                p_loss = abs(dc_branch_data[b, PF] + dc_branch_data[b, PT])
                            end
                            
                            dc_active_losses[b, time_idx+1] = p_loss
                            total_dc_active_losses[time_idx] += p_loss
                        end
                    end
                end
            end
        end
    end
    
    # 改进的时间轴标签处理
    time_labels = String[]
    
    # 根据总天数选择合适的标签格式
    if time_day <= 3  # 3天或更短
        # 每4小时一个标签
        for d in 1:time_day
            for h in 1:24
                if h % 4 == 0 || h == 1
                    push!(time_labels, "D$(d)-H$(h)")
                else
                    push!(time_labels, "")
                end
            end
        end
    elseif time_day <= 7  # 一周或更短
        # 每6小时一个标签
        for d in 1:time_day
            for h in 1:24
                if h % 6 == 0 || h == 1
                    push!(time_labels, "D$(d)-H$(h)")
                else
                    push!(time_labels, "")
                end
            end
        end
    elseif time_day <= 31  # 一个月或更短
        # 每天只显示一个标签（第一个小时）
        for d in 1:time_day
            for h in 1:24
                if h == 1
                    push!(time_labels, "D$(d)")
                else
                    push!(time_labels, "")
                end
            end
        end
    else  # 超过一个月
        # 每隔几天显示一个标签
        label_interval = max(1, round(Int, time_day / 20))  # 动态计算标签间隔，确保总标签数不超过20个左右
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
    
    # 创建时间点索引数组，用于指定哪些点显示标签
    time_points = 1:time_day*24
    xtick_indices = findall(x -> x != "", time_labels)
    xtick_labels = time_labels[xtick_indices]
    
    # 设置图表大小，根据时间跨度调整
    plot_size = time_day > 14 ? (900, 600) : (800, 500)
    
    # Create charts based on loss type and plot type
    if plot_type == "total"
        plot_result = plot(title="System Losses Time Series",
                          xlabel="Time",
                          ylabel="Losses",
                          grid=true,
                          size=plot_size,
                          legend=:topright)
        
        if loss_type == "active" || loss_type == "both"
            # Plot total active losses
            total_active_losses = total_ac_active_losses + total_dc_active_losses
            plot!(time_points, total_active_losses, 
                 label="Total Active Losses", 
                 linewidth=2, 
                 color=:red)
            
            # If there are DC branch losses, show them separately
            if any(total_dc_active_losses .!= 0)
                plot!(time_points, total_ac_active_losses, 
                     label="AC Active Losses", 
                     linewidth=2, 
                     linestyle=:dash,
                     color=:darkred)
                
                plot!(time_points, total_dc_active_losses, 
                     label="DC Active Losses", 
                     linewidth=2, 
                     linestyle=:dash,
                     color=:orange)
            end
            
            # Calculate active loss statistics
            avg_active_loss = mean(total_active_losses)
            max_active_loss = maximum(total_active_losses)
            min_active_loss = minimum(total_active_losses)
            max_active_loss_time_idx = argmax(total_active_losses)
            min_active_loss_time_idx = argmin(total_active_losses)
            
            # 获取最大/最小损耗时间的实际时间标签
            max_day = ceil(Int, max_active_loss_time_idx/24)
            max_hour = (max_active_loss_time_idx-1)%24+1
            max_time_label = "D$(max_day)-H$(max_hour)"
            
            min_day = ceil(Int, min_active_loss_time_idx/24)
            min_hour = (min_active_loss_time_idx-1)%24+1
            min_time_label = "D$(min_day)-H$(min_hour)"
            
            # Add active loss statistics annotation
            active_stats_text = "Average Active Loss: $(round(avg_active_loss, digits=2)) MW\n" *
                               "Maximum Active Loss: $(round(max_active_loss, digits=2)) MW ($(max_time_label))\n" *
                               "Minimum Active Loss: $(round(min_active_loss, digits=2)) MW ($(min_time_label))"
            
            if loss_type == "active"
                annotate!(0.7 * length(total_active_losses), 
                         min_active_loss + 0.8 * (max_active_loss - min_active_loss), 
                         text(active_stats_text, :left, 8))
            end
        end
        
        if loss_type == "reactive" || loss_type == "both"
            # Plot total reactive losses
            plot!(time_points, total_ac_reactive_losses, 
                 label="Total Reactive Losses", 
                 linewidth=2, 
                 color=:blue)
            
            # Calculate reactive loss statistics
            avg_reactive_loss = mean(total_ac_reactive_losses)
            max_reactive_loss = maximum(total_ac_reactive_losses)
            min_reactive_loss = minimum(total_ac_reactive_losses)
            max_reactive_loss_time_idx = argmax(total_ac_reactive_losses)
            min_reactive_loss_time_idx = argmin(total_ac_reactive_losses)
            
            # 获取最大/最小损耗时间的实际时间标签
            max_day = ceil(Int, max_reactive_loss_time_idx/24)
            max_hour = (max_reactive_loss_time_idx-1)%24+1
            max_time_label = "D$(max_day)-H$(max_hour)"
            
            min_day = ceil(Int, min_reactive_loss_time_idx/24)
            min_hour = (min_reactive_loss_time_idx-1)%24+1
            min_time_label = "D$(min_day)-H$(min_hour)"
            
            # Add reactive loss statistics annotation
            reactive_stats_text = "Average Reactive Loss: $(round(avg_reactive_loss, digits=2)) MVAr\n" *
                                "Maximum Reactive Loss: $(round(max_reactive_loss, digits=2)) MVAr ($(max_time_label))\n" *
                                "Minimum Reactive Loss: $(round(min_reactive_loss, digits=2)) MVAr ($(min_time_label))"
            
            if loss_type == "reactive"
                annotate!(0.7 * length(total_ac_reactive_losses), 
                         min_reactive_loss + 0.8 * (max_reactive_loss - min_reactive_loss), 
                         text(reactive_stats_text, :left, 8))
            end
        end
        
        # If showing both active and reactive losses, combine statistics
        if loss_type == "both"
            total_active_losses = total_ac_active_losses + total_dc_active_losses
            combined_stats_text = "Average Active Loss: $(round(mean(total_active_losses), digits=2)) MW\n" *
                                "Average Reactive Loss: $(round(mean(total_ac_reactive_losses), digits=2)) MVAr"
            
            annotate!(0.7 * length(total_active_losses), 
                     0.8 * maximum(total_active_losses), 
                     text(combined_stats_text, :left, 8))
        end
        
    else  # "branch" - Plot branch losses
        # This part needs to be modified since we now have multiple islands and different types of branches
        # We will only show the branches with the highest losses in the first island as an example
        
        # Get number of AC branches in the first island
        num_branch_AC_island1 = 0
        if size(results, 1) >= 1 && haskey(results[1, 1, 1], :branchAC) && !isempty(results[1, 1, 1].branchAC)
            num_branch_AC_island1 = size(results[1, 1, 1].branchAC, 1)
        end
        
        if num_branch_AC_island1 > 0
            # Create time series loss matrices - AC branches in first island
            island1_ac_active_losses = zeros(num_branch_AC_island1, time_day*24 + 1)
            island1_ac_reactive_losses = zeros(num_branch_AC_island1, time_day*24 + 1)
            island1_ac_active_losses[:, 1] = 1:num_branch_AC_island1
            island1_ac_reactive_losses[:, 1] = 1:num_branch_AC_island1
            
            # Fill in loss data for the first island
            for d in 1:time_day
                for hour in 1:24
                    time_idx = (d-1)*24 + hour
                    
                    if haskey(results[1, d, hour], :branchAC)
                        branch_data = results[1, d, hour].branchAC
                        
                        if size(branch_data, 2) >= QT
                            for b in 1:num_branch_AC_island1
                                # Calculate active power loss, considering power flow direction
                                if branch_data[b, PF] * branch_data[b, PT] <= 0  # Power flow direction consistent
                                    # Determine power flow direction
                                    if abs(branch_data[b, PF]) >= abs(branch_data[b, PT])
                                        # From start node to end node
                                        p_loss = abs(branch_data[b, PF]) - abs(branch_data[b, PT])
                                    else
                                        # From end node to start node
                                        p_loss = abs(branch_data[b, PT]) - abs(branch_data[b, PF])
                                    end
                                else  # Power flow direction inconsistent (possibly a ring network)
                                    p_loss = abs(branch_data[b, PF] + branch_data[b, PT])
                                end
                                
                                island1_ac_active_losses[b, time_idx+1] = p_loss
                                
                                # Calculate reactive power loss, considering power flow direction
                                if branch_data[b, QF] * branch_data[b, QT] <= 0  # Power flow direction consistent
                                    # Determine power flow direction
                                    if abs(branch_data[b, QF]) >= abs(branch_data[b, QT])
                                        # From start node to end node
                                        q_loss = abs(branch_data[b, QF]) - abs(branch_data[b, QT])
                                    else
                                        # From end node to start node
                                        q_loss = abs(branch_data[b, QT]) - abs(branch_data[b, QF])
                                    end
                                else  # Power flow direction inconsistent (possibly a ring network)
                                    q_loss = abs(branch_data[b, QF] + branch_data[b, QT])
                                end
                                
                                island1_ac_reactive_losses[b, time_idx+1] = q_loss
                            end
                        end
                    end
                end
            end
            
            # Select top 5 branches with highest losses for plotting
            if loss_type == "active" || loss_type == "both"
                avg_branch_active_losses = [mean(island1_ac_active_losses[b, 2:end]) for b in 1:num_branch_AC_island1]
                top_branches_idx_active = sortperm(avg_branch_active_losses, rev=true)[1:min(5, num_branch_AC_island1)]
                
                plot_result = plot(title="Major Branch Active Losses Time Series",
                                  xlabel="Time",
                                  ylabel="Active Losses (MW)",
                                  grid=true,
                                  size=plot_size,
                                  legend=:topright)
                
                # Plot active loss curves for each major branch
                colors = [:red, :blue, :green, :purple, :orange]
                for (i, b) in enumerate(top_branches_idx_active)
                    branch_losses = island1_ac_active_losses[b, 2:end]
                    
                    # Get branch from and to buses
                    from_bus = Int(results[1, 1, 1].branchAC[b, F_BUS])
                    to_bus = Int(results[1, 1, 1].branchAC[b, T_BUS])
                    
                    # Try to get bus names
                    from_name = haskey(case.id_to_bus_name, from_bus) ? case.id_to_bus_name[from_bus] : "Bus$from_bus"
                    to_name = haskey(case.id_to_bus_name, to_bus) ? case.id_to_bus_name[to_bus] : "Bus$to_bus"
                    
                    branch_name = "$from_name - $to_name"
                    avg_loss = round(mean(branch_losses), digits=2)
                    
                    plot!(time_points, branch_losses, 
                         label="$branch_name (Avg: $avg_loss MW)", 
                         linewidth=2, 
                         color=colors[i])
                end
            end
            
            # If reactive losses need to be plotted
            if loss_type == "reactive"
                avg_branch_reactive_losses = [mean(island1_ac_reactive_losses[b, 2:end]) for b in 1:num_branch_AC_island1]
                top_branches_idx_reactive = sortperm(avg_branch_reactive_losses, rev=true)[1:min(5, num_branch_AC_island1)]
                
                plot_result = plot(title="Major Branch Reactive Losses Time Series",
                                  xlabel="Time",
                                  ylabel="Reactive Losses (MVAr)",
                                  grid=true,
                                  size=plot_size,
                                  legend=:topright)
                
                # Plot reactive loss curves for each major branch
                colors = [:blue, :red, :green, :purple, :orange]
                for (i, b) in enumerate(top_branches_idx_reactive)
                    branch_losses = island1_ac_reactive_losses[b, 2:end]
                    
                    # Get branch from and to buses
                    from_bus = Int(results[1, 1, 1].branchAC[b, F_BUS])
                    to_bus = Int(results[1, 1, 1].branchAC[b, T_BUS])
                    
                    # Try to get bus names
                    from_name = haskey(case.id_to_bus_name, from_bus) ? case.id_to_bus_name[from_bus] : "Bus$from_bus"
                    to_name = haskey(case.id_to_bus_name, to_bus) ? case.id_to_bus_name[to_bus] : "Bus$to_bus"
                    
                    branch_name = "$from_name - $to_name"
                    avg_loss = round(mean(branch_losses), digits=2)
                    
                    plot!(time_points, branch_losses, 
                         label="$branch_name (Avg: $avg_loss MVAr)", 
                         linewidth=2, 
                         color=colors[i])
                end
            end
        else
            plot_result = plot(title="No Branch Data Available",
                              xlabel="Time",
                              ylabel="Losses",
                              grid=true,
                              size=plot_size)
        end
    end
    
    # 设置优化后的x轴刻度标签
    plot_result = plot!(plot_result, xticks=(xtick_indices, xtick_labels), xrotation=45)
    
    # 添加辅助网格线，使时间点更容易对应
    plot_result = plot!(plot_result, minorgrid=true, minorgridalpha=0.1)
    
    # 添加日期范围标注
    if time_day > 7
        day_range_text = "Time span: $(time_day) days"
        if plot_type == "total"
            if loss_type == "active"
                total_losses = total_ac_active_losses + total_dc_active_losses
                annotate!(0.5*length(time_points), minimum(total_losses) - 0.05*(maximum(total_losses)-minimum(total_losses)), 
                         text(day_range_text, :center, 8))
            elseif loss_type == "reactive"
                annotate!(0.5*length(time_points), minimum(total_ac_reactive_losses) - 0.05*(maximum(total_ac_reactive_losses)-minimum(total_ac_reactive_losses)), 
                         text(day_range_text, :center, 8))
            else # both
                total_losses = total_ac_active_losses + total_dc_active_losses
                annotate!(0.5*length(time_points), minimum(total_losses) - 0.05*(maximum(total_losses)-minimum(total_losses)), 
                         text(day_range_text, :center, 8))
            end
        else # branch
            # 为分支损耗图添加时间跨度标注
            annotate!(0.5*length(time_points), 0, text(day_range_text, :center, 8))
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
    
    # Return plot result and loss data
    return plot_result, total_ac_active_losses, total_ac_reactive_losses, total_dc_active_losses
end
