"""
绘制有功负荷的时间序列图
参数:
results - 结果数据集
bus_name - 母线名称
case - 系统案例
time_day - 天数
bus_type - 母线类型 (默认为 "AC")
"""
function plot_PD_time_series(results, bus_name, case, time_day, bus_type = "AC")
    default(fontfamily="Microsoft YaHei")
    # 创建一个时间序列有功负荷矩阵，第一列表示节点序号
    num_bus_AC = length(case.busesAC)
    PD_time_series_AC = zeros(num_bus_AC, time_day*24 + 1)  # 第一列为节点序号，后续列为每个时间点的有功负荷
    PD_time_series_AC[:, 1] = 1:num_bus_AC  # 第一列为节点序号
    
    for i in 1:size(results,1)
        for d in 1:time_day
            for hour in 1:24
                # 获取当前时间点的有功负荷
                PD_time_series_AC[Int.(results[i, d, hour].busAC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busAC[:, PD]  # 读取PD数据
            end
        end
    end

    # 如果需要DC母线的有功负荷，也可以类似处理
    num_bus_DC = length(case.busesDC)
    PD_time_series_DC = zeros(num_bus_DC, time_day*24 + 1)  # 第一列为节点序号，后续列为每个时间点的有功负荷
    PD_time_series_DC[:, 1] = 1:num_bus_DC  # 第一列为节点序号
    
    for i in 1:size(results,1)
        for d in 1:time_day
            for hour in 1:24
                # 获取当前时间点的有功负荷
                if size(results[i, d, hour].busDC, 2) >= PD  # 确保DC母线有PD列
                    PD_time_series_DC[Int.(results[i, d, hour].busDC[:, BUS_I]),(d-1)*24 + hour+1] .= results[i, d, hour].busDC[:, PD]
                end
            end
        end
    end

    # 创建时间轴标签
    time_labels = String[]
    for d in 1:time_day
        for h in 1:24
            push!(time_labels, "D$(d)-H$(h)")
        end
    end
    
    if bus_type == "AC"
        bus_name_to_index = case.bus_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        PD_series = PD_time_series_AC[voltage_row, 2:end]  # 获取对应母线的有功负荷

        # 创建时间点数组
        time_points = 1:length(PD_series)
        
        # 创建图表
        plot_result = plot(time_points, PD_series, 
                         label="", 
                         linewidth=2, 
                         color=:blue,
                         title="AC Bus \"$(bus_name)\" Active Load Time Series",
                         xlabel="Time",
                         ylabel="Active Load (MW)",
                         grid=true)
        
        # 设置x轴刻度标签
        if time_day * 24 <= 48  # 如果总时间点不超过48个，则全部显示
            xtick_indices = 1:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels))
        else  # 否则每6个小时显示一个标签
            xtick_indices = 1:6:length(PD_series)
            plot_result = plot!(plot_result, xticks=(xtick_indices, time_labels[xtick_indices]))
        end
        
    else  # DC bus
        bus_name_to_index = case.busdc_name_to_id
        voltage_row = bus_name_to_index[bus_name]
        PD_series = PD_time_series_DC[voltage_row, 2:end]  # 获取对应母线的有功负荷

        # 创建时间点数组
        time_points = 1:length(PD_series)
        
        # 绘制有功负荷时间序列图
        plot_result = plot(time_points, PD_series, 
                          label="", 
                          linewidth=2, 
                          color=:green,
                          title="DC Bus \"$(bus_name)\" Active Load Time Series",
                          xlabel="Time",
                          ylabel="Active Load (MW)",
                          grid=true)
        
        # 设置x轴刻度标签
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