function ac_power_flow_compared(mpc, dict_bus, dict_new, etap_file_path)
    # 第一部分：计算并导出IEEE格式的潮流结果
    # 获取索引常量
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,
        VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()
    
    # 获取母线数据
    bus = mpc["bus"]
    
    # 创建反向映射字典
    reversed_dict_new = Dict(value => key for (key, value) in dict_new)
    reversed_dict = Dict(value => key for (key, value) in dict_bus)
    
    # 将母线ID映射回原始ID
    bus[:, BUS_I] = map(k -> reversed_dict_new[k], bus[:, BUS_I])
    bus_ID = map(k -> reversed_dict[k], bus[:, BUS_I])
    
    # 输出到Excel文件
    ieee_output_path = "output.xlsx"
    XLSX.openxlsx(ieee_output_path, mode="w") do workbook
        sheet = workbook["Sheet1"]
        # 写入表头
        sheet[1, 1] = "ID"
        sheet[1, 2] = "Voltage"
        sheet[1, 3] = "Angle"
        
        # 写入数据
        for i in eachindex(bus_ID)
            sheet[i + 1, 1] = bus_ID[i]  # 第一列写入ID
            sheet[i + 1, 2] = bus[i, VM]  # 第二列写入电压幅值
            sheet[i + 1, 3] = bus[i, VA]  # 第三列写入电压角度
        end
    end
    
    println("IEEE格式潮流结果已保存至 '$ieee_output_path'")
    
    # 第二部分：读取ETAP结果并进行比较
    # 读取ETAP结果文件
    etap_data = Dict{String, DataFrame}()
    XLSX.openxlsx(etap_file_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            data = sheet[:]
            etap_data[sheet_name] = DataFrame(data, :auto)
        end
    end
    
    # 提取ETAP结果数据
    etap_result = etap_data["Sheet1"]
    etap_bus_ID = etap_result[:, 1]
    etap_voltage = etap_result[:, 2]
    etap_angle = etap_result[:, 3]
    
    # 创建ETAP结果字典
    etap_dict_voltage = Dict(etap_bus_ID .=> etap_voltage)
    etap_dict_angle = Dict(etap_bus_ID .=> etap_angle)
    
    # 读取IEEE格式结果
    ieee_data = Dict{String, DataFrame}()
    XLSX.openxlsx(ieee_output_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            data = sheet[:]
            ieee_data[sheet_name] = DataFrame(data, :auto)
        end
    end
    
    # 提取IEEE格式结果数据
    ieee_result = ieee_data["Sheet1"]
    ieee_bus_ID = ieee_result[2:end, 1]
    ieee_voltage = ieee_result[2:end, 2]
    ieee_angle = ieee_result[2:end, 3]
    
    # 计算误差
    voltage_diff = ieee_voltage .- map(k -> etap_dict_voltage[k], ieee_bus_ID)
    angle_diff = ieee_angle .- map(k -> etap_dict_angle[k], ieee_bus_ID)
    
    # 计算相对误差
    voltage_rel_error = voltage_diff ./ map(k -> etap_dict_voltage[k], ieee_bus_ID)
    angle_rel_error = angle_diff ./ map(k -> etap_dict_angle[k], ieee_bus_ID)
    
    # 计算平均绝对误差
    avg_abs_voltage_error = sum(abs.(voltage_diff)) / length(voltage_diff)
    avg_abs_angle_error = sum(abs.(angle_diff)) / length(angle_diff)
    
    # 计算最大绝对误差
    max_abs_voltage_error = maximum(abs.(voltage_diff))
    max_abs_angle_error = maximum(abs.(angle_diff))
    
    # 输出比较结果到Excel文件
    comparison_output_path = "comparison_results.xlsx"
    XLSX.openxlsx(comparison_output_path, mode="w") do workbook
        sheet = workbook["Sheet1"]
        
        # 写入表头
        sheet[1, 1] = "Bus ID"
        sheet[1, 2] = "IEEE Voltage"
        sheet[1, 3] = "ETAP Voltage"
        sheet[1, 4] = "Voltage Diff"
        sheet[1, 5] = "Voltage Rel Error"
        sheet[1, 6] = "IEEE Angle"
        sheet[1, 7] = "ETAP Angle"
        sheet[1, 8] = "Angle Diff"
        sheet[1, 9] = "Angle Rel Error"
        
        # 写入数据
        for i in eachindex(ieee_bus_ID)
            etap_v = etap_dict_voltage[ieee_bus_ID[i]]
            etap_a = etap_dict_angle[ieee_bus_ID[i]]
            
            sheet[i + 1, 1] = ieee_bus_ID[i]
            sheet[i + 1, 2] = ieee_voltage[i]
            sheet[i + 1, 3] = etap_v
            sheet[i + 1, 4] = voltage_diff[i]
            sheet[i + 1, 5] = voltage_rel_error[i]
            sheet[i + 1, 6] = ieee_angle[i]
            sheet[i + 1, 7] = etap_a
            sheet[i + 1, 8] = angle_diff[i]
            sheet[i + 1, 9] = angle_rel_error[i]
        end
        
        # 添加统计信息
        row_num = length(ieee_bus_ID) + 3
        sheet[row_num, 1] = "Statistics"
        sheet[row_num + 1, 1] = "Avg Abs Error"
        sheet[row_num + 1, 4] = avg_abs_voltage_error
        sheet[row_num + 1, 8] = avg_abs_angle_error
        sheet[row_num + 2, 1] = "Max Abs Error"
        sheet[row_num + 2, 4] = max_abs_voltage_error
        sheet[row_num + 2, 8] = max_abs_angle_error
    end
    
    println("比较结果已保存至 '$comparison_output_path'")
    
    # 同时绘制电压幅值和相角的相对误差图
    p = Plots.plot(
        title="Comparison with ETAP Results",
        xlabel="Bus Index",
        ylabel="Relative Error",
        legend=:topright,
        size=(800, 500)
    )
    
    # 绘制电压幅值相对误差曲线（不使用点标记）
    Plots.plot!(p, 1:length(voltage_rel_error), voltage_rel_error, 
        linewidth=2,
        color=:blue,
        label="Voltage Magnitude"
    )
    
    # 绘制相角相对误差曲线（不使用点标记）
    Plots.plot!(p, 1:length(angle_rel_error), angle_rel_error, 
        linewidth=2,
        color=:red,
        label="Voltage Angle"
    )
    
    # 显示图表
    display(p)
    
    # 返回相对误差
    return voltage_rel_error, angle_rel_error, avg_abs_voltage_error, avg_abs_angle_error
end
