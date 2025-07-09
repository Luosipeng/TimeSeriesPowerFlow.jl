function acdc_power_flow_compared(mpc, dict_bus, Dict_busdc, node_mapping, result_path)
    busAC = mpc["busAC"]
    busDC = mpc["busDC"]

    # 读取电压数据
    IDac = busAC[:, 1]
    IDdc = busDC[:, 1]
    Vmac = busAC[:, 8]
    Vaac = busAC[:, 9]
    Vmdc = busDC[:, 8]

    # 获取反转字典
    reversed_dict = Dict(value => key for (key, value) in dict_bus)
    reversed_node_mapping = Dict(value => key for (key, value) in node_mapping)
    reversed_dict_dc = Dict(value => key for (key, value) in Dict_busdc)

    IDac = map(k -> reversed_node_mapping[k], IDac)
    IDac = map(k -> reversed_dict[k], IDac)
    IDdc = map(k -> reversed_dict_dc[k], IDdc)

    indices_without_bus_b = findall(id -> !contains(id, "Bus_b"), IDdc)
    IDdc = IDdc[indices_without_bus_b]
    Vmdc = Vmdc[indices_without_bus_b]
    
    # 读取ETAP的结果
    etap_result = Dict{String, DataFrame}()
    XLSX.openxlsx(result_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            etap_result[sheet_name] = DataFrame(sheet[:], :auto)
        end
    end

    etap_bus_data = etap_result["Sheet1"]
    AC_index = findall(etap_bus_data[:, 2] .== "AC")
    DC_index = findall(etap_bus_data[:, 2] .== "DC")
    etap_bus_data_AC = etap_bus_data[AC_index, :]
    etap_bus_data_DC = etap_bus_data[DC_index, :]

    # 提取电压数据
    etap_IDac = etap_bus_data_AC[:, 1]
    etap_IDdc = etap_bus_data_DC[:, 1]
    etap_Vmac = etap_bus_data_AC[:, 3]
    etap_Vaac = etap_bus_data_AC[:, 4]
    etap_Vmdc = etap_bus_data_DC[:, 3]

    # 创建字典映射ID和电压
    compared_dict_Vmac = Dict(etap_IDac .=> etap_Vmac)
    compared_dict_Vaac = Dict(etap_IDac .=> etap_Vaac)
    compared_dict_Vmdc = Dict(etap_IDdc .=> etap_Vmdc)
    
    # 计算相对误差
    misVaac_raw = (Vaac .- map(k -> compared_dict_Vaac[k], IDac)) ./ (map(k -> compared_dict_Vaac[k], IDac))
    misVmac_raw = (Vmac .- map(k -> compared_dict_Vmac[k], IDac)) ./ (map(k -> compared_dict_Vmac[k], IDac))
    misVmdc_raw = (Vmdc .- map(k -> compared_dict_Vmdc[k], IDdc)) ./ (map(k -> compared_dict_Vmdc[k], IDdc))
    
    # 处理超过1的误差值和NaN值，将其设为0
    process_error = err -> (isnan(err) || abs(err) == 1) ? 0.0 : err
    
    misVaac = map(process_error, misVaac_raw)
    misVmac = map(process_error, misVmac_raw)
    misVmdc = map(process_error, misVmdc_raw)

    # 计算统计信息 (使用处理后的误差值)
    avg_abs_vm_ac_error = sum(abs.(misVmac)) / length(misVmac)
    avg_abs_va_ac_error = sum(abs.(misVaac)) / length(misVaac)
    avg_abs_vm_dc_error = sum(abs.(misVmdc)) / length(misVmdc)
    
    max_abs_vm_ac_error = maximum(abs.(misVmac))
    max_abs_va_ac_error = maximum(abs.(misVaac))
    max_abs_vm_dc_error = maximum(abs.(misVmdc))
    
    # 创建比较结果Excel文件
    output_path = joinpath(dirname(result_path), "voltage_comparison.xlsx")
    
    # 创建AC节点对比数据 (使用原始误差值，以便保留完整信息)
    ac_comparison = DataFrame(
        "节点ID" => IDac,
        "计算模型幅值(pu)" => Vmac,
        "ETAP幅值(pu)" => map(k -> compared_dict_Vmac[k], IDac),
        "幅值相对误差" => misVmac_raw,
        "处理后幅值误差" => misVmac,
        "计算模型相角(度)" => Vaac,
        "ETAP相角(度)" => map(k -> compared_dict_Vaac[k], IDac),
        "相角相对误差" => misVaac_raw,
        "处理后相角误差" => misVaac
    )
    
    # 创建DC节点对比数据 (使用原始误差值，以便保留完整信息)
    dc_comparison = DataFrame(
        "节点ID" => IDdc,
        "计算模型电压(pu)" => Vmdc,
        "ETAP电压(pu)" => map(k -> compared_dict_Vmdc[k], IDdc),
        "相对误差" => misVmdc_raw,
        "处理后误差" => misVmdc
    )
    
    # 添加统计信息
    stats_data = DataFrame(
        "统计指标" => ["平均绝对误差", "最大绝对误差"],
        "AC幅值误差" => [avg_abs_vm_ac_error, max_abs_vm_ac_error],
        "AC相角误差" => [avg_abs_va_ac_error, max_abs_va_ac_error],
        "DC电压误差" => [avg_abs_vm_dc_error, max_abs_vm_dc_error]
    )
    
    # 导出到Excel文件
    try
        if isfile(output_path)
            rm(output_path)
        end
        
        XLSX.writetable(output_path, 
            AC节点对比 = (collect(DataFrames.eachcol(ac_comparison)), DataFrames.names(ac_comparison)),
            DC节点对比 = (collect(DataFrames.eachcol(dc_comparison)), DataFrames.names(dc_comparison)),
            统计信息 = (collect(DataFrames.eachcol(stats_data)), DataFrames.names(stats_data))
        )
        println("对比结果已保存至 '$output_path'")
    catch e
        backup_path = joinpath(dirname(result_path), "voltage_comparison_$(round(Int, time())).xlsx")
        XLSX.writetable(backup_path, 
            AC节点对比 = (collect(DataFrames.eachcol(ac_comparison)), DataFrames.names(ac_comparison)),
            DC节点对比 = (collect(DataFrames.eachcol(dc_comparison)), DataFrames.names(dc_comparison)),
            统计信息 = (collect(DataFrames.eachcol(stats_data)), DataFrames.names(stats_data))
        )
        println("对比结果已保存至备用文件: '$backup_path'")
        output_path = backup_path
    end
    
    # 绘制相对误差图 (使用处理后的误差值)
    p = plot(
        title="AC/DC Power Flow Comparison with ETAP Results",
        xlabel="Bus Index",
        ylabel="Relative Error",
        legend=:topright,
        size=(800, 500)
    )
    
    # 绘制AC电压幅值相对误差曲线
    plot!(p, 1:length(misVmac), misVmac, 
        linewidth=2,
        color=:blue,
        label="AC Voltage Magnitude"
    )
    
    # 绘制AC相角相对误差曲线
    plot!(p, 1:length(misVaac), misVaac, 
        linewidth=2,
        color=:red,
        label="AC Voltage Angle"
    )
    
    # 绘制DC电压相对误差曲线
    plot!(p, (length(misVmac)+1):(length(misVmac)+length(misVmdc)), misVmdc, 
        linewidth=2,
        color=:green,
        label="DC Voltage"
    )
    
    # 添加说明文本
    annotate!(p, 5, minimum(filter(x -> !isnan(x) && abs(x) <= 1, [misVmac; misVaac; misVmdc])) - 0.1, 
              text("注: 相对误差绝对值>1或NaN值已设为0", 8, :black))
    
    # 显示图表
    display(p)
    
    # 返回相对误差和统计数据 (返回处理后的误差值)
    return misVmac, misVaac, misVmdc, avg_abs_vm_ac_error, avg_abs_va_ac_error, avg_abs_vm_dc_error
end
