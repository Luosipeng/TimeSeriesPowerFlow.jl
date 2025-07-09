function ac_element_validate(elementtype, field, dict_bus, node_mapping, result_path, opt, mpc; connectedname=nothing, Fbus=nothing, Tbus=nothing)
    bus_ID = []
    voltages = []
    angles = []  # 新增存储相角的数组
    value = []
    success_status = [] # 存储成功状态

    if elementtype == "LUMPLOAD" && field == "MW"
        if connectedname !== nothing
            aim_bus = findall(i -> mpc["bus"][i, 1] == connectedname, 1:size(mpc["bus"], 1))
            aim_load = findall(i -> mpc["load"][i, 2] == connectedname, 1:size(mpc["load"], 1))
        end
        start_value = 0.04
        for i in 0:9
            current_value = start_value + i * 0.01
            mpc["load"][aim_load, 4] .= current_value*0.85
            mpc["load"][aim_load, 5] .= current_value*sqrt(1-0.85^2)
            mpc["bus"][aim_bus, 3] .= current_value*0.85
            mpc["bus"][aim_bus, 4] .= current_value*sqrt(1-0.85^2)
            result = PowerFlow.runpf(mpc, opt)
            
            # 输出并存储潮流计算是否成功
            is_success = result["success"]
            push!(success_status, is_success)
            println("计算 $(i+1)/10: 负荷值 = $(current_value*1000) KW, 潮流计算 $(is_success ? "成功" : "失败")")
            
            append!(bus_ID, result["bus"][aim_bus, 1])
            append!(voltages, result["bus"][aim_bus, 8])
            append!(angles, result["bus"][aim_bus, 9])  # 存储相角
            append!(value, current_value * 1000)
        end

    elseif elementtype == "XFORM3W" && field == "PrimPercentTap"
        start_value = 0.0
        for i in 0:9
            current_value = start_value + i * 0.2
            mpc["branch"][1, 9] = current_value/100 + 1.0
            result = PowerFlow.runpf(mpc, opt)
            
            # 输出并存储潮流计算是否成功
            is_success = result["success"]
            push!(success_status, is_success)
            println("计算 $(i+1)/10: 变压器分接头电压比 = $(current_value), 潮流计算 $(is_success ? "成功" : "失败")")
            
            append!(bus_ID, result["bus"][3, 1])
            append!(voltages, result["bus"][3, 8])
            append!(angles, result["bus"][3, 9])  # 存储相角
            append!(value, current_value )
        end
    end

    reversed_dict_bus = Dict(value => key for (key, value) in dict_bus)
    reverse_node_mapping = Dict(value => key for (key, value) in node_mapping)
    bus_ID = map(k -> reverse_node_mapping[k], bus_ID)
    bus_ID = map(k -> reversed_dict_bus[k], bus_ID)

    # 输出到Excel文件
    df = DataFrame(
        Value_KW = value,
        Bus_ID = bus_ID,
        Voltage = voltages,
        Angle = angles,  # 添加相角到数据框
        Success = success_status
    )
    
    # 输出到Excel文件（覆盖模式）
    filename = "$(elementtype)_$(field)_validation.xlsx"
    XLSX.writetable(filename, collect(eachcol(df)), names(df), overwrite=true)
    println("数据已成功导出到 $filename")

    # 读取 AC 系统数据
    ETAP_data = Dict{String, DataFrame}()
    XLSX.openxlsx(result_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            ETAP_data[sheet_name] = DataFrame(sheet[:], :auto)
        end
    end
     
    ETAP_ac_bus = ETAP_data["Sheet1"]
    ETAP_Vm = ETAP_ac_bus[2:end, 5]
    ETAP_Va = ETAP_ac_bus[2:end, 7]

    misVm = (voltages .- ETAP_Vm)./ETAP_Vm
    misVa = (angles .- ETAP_Va)./ETAP_Va  # 计算相角差异

    # 创建负载值作为x轴（从40开始，每次增加10）
    load_values = [40 + i*10 for i in 0:9]
    
    # 创建图表，显示电压幅值和相角差异
    p = Plots.plot(
        title="Voltage Comparison with ETAP",
        xlabel="Load (KVA)",
        ylabel="Relative Error",
        legend=:topright,
        size=(800, 600)
    )
    
    # 绘制电压幅值差异
    Plots.plot!(p, load_values, misVm, 
        linewidth=2, 
        marker=:circle,
        markersize=6,
        color=:blue,
        label="Voltage Magnitude"
    )
    
    # 绘制相角差异
    Plots.plot!(p, load_values, misVa, 
        linewidth=2, 
        marker=:square,
        markersize=6,
        color=:purple,
        label="Voltage Angle"
    )
    
    # 添加成功/失败状态标记
    success_indices = findall(success_status)
    failure_indices = findall(.!success_status)
    
    # 标记成功的点（幅值）
    if !isempty(success_indices)
        Plots.scatter!(p, load_values[success_indices], misVm[success_indices], 
            marker=:circle, 
            markersize=8, 
            color=:green, 
            label="Success (Magnitude)"
        )
    end
    
    # 标记失败的点（幅值）
    if !isempty(failure_indices)
        Plots.scatter!(p, load_values[failure_indices], misVm[failure_indices], 
            marker=:xcross, 
            markersize=8, 
            color=:red, 
            label="Failure (Magnitude)"
        )
    end
    
    # 标记成功的点（相角）
    if !isempty(success_indices)
        Plots.scatter!(p, load_values[success_indices], misVa[success_indices], 
            marker=:square, 
            markersize=8, 
            color=:green, 
            label="Success (Angle)"
        )
    end
    
    # 标记失败的点（相角）
    if !isempty(failure_indices)
        Plots.scatter!(p, load_values[failure_indices], misVa[failure_indices], 
            marker=:xcross, 
            markersize=8, 
            color=:red, 
            label="Failure (Angle)"
        )
    end
    
    display(p)
    
    # 输出成功率统计
    success_rate = sum(success_status) / length(success_status) * 100
    println("潮流计算成功率: $(success_rate)% ($(sum(success_status))/$(length(success_status)))")
    
    # 输出平均相对误差
    if !isempty(success_indices)
        avg_vm_error = mean(abs.(misVm[success_indices])) * 100
        avg_va_error = mean(abs.(misVa[success_indices])) * 100
        println("成功计算点的平均相对误差:")
        println("  - 电压幅值: $(round(avg_vm_error, digits=2))%")
        println("  - 电压相角: $(round(avg_va_error, digits=2))%")
    end
end
