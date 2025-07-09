function dc_element_validate(elementtype, field, Dict_busdc, result_path, opt, mpc; connectedname=nothing, Fbus=nothing, Tbus=nothing)
    bus_ID = []
    voltages = []
    value = []
    success_status = [] # 存储成功状态
    
    if elementtype == "DCLUMPLOAD" && field == "RatedKW"
        if connectedname !== nothing
            aim_bus = findall(i -> mpc["bus"][i, 1] == connectedname, 1:size(mpc["bus"], 1))
            aim_load = findall(i -> mpc["load"][i, 2] == connectedname, 1:size(mpc["load"], 1))
        end
        start_value = 0.04
        for i in 0:9
            current_value = start_value + i * 0.01
            mpc["load"][aim_load, 4] .= current_value
            mpc["bus"][aim_bus, 3] .= current_value
            result = PowerFlow.rundcpf(mpc, opt)
            
            # 输出并存储潮流计算是否成功
            is_success = result["success"]
            push!(success_status, is_success)
            println("计算 $(i+1)/10: 负荷值 = $(current_value*1000) KW, 潮流计算 $(is_success ? "成功" : "失败")")
            
            append!(bus_ID, result["bus"][aim_bus, 1])
            append!(voltages, result["bus"][aim_bus, 8])
            append!(value, current_value * 1000)
        end
    end
    
    reversed_Dict_busdc = Dict(value => key for (key, value) in Dict_busdc)
    bus_ID = map(k -> reversed_Dict_busdc[k], bus_ID)
    
    # 输出到Excel文件
    df = DataFrame(
        Value_KW = value,
        Bus_ID = bus_ID,
        Voltage = voltages,
        Success = success_status
    )
    
    # 输出到Excel文件（覆盖模式）
    filename = "$(elementtype)_$(field)_validation.xlsx"
    XLSX.writetable(filename, collect(eachcol(df)), names(df), overwrite=true)
    println("数据已成功导出到 $filename")

    # 读取 DC 系统数据
    ETAP_data = Dict{String, DataFrame}()
    XLSX.openxlsx(result_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            ETAP_data[sheet_name] = DataFrame(sheet[:], :auto)
        end
    end
    
    ETAP_dc_bus = ETAP_data["Sheet1"]
    ETAP_bus_ID = ETAP_dc_bus[2:end, 2]
    ETAP_voltages = ETAP_dc_bus[2:end, 3]

    # 计算相对误差
    misV = (voltages .- ETAP_voltages)./ETAP_voltages
    
    # 创建负载值作为x轴（从40开始，每次增加10）
    load_values = [40 + i*10 for i in 0:9]
    
    # 绘制图表，显示电压差异
    p = Plots.plot(
        title="DC Voltage Comparison with ETAP",
        xlabel="Load (KW)",
        ylabel="Relative Error",
        legend=:topright,
        size=(800, 600)
    )
    
    # 绘制电压差异
    Plots.plot!(p, load_values, misV, 
        linewidth=2, 
        marker=:circle,
        markersize=6,
        color=:blue,
        label="Voltage"
    )
    
    # 添加成功/失败状态标记
    success_indices = findall(success_status)
    failure_indices = findall(.!success_status)
    
    # 标记成功的点
    if !isempty(success_indices)
        Plots.scatter!(p, load_values[success_indices], misV[success_indices], 
            marker=:circle, 
            markersize=8, 
            color=:green, 
            label="Success"
        )
    end
    
    # 标记失败的点
    if !isempty(failure_indices)
        Plots.scatter!(p, load_values[failure_indices], misV[failure_indices], 
            marker=:xcross, 
            markersize=8, 
            color=:red, 
            label="Failure"
        )
    end
    
    display(p)
    
    # 输出成功率统计
    success_rate = sum(success_status) / length(success_status) * 100
    println("潮流计算成功率: $(success_rate)% ($(sum(success_status))/$(length(success_status)))")
    
    # 输出平均相对误差
    if !isempty(success_indices)
        avg_v_error = mean(abs.(misV[success_indices])) * 100
        println("成功计算点的平均相对误差:")
        println("  - 电压: $(round(avg_v_error, digits=2))%")
    end
    
    return bus_ID, voltages, value, misV, success_status
end
