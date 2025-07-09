"""
    compare_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String; 
                           tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3)

比较JPC格式的潮流计算结果与参考文件中的电压结果，使用JuliaPowerCase提供节点映射。
支持原始格式和包含交直流混合数据的新格式。
参数:
- results: 包含JPC格式结果的NamedTuple (通常是潮流计算的返回值)
- case: 原始的JuliaPowerCase对象，用于获取节点名称和ID映射
- reference_file: 包含参考电压值的Excel文件路径
- tolerance_mag: 电压幅值比较的容差 (默认: 1e-4)
- tolerance_ang: 电压相角比较的容差 (默认: 1e-3)

返回:
- 包含比较结果的NamedTuple，包含ac和dc两个DataFrame
"""
function compare_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String; 
                                tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3)
    # 检查参考文件是否存在
    if !isfile(reference_file)
        error("Reference file $reference_file does not exist")
    end
    
    # 读取参考文件
    reference_data = DataFrame(XLSX.readtable(reference_file, "Sheet1"))
    
    # 检查参考文件格式，确定是否包含type列（交直流混合格式）
    has_type_column = "type" in lowercase.(names(reference_data))
    
    # 创建AC节点ID到名称的映射
    ac_bus_id_to_name = Dict{Int, String}()
    ac_bus_name_to_id = Dict{String, Int}()
    for bus in case.busesAC
        if bus.in_service
            ac_bus_id_to_name[bus.bus_id] = bus.name
            ac_bus_name_to_id[bus.name] = bus.bus_id
        end
    end
    
    # 创建DC节点ID到名称的映射（如果case中有DC节点）
    dc_bus_id_to_name = Dict{Int, String}()
    dc_bus_name_to_id = Dict{String, Int}()
    if hasproperty(case, :busesDC)
        for bus in case.busesDC
            if bus.in_service
                dc_bus_id_to_name[bus.bus_id] = bus.name
                dc_bus_name_to_id[bus.name] = bus.bus_id
            end
        end
    end
    
    # 创建AC结果比较数据框
    ac_comparison_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Calc_Volt_Mag = Float64[],
        Ref_Volt_Mag = Float64[],
        Mag_Diff = Float64[],
        Mag_Error_Percent = Float64[],
        Mag_Within_Tolerance = Bool[],
        Calc_Volt_Ang = Float64[],
        Ref_Volt_Ang = Float64[],
        Ang_Diff = Float64[],
        Ang_Within_Tolerance = Bool[]
    )
    
    # 创建DC结果比较数据框
    dc_comparison_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Calc_Volt_Mag = Float64[],
        Ref_Volt_Mag = Float64[],
        Mag_Diff = Float64[],
        Mag_Error_Percent = Float64[],
        Mag_Within_Tolerance = Bool[]
    )
    
    # 获取计算结果中的电压值
    voltage_results = get_bus_voltage_results_acdc(results, case)
    ac_results = voltage_results.ac
    dc_results = voltage_results.dc
    
    # 创建AC结果中节点名称到结果的映射
    ac_name_to_result = Dict{String, NamedTuple}()
    for row in eachrow(ac_results)
        ac_name_to_result[row.Bus_Name] = (vm = row.Volt_Mag, va = row.Volt_Ang)
    end
    
    # 创建DC结果中节点名称到结果的映射
    dc_name_to_result = Dict{String, Float64}()
    for row in eachrow(dc_results)
        dc_name_to_result[row.Bus_Name] = row.Volt_Mag
    end
    
    # 遍历参考文件中的每个节点
    for i in 1:size(reference_data, 1)
        bus_name = reference_data.bus_ID[i]
        
        # 确定节点类型
        if has_type_column
            # 新格式：使用type列确定节点类型
            bus_type = uppercase(reference_data.type[i])
        else
            # 原始格式：假设所有节点都是AC
            bus_type = "AC"
        end
        
        ref_mag = reference_data.volt_mag[i]
        
        if bus_type == "AC"
            # 处理AC节点
            ref_ang = reference_data.volt_ang[i]
            
            # 在计算结果中查找对应的节点
            if haskey(ac_name_to_result, bus_name)
                result = ac_name_to_result[bus_name]
                vm = result.vm
                va = result.va
                
                # 计算差异
                mag_diff = abs(vm - ref_mag)
                mag_error_percent = ref_mag != 0 ? (mag_diff / ref_mag) * 100 : 0.0
                mag_within_tol = mag_diff <= tolerance_mag
                
                ang_diff = abs(va - ref_ang)
                ang_within_tol = ang_diff <= tolerance_ang
                
                # 获取节点ID
                bus_id = get(ac_bus_name_to_id, bus_name, -1)
                
                push!(ac_comparison_df, [
                    bus_id,
                    bus_name,
                    vm,
                    ref_mag,
                    mag_diff,
                    mag_error_percent,
                    mag_within_tol,
                    va,
                    ref_ang,
                    ang_diff,
                    ang_within_tol
                ])
            else
                @warn "AC bus $bus_name not found in calculation results, skipping comparison"
            end
        elseif bus_type == "DC"
            # 处理DC节点
            # 在计算结果中查找对应的节点
            if haskey(dc_name_to_result, bus_name)
                vm = dc_name_to_result[bus_name]
                
                # 计算差异
                mag_diff = abs(vm - ref_mag)
                mag_error_percent = ref_mag != 0 ? (mag_diff / ref_mag) * 100 : 0.0
                mag_within_tol = mag_diff <= tolerance_mag
                
                # 获取节点ID
                bus_id = get(dc_bus_name_to_id, bus_name, -1)
                
                push!(dc_comparison_df, [
                    bus_id,
                    bus_name,
                    vm,
                    ref_mag,
                    mag_diff,
                    mag_error_percent,
                    mag_within_tol
                ])
            else
                @warn "DC bus $bus_name not found in calculation results, skipping comparison"
            end
        else
            @warn "Unknown bus type '$bus_type' for bus $bus_name, skipping comparison"
        end
    end
    
    # 添加统计信息
    # AC统计
    total_ac_buses = nrow(ac_comparison_df)
    if total_ac_buses > 0
        ac_mag_match_count = count(ac_comparison_df.Mag_Within_Tolerance)
        ac_ang_match_count = count(ac_comparison_df.Ang_Within_Tolerance)
        
        println("\nAC Buses Statistics:")
        println("Total AC buses: $total_ac_buses")
        println("AC buses with voltage magnitude within tolerance: $ac_mag_match_count ($(round(ac_mag_match_count/total_ac_buses*100, digits=2))%)")
        println("AC buses with voltage angle within tolerance: $ac_ang_match_count ($(round(ac_ang_match_count/total_ac_buses*100, digits=2))%)")
        
        # 输出不匹配的AC节点
        if ac_mag_match_count < total_ac_buses
            println("\nAC buses with voltage magnitude mismatch:")
            for row in eachrow(filter(row -> !row.Mag_Within_Tolerance, ac_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Mag), Ref=$(row.Ref_Volt_Mag), Diff=$(row.Mag_Diff), Error=$(round(row.Mag_Error_Percent, digits=2))%")
            end
        end
        
        if ac_ang_match_count < total_ac_buses
            println("\nAC buses with voltage angle mismatch:")
            for row in eachrow(filter(row -> !row.Ang_Within_Tolerance, ac_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Ang), Ref=$(row.Ref_Volt_Ang), Diff=$(row.Ang_Diff)")
            end
        end
    else
        println("No AC buses found for comparison!")
    end
    
    # DC统计
    total_dc_buses = nrow(dc_comparison_df)
    if total_dc_buses > 0
        dc_mag_match_count = count(dc_comparison_df.Mag_Within_Tolerance)
        
        println("\nDC Buses Statistics:")
        println("Total DC buses: $total_dc_buses")
        println("DC buses with voltage magnitude within tolerance: $dc_mag_match_count ($(round(dc_mag_match_count/total_dc_buses*100, digits=2))%)")
        
        # 输出不匹配的DC节点
        if dc_mag_match_count < total_dc_buses
            println("\nDC buses with voltage magnitude mismatch:")
            for row in eachrow(filter(row -> !row.Mag_Within_Tolerance, dc_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Mag), Ref=$(row.Ref_Volt_Mag), Diff=$(row.Mag_Diff), Error=$(round(row.Mag_Error_Percent, digits=2))%")
            end
        end
    end
    
    return (ac = ac_comparison_df, dc = dc_comparison_df)
end

"""
    save_comparison_results(comparison_results, output_file::String)

将比较结果保存到Excel文件中。支持原始格式和交直流混合格式。
参数:
- comparison_results: 可以是DataFrame或包含ac和dc两个DataFrame的NamedTuple
- output_file: 输出文件路径
"""
function save_comparison_results(comparison_results, output_file::String)
    # 如果文件已存在，先删除它
    if isfile(output_file)
        rm(output_file)
    end
    
    # 判断结果类型
    if isa(comparison_results, DataFrame)
        # 原始格式：单个DataFrame
        XLSX.writetable(output_file, 
            Comparison = (collect(eachcol(comparison_results)), names(comparison_results))
        )
    else
        # 交直流混合格式：包含ac和dc的NamedTuple
        # 创建要写入的工作表参数
        sheets = Dict()
        
        # 添加AC结果（如果有）
        if !isempty(comparison_results.ac)
            sheets[:AC_Comparison] = (collect(eachcol(comparison_results.ac)), names(comparison_results.ac))
        end
        
        # 添加DC结果（如果有）
        if !isempty(comparison_results.dc)
            sheets[:DC_Comparison] = (collect(eachcol(comparison_results.dc)), names(comparison_results.dc))
        end
        
        # 写入Excel文件
        XLSX.writetable(output_file; sheets...)
    end
    
    println("Comparison results saved to $output_file")
end



"""
    get_bus_voltage_results(results::NamedTuple, case::JuliaPowerCase)

从JPC格式的潮流计算结果中提取节点电压信息，并与JuliaPowerCase中的节点名称对应。
参数:
- results: 包含JPC格式结果的NamedTuple (通常是潮流计算的返回值)
- case: 原始的JuliaPowerCase对象，用于获取节点名称和ID映射

返回:
- 包含节点电压结果的DataFrame
"""
function get_bus_voltage_results(results::NamedTuple, case::JuliaPowerCase)
    # 从结果中提取JPC对象
    jpc = results.value[1]
    
    # 创建节点ID到名称的映射
    bus_id_to_name = Dict{Int, String}()
    for bus in case.busesAC
        if bus.in_service
            bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # 创建结果数据框
    voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Volt_Ang = Float64[],
        Bus_Type = Int[]
    )
    
    # 遍历JPC中的所有节点
    for i in 1:size(jpc.busAC, 1)
        bus_id = Int(jpc.busAC[i, 1])
        vm = jpc.busAC[i, 8]  # 电压幅值通常在第8列
        va = jpc.busAC[i, 9]  # 电压相角通常在第9列
        bus_type = Int(jpc.busAC[i, 2])  # 节点类型
        
        # 获取节点名称
        bus_name = get(bus_id_to_name, bus_id, "Unknown_Bus_$bus_id")
        
        push!(voltage_df, [
            bus_id,
            bus_name,
            vm,
            va,
            bus_type
        ])
    end
    
    # 按节点ID排序
    sort!(voltage_df, :Bus_ID)
    
    return voltage_df
end

"""
    get_bus_voltage_results_acdc(results::NamedTuple, case::JuliaPowerCase)

从JPC格式的潮流计算结果中同时提取交流和直流节点的电压信息，并与JuliaPowerCase中的节点名称对应。
参数:
- results: 包含JPC格式结果的NamedTuple (通常是潮流计算的返回值)
- case: 原始的JuliaPowerCase对象，用于获取节点名称和ID映射

返回:
- 包含交流和直流节点电压结果的NamedTuple，包括两个DataFrame：ac和dc
"""
function get_bus_voltage_results_acdc(results::NamedTuple, case::JuliaPowerCase)
    # 创建交流节点ID到名称的映射
    ac_bus_id_to_name = Dict{Int, String}()
    for bus in case.busesAC
        if bus.in_service
            ac_bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # 创建直流节点ID到名称的映射
    dc_bus_id_to_name = Dict{Int, String}()
    for bus in case.busesDC
        if bus.in_service
            dc_bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # 创建交流结果数据框
    ac_voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Volt_Ang = Float64[],
        Bus_Type = Int[],
        Island = Int[]  # 添加岛屿标识
    )
    
    # 创建直流结果数据框
    dc_voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Bus_Type = Int[],
        Island = Int[]  # 添加岛屿标识
    )
    
    # 遍历所有岛屿
    n_islands = length(results.value)
    for island_idx in 1:n_islands
        jpc = results.value[island_idx]
        
        # 处理交流节点
        if size(jpc.busAC, 1) > 0
            # 遍历JPC中的所有交流节点
            for i in 1:size(jpc.busAC, 1)
                bus_id = Int(jpc.busAC[i, 1])
                vm = jpc.busAC[i, 8]  # 电压幅值通常在第8列
                va = jpc.busAC[i, 9]  # 电压相角通常在第9列
                bus_type = Int(jpc.busAC[i, 2])  # 节点类型
                
                # 获取节点名称
                bus_name = get(ac_bus_id_to_name, bus_id, "Unknown_AC_Bus_$bus_id")
                
                push!(ac_voltage_df, [
                    bus_id,
                    bus_name,
                    vm,
                    va,
                    bus_type,
                    island_idx  # 添加岛屿编号
                ])
            end
        end
        
        # 处理直流节点
        if size(jpc.busDC, 1) > 0
            # 遍历JPC中的所有直流节点
            for i in 1:size(jpc.busDC, 1)
                bus_id = Int(jpc.busDC[i, 1])
                vm = jpc.busDC[i, 8]  # 直流电压值通常在第8列
                bus_type = Int(jpc.busDC[i, 2])  # 节点类型
                
                # 获取节点名称
                bus_name = get(dc_bus_id_to_name, bus_id, "Unknown_DC_Bus_$bus_id")
                
                push!(dc_voltage_df, [
                    bus_id,
                    bus_name,
                    vm,
                    bus_type,
                    island_idx  # 添加岛屿编号
                ])
            end
        end
    end
    
    # 只按照Bus_ID排序
    sort!(ac_voltage_df, :Bus_ID)
    sort!(dc_voltage_df, :Bus_ID)
    
    return (ac = ac_voltage_df, dc = dc_voltage_df)
end



"""
    plot_voltage_errors(comparison_results, output_file::String="voltage_errors.png";
                       show_plot::Bool=true)

绘制电压幅值误差和相角误差的曲线图。支持原始格式和交直流混合格式。
参数:
- comparison_results: 可以是DataFrame或包含ac和dc两个DataFrame的NamedTuple
- output_file: 保存图表的文件路径或目录 (默认: "voltage_errors.png")
- show_plot: 是否显示图表 (默认: true)

返回:
- 生成的图表对象或包含图表对象的NamedTuple
"""
function plot_voltage_errors(comparison_results, output_file::String="voltage_errors.png";
                            show_plot::Bool=true)
    # 设置字体，使用系统默认字体
    default_font = "Arial"  # 使用一个通常可用的西文字体
    
    # 判断结果类型
    if isa(comparison_results, DataFrame)
        # 原始格式：单个DataFrame
        # 对节点按ID排序
        sorted_df = sort(comparison_results, :Bus_ID)
        
        # 创建一个图表布局，包含2个子图
        plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                   fontfamily=default_font)
        
        # 1. 电压幅值误差百分比曲线图
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Mag_Error_Percent, 
              label="Magnitude Error %", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Magnitude Error Percentage", 
              xlabel="Bus ID", ylabel="Error Percentage (%)")
        # 添加零误差参考线
        hline!(plt[1], [0], linestyle=:dash, color=:black, label="Zero Error")
        
        # 2. 电压相角误差曲线图
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Ang_Diff, 
              label="Angle Error", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Angle Error", 
              xlabel="Bus ID", ylabel="Error (degrees)")
        # 添加零误差参考线
        hline!(plt[2], [0], linestyle=:dash, color=:black, label="Zero Error")
        
        # 添加总标题
        plot!(plt, title="Voltage Calculation Error Analysis", titlefontsize=14)
        
        # 保存图表
        savefig(plt, output_file)
        
        if show_plot
            display(plt)
        end
        
        return plt
    else
        # 交直流混合格式：包含ac和dc的NamedTuple
        # 确保输出目录存在
        output_dir = dirname(output_file)
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        
        results = Dict()
        
        # 处理AC部分
        if !isempty(comparison_results.ac)
            # 对节点按ID排序
            sorted_ac_df = sort(comparison_results.ac, :Bus_ID)
            
            # 创建一个图表布局，包含2个子图
            ac_plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 1. 电压幅值误差百分比曲线图
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Mag_Error_Percent, 
                  label="Magnitude Error %", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Magnitude Error Percentage", 
                  xlabel="Bus ID", ylabel="Error Percentage (%)")
            # 添加零误差参考线
            hline!(ac_plt[1], [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # 2. 电压相角误差曲线图
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Ang_Diff, 
                  label="Angle Error", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Angle Error", 
                  xlabel="Bus ID", ylabel="Error (degrees)")
            # 添加零误差参考线
            hline!(ac_plt[2], [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # 添加总标题
            plot!(ac_plt, title="AC Voltage Calculation Error Analysis", titlefontsize=14)
            
            # 保存图表
            ac_error_plot_file = joinpath(dirname(output_file), "ac_voltage_errors.png")
            savefig(ac_plt, ac_error_plot_file)
            
            if show_plot
                display(ac_plt)
            end
            
            results[:ac] = (plot = ac_plt, file = ac_error_plot_file)
        end
        
        # 处理DC部分
        if !isempty(comparison_results.dc)
            # 对节点按ID排序
            sorted_dc_df = sort(comparison_results.dc, :Bus_ID)
            
            # 创建一个图表
            dc_plt = plot(size=(800, 500), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 电压幅值误差百分比曲线图
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Mag_Error_Percent, 
                  label="Magnitude Error %", marker=:circle, markersize=4, 
                  linewidth=2, title="DC Voltage Magnitude Error Percentage", 
                  xlabel="Bus ID", ylabel="Error Percentage (%)")
            # 添加零误差参考线
            hline!(dc_plt, [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # 保存图表
            dc_error_plot_file = joinpath(dirname(output_file), "dc_voltage_errors.png")
            savefig(dc_plt, dc_error_plot_file)
            
            if show_plot
                display(dc_plt)
            end
            
            results[:dc] = (plot = dc_plt, file = dc_error_plot_file)
        end
        
        return (ac = get(results, :ac, nothing), dc = get(results, :dc, nothing))
    end
end


"""
    plot_voltage_comparison(comparison_results, output_file::String="voltage_comparison.png";
                           show_plot::Bool=true)

绘制电压幅值和相角的计算值与参考值对比曲线图。支持原始格式和交直流混合格式。
参数:
- comparison_results: 可以是DataFrame或包含ac和dc两个DataFrame的NamedTuple
- output_file: 保存图表的文件路径或目录 (默认: "voltage_comparison.png")
- show_plot: 是否显示图表 (默认: true)

返回:
- 生成的图表对象或包含图表对象的NamedTuple
"""
function plot_voltage_comparison(comparison_results, output_file::String="voltage_comparison.png";
                                show_plot::Bool=true)
    # 设置字体，使用系统默认字体
    default_font = "Arial"  # 使用一个通常可用的西文字体
    
    # 判断结果类型
    if isa(comparison_results, DataFrame)
        # 原始格式：单个DataFrame
        # 对节点按ID排序
        sorted_df = sort(comparison_results, :Bus_ID)
        
        # 创建一个图表布局，包含2个子图
        plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                   fontfamily=default_font)
        
        # 1. 电压幅值比较图 (计算值 vs 参考值)
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Calc_Volt_Mag, 
              label="Calculated", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Magnitude Comparison", 
              xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Ref_Volt_Mag, 
              label="Reference", marker=:square, markersize=4, 
              linewidth=2, linestyle=:dash)
        
        # 2. 电压相角比较图 (计算值 vs 参考值)
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Calc_Volt_Ang, 
              label="Calculated", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Angle Comparison", 
              xlabel="Bus ID", ylabel="Voltage Angle (degrees)")
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Ref_Volt_Ang, 
              label="Reference", marker=:square, markersize=4, 
              linewidth=2, linestyle=:dash)
        
        # 添加总标题
        plot!(plt, title="Voltage Calculation vs Reference Comparison", titlefontsize=14)
        
        # 保存图表
        savefig(plt, output_file)
        
        if show_plot
            display(plt)
        end
        
        return plt
    else
        # 交直流混合格式：包含ac和dc的NamedTuple
        # 确保输出目录存在
        output_dir = dirname(output_file)
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        
        results = Dict()
        
        # 处理AC部分
        if !isempty(comparison_results.ac)
            # 对节点按ID排序
            sorted_ac_df = sort(comparison_results.ac, :Bus_ID)
            
            # 创建一个图表布局，包含2个子图
            ac_plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 1. 电压幅值比较图 (计算值 vs 参考值)
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Calc_Volt_Mag, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Magnitude Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Ref_Volt_Mag, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # 2. 电压相角比较图 (计算值 vs 参考值)
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Calc_Volt_Ang, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Angle Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Angle (degrees)")
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Ref_Volt_Ang, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # 添加总标题
            plot!(ac_plt, title="AC Voltage Calculation vs Reference Comparison", titlefontsize=14)
            
            # 保存图表
            ac_comparison_plot_file = joinpath(dirname(output_file), "ac_voltage_comparison.png")
            savefig(ac_plt, ac_comparison_plot_file)
            
            if show_plot
                display(ac_plt)
            end
            
            results[:ac] = (plot = ac_plt, file = ac_comparison_plot_file)
        end
        
        # 处理DC部分
        if !isempty(comparison_results.dc)
            # 对节点按ID排序
            sorted_dc_df = sort(comparison_results.dc, :Bus_ID)
            
            # 创建一个图表
            dc_plt = plot(size=(800, 500), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 电压幅值比较图 (计算值 vs 参考值)
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Calc_Volt_Mag, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="DC Voltage Magnitude Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Ref_Volt_Mag, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # 保存图表
            dc_comparison_plot_file = joinpath(dirname(output_file), "dc_voltage_comparison.png")
            savefig(dc_plt, dc_comparison_plot_file)
            
            if show_plot
                display(dc_plt)
            end
            
            results[:dc] = (plot = dc_plt, file = dc_comparison_plot_file)
        end
        
        return (ac = get(results, :ac, nothing), dc = get(results, :dc, nothing))
    end
end


"""
    analyze_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String;
                           tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3,
                           output_dir::String="./results")

分析潮流计算结果与参考文件的电压差异，生成比较报告和图表。
支持原始格式和交直流混合格式。
参数:
- results: 包含JPC格式结果的NamedTuple (通常是潮流计算的返回值)
- case: 原始的JuliaPowerCase对象，用于获取节点名称和ID映射
- reference_file: 包含参考电压值的Excel文件路径
- tolerance_mag: 电压幅值比较的容差 (默认: 1e-4)
- tolerance_ang: 电压相角比较的容差 (默认: 1e-3)
- output_dir: 输出目录 (默认: "./results")

返回:
- 包含比较结果的DataFrame或NamedTuple
"""
function analyze_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String;
                                tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3,
                                output_dir::String="./results")
    # 确保输出目录存在
    mkpath(output_dir)
    
    # 比较电压结果
    comparison_results = compare_voltage_results(results, case, reference_file, 
                                               tolerance_mag=tolerance_mag, 
                                               tolerance_ang=tolerance_ang)
    
    # 保存比较结果到Excel文件
    excel_output = joinpath(output_dir, "voltage_comparison_results.xlsx")
    save_comparison_results(comparison_results, excel_output)
    
    # 绘制电压误差曲线图
    error_plot_file = joinpath(output_dir, "voltage_errors.png")
    error_plots = plot_voltage_errors(comparison_results, error_plot_file)
    
    # 绘制电压计算值与参考值对比曲线图
    comparison_plot_file = joinpath(output_dir, "voltage_comparison.png")
    comparison_plots = plot_voltage_comparison(comparison_results, comparison_plot_file)
    
    println("\nAnalysis complete!")
    println("Comparison results saved to: $excel_output")
    
    # 根据结果类型输出不同的信息
    if isa(comparison_results, DataFrame)
        # 原始格式
        println("Voltage error plots saved to: $error_plot_file")
        println("Voltage comparison plots saved to: $comparison_plot_file")
    else
        # 交直流混合格式
        if !isempty(comparison_results.ac)
            println("AC voltage error plots saved to: $(error_plots.ac.file)")
            println("AC voltage comparison plots saved to: $(comparison_plots.ac.file)")
        end
        
        if !isempty(comparison_results.dc)
            println("DC voltage error plots saved to: $(error_plots.dc.file)")
            println("DC voltage comparison plots saved to: $(comparison_plots.dc.file)")
        end
    end
    
    return comparison_results
end