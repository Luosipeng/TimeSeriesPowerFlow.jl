module demo
    # using Gtk

    # 定义 ROOT_DIR 为当前项目的根目录
    const ROOT_DIR = dirname(dirname(@__FILE__))

    if Sys.iswindows()
        push!(LOAD_PATH, joinpath(ROOT_DIR, "src"))
    else
        push!(LOAD_PATH, joinpath(ROOT_DIR, "src"))
    end
    
    using PowerFlow
    using DataFrames
    using XLSX

    # 修改清屏函数，增加跨平台支持
    function clear_screen()
        if Sys.iswindows()
            try
                run(`cmd /c cls`)
            catch
                println("\n" ^ 50)  # 备选方案：打印多个换行
            end
        else
            try
                run(`clear`)
            catch
                println("\n" ^ 50)
            end
        end
    end
    
    function parse_matlab_case_file(filepath)
        # 读取文件内容
        content = read(filepath, String)
        
        # 创建空字典存储结果
        mpc = Dict{String, Any}()
        
        # 解析 baseMVA
        if occursin(r"mpc\.baseMVA\s*=\s*(\d+)", content)
            basemva_match = match(r"mpc\.baseMVA\s*=\s*(\d+)", content)
            mpc["baseMVA"] = parse(Float64, basemva_match[1])
        end
        
        # 解析 version
        if occursin(r"mpc\.version\s*=\s*'(\d+)'", content)
            version_match = match(r"mpc\.version\s*=\s*'(\d+)'", content)
            mpc["version"] = version_match[1]
        end
        
        # 解析矩阵或字符串数据的函数
        function extract_data(content, key)
            # 分割内容为行
            lines = split(content, '\n')
            
            # 找到矩阵开始的行
            start_pattern = "mpc.$key = ["
            end_pattern = "];"
            
            start_idx = 0
            end_idx = 0
            
            # 查找矩阵的开始和结束位置
            for (i, line) in enumerate(lines)
                if occursin(start_pattern, line)
                    start_idx = i
                elseif start_idx > 0 && occursin(end_pattern, line)
                    end_idx = i
                    break
                end
            end
            
            # 如果找到了矩阵
            if start_idx > 0 && end_idx > 0
                # 提取矩阵内容
                matrix_lines = lines[start_idx+1:end_idx-1]
                
                # 过滤掉空行和注释行
                matrix_lines = filter(line -> !isempty(strip(line)) && !startswith(strip(line), '%'), matrix_lines)
                
                # 检查是否包含字符串数据
                contains_strings = any(line -> occursin("'", line) || occursin("\"", line), matrix_lines)
                
                if contains_strings
                    # 处理字符串数据
                    matrix = []
                    for line in matrix_lines
                        # 移除行尾的分号和注释
                        line = replace(line, r";.*$" => "")
                        # 提取引号中的内容和数字
                        parts = String[]
                        current_str = ""
                        in_quotes = false
                        quote_char = nothing
                        
                        for char in line
                            if char in ['\'', '"']
                                if !in_quotes
                                    in_quotes = true
                                    quote_char = char
                                elseif char == quote_char
                                    in_quotes = false
                                    if !isempty(current_str)
                                        push!(parts, current_str)
                                        current_str = ""
                                    end
                                else
                                    current_str *= char
                                end
                            elseif in_quotes
                                current_str *= char
                            elseif !isspace(char) && char != ';'
                                current_str *= char
                            elseif !isempty(current_str)
                                push!(parts, current_str)
                                current_str = ""
                            end
                        end
                        
                        if !isempty(current_str)
                            push!(parts, current_str)
                        end
                        
                        # 过滤掉空字符串
                        parts = filter(!isempty, parts)
                        
                        if !isempty(parts)
                            push!(matrix, parts)
                        end
                    end
                    return length(matrix) > 0 ? reduce(vcat, transpose.(matrix)) : nothing
                else
                    # 处理数值数据
                    matrix = []
                    for line in matrix_lines
                        # 移除行尾的分号和注释
                        line = replace(line, r";.*$" => "")
                        # 分割并转换为数值
                        try
                            row = parse.(Float64, split(strip(line)))
                            if !isempty(row)
                                push!(matrix, row)
                            end
                        catch
                            @warn "无法解析行: $line"
                            continue
                        end
                    end
                    return length(matrix) > 0 ? reduce(vcat, transpose.(matrix)) : nothing
                end
            end
            return nothing
        end
        
        # 查找所有可能的矩阵名称
        matrix_names = String[]
        for line in split(content, '\n')
            m = match(r"mpc\.(\w+)\s*=\s*\[", line)
            if m !== nothing
                push!(matrix_names, m[1])
            end
        end
        
        # 解析每个找到的矩阵
        for name in matrix_names
            if name ∉ ["version", "baseMVA"]  # 跳过已处理的特殊字段
                matrix = extract_data(content, name)
                if matrix !== nothing
                    mpc[name] = matrix
                end
            end
        end
        
        return mpc
    end
    
    function save_to_julia_file(mpc, output_filepath)
        open(output_filepath, "w") do f
            write(f, """function case_data()
        mpc = Dict{String, Any}()
    
    """)
            
            # 写入version
            if haskey(mpc, "version")
                write(f, "    mpc[\"version\"] = \"$(mpc["version"])\"\n\n")
            end
            
            # 写入baseMVA
            if haskey(mpc, "baseMVA")
                write(f, "    mpc[\"baseMVA\"] = $(mpc["baseMVA"])\n\n")
            end
            
            for (key, matrix) in mpc
                if key ∉ ["version", "baseMVA"]
                    write(f, "    mpc[\"$key\"] = [\n")
                    for idx in CartesianIndices(matrix)
                        i, j = Tuple(idx)
                        if j == 1  # 每行开始
                            write(f, "        ")
                        end
                        value = matrix[i,j]
                        if typeof(value) <: AbstractString
                            write(f, "\"$value\" ")
                        else
                            write(f, "$(value) ")
                        end
                        if j == size(matrix, 2)  # 每行结束
                            write(f, ";\n")
                        end
                    end
                    write(f, "    ]\n\n")
                end
            end

            write(f, "    return mpc\nend")
        end
    end
    
    function convert_matpower_case(input_filepath, output_filepath)
        try
            println("正在解析MATLAB文件...")
            mpc = parse_matlab_case_file(input_filepath)
            
            println("正在保存为Julia文件...")
            save_to_julia_file(mpc, output_filepath)
            
            println("转换完成！")
            println("输入文件：$input_filepath")
            println("输出文件：$output_filepath")
            
            # 打印数据统计
            println("\n数据统计：")
            for (key, value) in mpc
                if value isa Array
                    println("$key 矩阵大小: $(size(value))")
                else
                    println("$key: $(value)")
                end
            end
            
            return mpc
        catch e
            println("转换过程中出现错误：")
            println(e)
            return nothing
        end
    end

    function Poweranalysis(mpc)
        # 直接在当前环境中执行文件内容
        # eval(Meta.parse(read(input_filepath, String)))
        
        # 调用 case_data 函数获取数据
        # mpc = Main.case_data()
        
        # 设置潮流计算选项
        opt = PowerFlow.options()
        opt["PF"]["NR_ALG"] = "bicgstab"
        opt["PF"]["ENFORCE_Q_LIMS"] = 0
        
        println("\n开始执行潮流计算...")
        result = PowerFlow.runpf(mpc, opt)
        
        return result
    end


    # 添加一个用于显示结果的辅助函数
    function display_results(result)
        if result === nothing
            println("计算失败或没有结果")
            return
        end
        
        try
            # 创建输出文件名
            output_file = "power_flow_results.xlsx"
            
            # 如果文件已存在，先删除它
            isfile(output_file) && rm(output_file)
            
            # 提取所需数据
            bus_data = result["bus"]
            gen_data = result["gen"]
            branch_data = result["branch"]
            baseMVA = result["baseMVA"]
            
            # 定义基础列名（标准列）
            base_bus_headers = [
                "bus_i", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va",
                "baseKV", "zone", "Vmax", "Vmin"
            ]
            
            base_gen_headers = [
                "bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status",
                "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min", "Qc1max", "Qc2min",
                "Qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf"
            ]
            
            base_branch_headers = [
                "fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC",
                "ratio", "angle", "status", "angmin", "angmax", "Pf", "Qf", "Pt", "Qt"
            ]
            
            # 辅助函数：根据实际列数生成列名
            function generate_headers(base_headers, actual_cols)
                if actual_cols <= length(base_headers)
                    return base_headers[1:actual_cols]
                else
                    extra_cols = ["Col_$(i)" for i in (length(base_headers)+1):actual_cols]
                    return vcat(base_headers, extra_cols)
                end
            end
            
            # 获取实际列数
            bus_cols = size(bus_data, 2)
            gen_cols = size(gen_data, 2)
            branch_cols = size(branch_data, 2)
            
            # 生成适配的列名
            bus_headers = generate_headers(base_bus_headers, bus_cols)
            gen_headers = generate_headers(base_gen_headers, gen_cols)
            branch_headers = generate_headers(base_branch_headers, branch_cols)
            
            # 创建DataFrames
            bus_df = DataFrame(bus_data, bus_headers)
            gen_df = DataFrame(gen_data, gen_headers)
            branch_df = DataFrame(branch_data, branch_headers)
            baseMVA_df = DataFrame(baseMVA = [baseMVA])
            
            # 创建一个新的XLSX文件
            XLSX.openxlsx(output_file, mode="w") do xf
                # 重命名第一个工作表为Bus_Data并写入数据
                XLSX.rename!(XLSX.getsheet(xf, 1), "Bus_Data")
                XLSX.writetable!(XLSX.getsheet(xf, "Bus_Data"), bus_df)
                
                # 添加其他工作表并写入数据
                gen_sheet = XLSX.addsheet!(xf, "Gen_Data")
                branch_sheet = XLSX.addsheet!(xf, "Branch_Data")
                base_sheet = XLSX.addsheet!(xf, "BaseMVA")
                
                XLSX.writetable!(gen_sheet, gen_df)
                XLSX.writetable!(branch_sheet, branch_df)
                XLSX.writetable!(base_sheet, baseMVA_df)
            end
            
            println("\n结果已保存到文件: $output_file")
            println("包含以下表格：")
            println("Bus_Data：母线数据")
            println("Gen_Data：发电机数据")
            println("Branch_Data：支路数据")
            println("BaseMVA：基准功率")
            
            # 打印详细的数据维度信息
            println("\n数据维度信息：")
            println("Bus data: $(size(bus_data)) ($(bus_cols) 列)")
            println("Gen data: $(size(gen_data)) ($(gen_cols) 列)")
            println("Branch data: $(size(branch_data)) ($(branch_cols) 列)")
            println("基准功率: $(baseMVA) MVA")
            
            # 打印每个数据表实际使用的列数
            println("\n实际列数信息：")
            println("Bus 标准列数: $(length(base_bus_headers)), 实际列数: $bus_cols")
            println("Gen 标准列数: $(length(base_gen_headers)), 实际列数: $gen_cols")
            println("Branch 标准列数: $(length(base_branch_headers)), 实际列数: $branch_cols")
            
        catch e
            println("保存结果时出现错误：")
            println(e)
            println("\n错误堆栈跟踪：")
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
        end
    end
    
    function julia_main()
        while true
            clear_screen()
            
            println("\n请输入MATPOWER算例文件的完整路径（输入'e'退出）：")
            println("(例如: C:/Users/DELL/Desktop/matpower8.0/data/case9.m)")
            input_filepath = strip(readline())
            
            if lowercase(input_filepath) == "e"
                println("\n程序已退出。")
                break
            end
            
            if !isfile(input_filepath)
                println("\n错误：文件不存在！")
                println("路径：$input_filepath")
                println("\n按回车键继续...")
                readline()
                continue
            end
            
            mpc = convert_matpower_case(input_filepath, "case_data.jl")
            # 使用实际输入的文件路径，而不是硬编码的路径
            try
                # 执行潮流分析并计时
                time_result = @timed Poweranalysis(mpc)
                result = time_result.value
                
                # 显示性能信息
                println("\n计算耗时: $(round(time_result.time, digits=3)) 秒")
                println("内存分配: $(round(time_result.bytes/1024/1024, digits=2)) MB")
                
                # 显示结果
                println("潮流分析完成！")
                display_results(result)
            catch e
                println("\n执行过程中出现错误：")
                println(e)
                println("\n错误堆栈跟踪：")
                for (exc, bt) in Base.catch_stack()
                    showerror(stdout, exc, bt)
                    println()
                end
            end
            
            println("\n按回车键继续下一次计算...")
            readline()
        end
        return nothing
    end
end
