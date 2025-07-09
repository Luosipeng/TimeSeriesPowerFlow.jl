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
        
        # 写入所有矩阵数据
        for key in keys(mpc)
            if key ∉ ["version", "baseMVA"]
                matrix = mpc[key]
                write(f, "    mpc[\"$key\"] = [\n")
                
                # 检查矩阵维度
                if ndims(matrix) == 1
                    # 处理一维数组
                    write(f, "        ")
                    for value in matrix
                        if typeof(value) <: AbstractString
                            write(f, "\"$value\" ")  # 字符串值用引号包围
                        else
                            write(f, "$(value) ")    # 数值直接写入
                        end
                    end
                    write(f, "\n")
                else
                    # 处理二维数组
                    for i in 1:size(matrix, 1)
                        write(f, "        ")
                        for j in 1:size(matrix, 2)
                            value = matrix[i,j]
                            if typeof(value) <: AbstractString
                                write(f, "\"$value\" ")  # 字符串值用引号包围
                            else
                                write(f, "$(value) ")    # 数值直接写入
                            end
                        end
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