"""
读取负荷数据文件
参数:
    file_path: 负荷数据Excel文件的路径
返回:
    time_column: 时间列
    time_str_column: 时间字符串列
    load_names: 负荷名称列表
    load_data: 包含所有数据的DataFrame
"""
function read_load_data(file_path)
    # 读取Excel文件到DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # 提取时间列和时间字符串列
    time_column = data[!, 1]  # 第一列为时间
    time_str_column = data[!, 2]  # 第二列为时间字符串
    
    # 获取负荷名称（列名）
    load_names = names(data)[3:end]
    
    return time_column, time_str_column, load_names, data
end

function read_price_data(file_path)
    # 读取Excel文件到DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # 提取时间列和时间字符串列
    time_column = data[!, 1]  # 第一列为时间
    time_str_column = data[!, 2]  # 第二列为时间字符串
    
    
    return time_column, time_str_column, data
end

function read_irradiance_data(file_path)
    # 读取Excel文件到DataFrame
    data = DataFrame(XLSX.readtable(file_path, 1, header=true))
    
    # 提取小时列
    hour_column = data[!, 1]  # 第一列为小时
    # 提取时间字符串列
    time_str_column = data[!, 2]  # 第二列为时间字符串
    
    return hour_column, time_str_column, data
end

function create_time_series_loads(case::Utils.JuliaPowerCase, data, load_names, num_days)
    # 获取时间点数量
    num_timepoints = size(data, 1)
    hours_per_day = 24
    
    # 检查时间点是否能被24整除
    if num_timepoints % hours_per_day != 0
        @warn "时间点数量 ($num_timepoints) 不是24的整数倍，将只处理完整的天数"
        num_days = div(num_timepoints, hours_per_day)  # 调整为完整的天数
        num_timepoints = num_days * hours_per_day
    else
        num_days = div(num_timepoints, hours_per_day)
    end
    
    # 过滤出投运的负荷
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # 如果没有投运的负荷，返回空结果
    if isempty(in_service_loads)
        @warn "没有投运的负荷"
        return Dict{Int, Array{Float64, 2}}()
    end
    
    # 创建一个字典，用于存储数据表中每个负荷名称对应的列索引
    load_column_indices = Dict{String, Int}()
    
    # 获取数据表的列名作为字符串
    data_column_names = string.(names(data))
    
    # 查找每个负荷名称在数据表中的列索引
    for name in load_names
        if name in data_column_names
            col_idx = findfirst(x -> string(x) == name, names(data))
            if col_idx !== nothing
                load_column_indices[name] = col_idx
            end
        end
    end
    
    # 创建一个字典，按天存储负荷时间序列数据
    # 格式: day_loads[day] = [timepoint, bus_id, p_mw, q_mvar, const_z, const_i, const_p]
    day_loads = Dict{Int, Array{Float64, 2}}()
    
    # 获取系统中所有唯一的母线ID
    unique_bus_ids = unique([load.bus for load in in_service_loads])
    num_buses = length(unique_bus_ids)
    
    # 创建母线ID到索引的映射
    bus_id_to_index = Dict(bus_id => i for (i, bus_id) in enumerate(unique_bus_ids))
    
    # 按天处理
    for day in 1:num_days
        @info "处理第 $day 天的负荷数据，共 $num_days 天"
        
        # 计算当天的时间点范围
        start_t = (day - 1) * hours_per_day + 1
        end_t = day * hours_per_day
        
        # 创建当天的负荷矩阵
        # 每个时间点每个母线一行，列为 [时间点索引, 母线ID, P, Q, const_z, const_i, const_p]
        day_load_matrix = zeros(num_buses * hours_per_day, 7)
        
        # 处理当天的每个时间点
        for hour in 1:hours_per_day
            t = start_t + hour - 1  # 全局时间点索引
            
            # 创建一个字典，用于累加连接到同一母线的负荷
            bus_load_sum = Dict{Int, Vector{Float64}}()
            
            # 初始化母线负荷特性累加器
            # [p_mw, q_mvar, p_mw*const_z, q_mvar*const_z, p_mw*const_i, q_mvar*const_i, p_mw*const_p, q_mvar*const_p]
            for bus_id in unique_bus_ids
                bus_load_sum[bus_id] = zeros(8)
            end
            
            # 为每个负荷创建一个字典，用于累加匹配的负荷值
            load_power_dict = Dict{String, Float64}()
            
            # 从数据表中获取指定时间点的负荷值
            for (name, col_idx) in load_column_indices
                load_power_dict[name] = data[t, col_idx]
            end
            
            # 处理每个负荷
            for load in in_service_loads
                load_name = load.name
                bus_id = load.bus
                
                # 查找所有与当前负荷名称匹配的负荷
                matching_loads = [name for name in keys(load_power_dict) if startswith(name, load_name)]
                
                # 计算实际的有功和无功负荷
                if !isempty(matching_loads)
                    # 找到匹配的负荷，累加所有匹配负荷的视在功率值
                    total_apparent_power = sum(load_power_dict[name] for name in matching_loads)
                    
                    # 获取原始功率值
                    old_p = load.p_mw
                    old_q = load.q_mvar
                    
                    # 计算原始功率因数（保持不变）
                    if old_p != 0 || old_q != 0
                        # 计算原始视在功率和功率因数
                        old_s = sqrt(old_p^2 + old_q^2)
                        if old_s > 0
                            pf = old_p / old_s  # 功率因数 = P/S
                            # 保持功率因数不变，根据新的视在功率计算有功和无功功率
                            new_s = total_apparent_power / 1000.0  # 假设数据单位为kVA，转换为MVA
                            new_p = new_s * pf  # P = S * cos(φ)
                            new_q = new_s * sqrt(1 - pf^2)  # Q = S * sin(φ)
                            
                            # 保持无功功率的符号（感性或容性）
                            if old_q < 0
                                new_q = -new_q
                            end
                        else
                            # 原始视在功率为0的情况，假设功率因数为0.9
                            new_s = total_apparent_power / 1000.0
                            new_p = new_s * 0.9
                            new_q = new_s * sqrt(1 - 0.9^2)
                        end
                    else
                        # 原始功率都为0的情况，假设功率因数为0.9
                        new_s = total_apparent_power / 1000.0
                        new_p = new_s * 0.9
                        new_q = new_s * sqrt(1 - 0.9^2)
                    end
                else
                    # 不在数据表中的负荷保持原值不变
                    new_p = load.p_mw * load.scaling
                    new_q = load.q_mvar * load.scaling
                end
                
                # 获取负荷特性
                const_z = load.const_z_percent / 100.0
                const_i = load.const_i_percent / 100.0
                const_p = load.const_p_percent / 100.0
                
                # 累加到对应母线的负荷
                bus_load_sum[bus_id][1] += new_p  # 总有功功率
                bus_load_sum[bus_id][2] += new_q  # 总无功功率
                bus_load_sum[bus_id][3] += new_p * const_z  # 恒阻抗有功
                bus_load_sum[bus_id][4] += new_q * const_z  # 恒阻抗无功
                bus_load_sum[bus_id][5] += new_p * const_i  # 恒电流有功
                bus_load_sum[bus_id][6] += new_q * const_i  # 恒电流无功
                bus_load_sum[bus_id][7] += new_p * const_p  # 恒功率有功
                bus_load_sum[bus_id][8] += new_q * const_p  # 恒功率无功
            end
            
            # 填充负荷矩阵，合并相同母线的负荷
            for bus_id in unique_bus_ids
                idx = bus_id_to_index[bus_id]
                row_idx = (hour - 1) * num_buses + idx
                
                # 获取累加的负荷数据
                load_data = bus_load_sum[bus_id]
                total_p = load_data[1]
                total_q = load_data[2]
                
                # 计算合并后的负荷特性百分比
                const_z_percent = total_p > 0 ? load_data[3] / total_p : 0.0
                const_i_percent = total_p > 0 ? load_data[5] / total_p : 0.0
                const_p_percent = total_p > 0 ? load_data[7] / total_p : 0.0
                
                # 确保百分比之和为1
                sum_percent = const_z_percent + const_i_percent + const_p_percent
                if sum_percent > 0
                    const_z_percent /= sum_percent
                    const_i_percent /= sum_percent
                    const_p_percent /= sum_percent
                else
                    # 默认为恒功率负荷
                    const_p_percent = 1.0
                    const_z_percent = 0.0
                    const_i_percent = 0.0
                end
                
                # 填充负荷矩阵的对应行
                day_load_matrix[row_idx, :] = [
                    hour,           # 时间点索引（小时）
                    bus_id,         # 母线ID
                    total_p,        # 有功负荷(MW)
                    total_q,        # 无功负荷(MVAr)
                    const_z_percent,  # 恒阻抗负荷百分比
                    const_i_percent,  # 恒电流负荷百分比
                    const_p_percent   # 恒功率负荷百分比
                ]
            end
        end
        
        # 将当天的负荷数据存储到字典中
        day_loads[day] = day_load_matrix
    end
    
    return day_loads
end

function create_time_series_prices(price_profiles, num_days=nothing)
    # 获取天数
    total_days = size(price_profiles, 1)
    hours_per_day = 24
    
    # 如果指定了天数，确保不超过可用的天数
    if num_days === nothing
        num_days = total_days
    else
        num_days = min(num_days, total_days)
    end
    
    @info "处理电价数据，共 $num_days 天"
    
    # 创建一个字典，按天存储电价时间序列数据
    # 格式: day_prices[day] = [hour, price]
    day_prices = Dict{Int, Array{Float64, 2}}()
    
    # 按天处理
    for day in 1:num_days
        @info "处理第 $day 天的电价数据，共 $num_days 天"
        
        # 创建当天的电价矩阵 [hour, price]
        day_price_matrix = zeros(hours_per_day, 2)
        
        # 处理当天的每个小时
        for hour in 1:hours_per_day
            # 列索引（第一列是日期，所以从第二列开始）
            col_idx = hour + 1
            
            # 获取电价
            price = price_profiles[day, col_idx]
            
            # 填充电价矩阵
            day_price_matrix[hour, :] = [hour, price]
        end
        
        # 将当天的电价数据存储到字典中
        day_prices[day] = day_price_matrix
    end
    
    return day_prices
end

function create_time_series_irradiance(irradiance_profiles, num_days=nothing)
    # 获取天数
    total_days = size(irradiance_profiles, 1)
    hours_per_day = 24
    
    # 如果指定了天数，确保不超过可用的天数
    if num_days === nothing
        num_days = total_days
    else
        num_days = min(num_days, total_days)
    end
    
    @info "处理光照强度数据，共 $num_days 天"
    
    # 创建一个字典，按天存储光照强度时间序列数据
    # 格式: day_irradiance[day] = [hour, irradiance]
    day_irradiance = Dict{Int, Array{Float64, 2}}()
    
    # 按天处理
    for day in 1:num_days
        @info "处理第 $day 天的光照强度数据，共 $num_days 天"
        
        # 创建当天的光照强度矩阵 [hour, irradiance]
        day_irradiance_matrix = zeros(hours_per_day, 2)
        
        # 处理当天的每个小时
        for hour in 1:hours_per_day
            # 列索引（第一列是日期，所以从第二列开始）
            col_idx = hour + 1
            
            # 获取光照强度
            irradiance = irradiance_profiles[day, col_idx]
            
            # 填充光照强度矩阵
            day_irradiance_matrix[hour, :] = [hour, irradiance]
        end
        
        # 将当天的光照强度数据存储到字典中
        day_irradiance[day] = day_irradiance_matrix
    end
    
    return day_irradiance
end


"""
处理所有时间点的负荷数据，并进行潮流计算
参数:
    case: powerflow系统案例
    data: 负荷数据DataFrame
    load_names: 负荷名称列表
    price_profiles: 价格数据
    opt: 潮流计算选项
返回:
    results: 所有时间点的潮流计算结果
"""
function runtdpf(case, data, load_names, price_profiles, irradiance_profiles, opt)
    # 获取时间点数量和天数
    num_timepoints = size(data, 1)
    hours_per_day = 24
    num_days = div(num_timepoints, hours_per_day)
    
    # 检查时间点是否能被24整除
    if num_timepoints % hours_per_day != 0
        @warn "时间点数量 ($num_timepoints) 不是24的整数倍，将只处理完整的天数"
        num_timepoints = num_days * hours_per_day  # 调整为完整的天数
    end
    
    # 生成按天组织的负荷数据（已合并相同母线的负荷）
    day_loads = TimeSeriesPowerFlow.create_time_series_loads(case, data, load_names, num_days)
    day_price = TimeSeriesPowerFlow.create_time_series_prices(price_profiles, num_days)
    day_irradiance = TimeSeriesPowerFlow.create_time_series_irradiance(irradiance_profiles, num_days)

    # 将powerflow案例转换为JPC格式
    jpc = TimeSeriesPowerFlow.JuliaPowerCase2Jpc(case)

    # 提取电网孤岛
    jpc_list, isolated = TimeSeriesPowerFlow.extract_islands_acdc(jpc)
    n_islands = length(jpc_list)

    # 创建结果数组
    results = Array{Any}(undef, n_islands, num_days, 24)
    
    # 使用并行处理每一天
    @threads for d in 1:num_days
        @info "处理第 $d 天，共 $num_days 天"
        
        # 获取当天的负荷数据
        if !haskey(day_loads, d)
            @error "没有第 $d 天的负荷数据"
            continue
        end
        
        day_load_matrix = day_loads[d]
        day_price_line = day_price[d]
        day_irradiance_line = day_irradiance[d]
            
        # 提取当天的负荷矩阵
        load_matrix_list, isolated_load_matrix = TimeSeriesPowerFlow.extract_load_matrix_by_islands(day_load_matrix, jpc_list)
        
        # 对每个孤岛执行潮流计算
        # 注意：这里不再使用@threads，因为外层已经并行处理每一天
        for i in 1:n_islands
            results[i,d,:] = run_single_day(jpc_list[i], opt, load_matrix_list[i], day_price_line, day_irradiance_line)
        end

        @info "完成第 $d 天的处理"
    end
    
    return results
end


function extract_load_matrix_by_islands(day_load_matrix, jpc_list)
    # 初始化结果列表，用于存储每个岛屿对应的负荷矩阵
    load_matrix_list = []
    
    # 假设day_load_matrix的第一列是时间，第二列是母线编号
    # 从day_load_matrix中提取母线编号
    bus_ids_in_load = day_load_matrix[:, 2]
    
    # 对每个岛屿进行处理
    for i in eachindex(jpc_list)
        jpck = jpc_list[i]
        
        # 获取当前岛屿中的AC母线编号
        island_ac_buses = jpck.busAC[:, BUS_I]
        
        # 找出day_load_matrix中属于当前岛屿的行
        load_indices = findall(bus_id -> bus_id in island_ac_buses, bus_ids_in_load)
        
        if !isempty(load_indices)
            # 提取对应的负荷矩阵
            island_load_matrix = day_load_matrix[load_indices, :]
            push!(load_matrix_list, island_load_matrix)
        else
            # 如果当前岛屿没有对应的负荷，则添加一个空矩阵
            push!(load_matrix_list, similar(day_load_matrix, 0, size(day_load_matrix, 2)))
        end
    end
    
    # 处理孤立的母线（不在任何有效岛屿中的母线）
    isolated_load_indices = findall(bus_id -> !any(bus_id in jpck.busAC[:, BUS_I] for jpck in jpc_list), bus_ids_in_load)
    isolated_load_matrix = isempty(isolated_load_indices) ? similar(day_load_matrix, 0, size(day_load_matrix, 2)) : day_load_matrix[isolated_load_indices, :]
    
    return load_matrix_list, isolated_load_matrix
end


