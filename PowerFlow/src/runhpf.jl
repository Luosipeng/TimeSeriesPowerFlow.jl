function runhpf(jpc, opt)
    
    # 设置迭代参数
    max_iterations = 20  # 最大迭代次数
    convergence_tolerance = 1e-4  # 收敛容差
    
    # 创建结果的深拷贝
    result_jpc = deepcopy(jpc)
    
    # 保存初始负荷数据，避免无限叠加
    initial_busAC = deepcopy(jpc.busAC)
    initial_busDC = deepcopy(jpc.busDC)
    initial_loadAC = deepcopy(jpc.loadAC)
    initial_loadDC = deepcopy(jpc.loadDC)
    
    # 保存上一次迭代的功率值，用于判断收敛
    prev_power_values = Dict()
    
    # 主迭代循环
    for iter in 1:max_iterations
        # 重置为初始负荷状态，避免累积叠加
        result_jpc.busAC = deepcopy(initial_busAC)
        result_jpc.busDC = deepcopy(initial_busDC)
        result_jpc.loadAC = deepcopy(initial_loadAC)
        result_jpc.loadDC = deepcopy(initial_loadDC)
        
        # 准备AC系统数据
        jpc1 = PowerFlow.JPC()
        jpc1.baseMVA = deepcopy(result_jpc.baseMVA)
        jpc1.busAC = deepcopy(result_jpc.busAC)
        jpc1.genAC = deepcopy(result_jpc.genAC)
        jpc1.loadAC = deepcopy(result_jpc.loadAC)
        jpc1.branchAC = deepcopy(result_jpc.branchAC)
        jpc1.pv_acsystem = deepcopy(result_jpc.pv_acsystem)
        jpc1.version = "2"
        
        # 准备DC系统数据
        jpc2 = PowerFlow.JPC()
        jpc2.baseMVA = deepcopy(result_jpc.baseMVA)
        jpc2.busDC = deepcopy(result_jpc.busDC)
        jpc2.genDC = deepcopy(result_jpc.genDC)
        jpc2.loadDC = deepcopy(result_jpc.loadDC)
        jpc2.branchDC = deepcopy(result_jpc.branchDC)
        jpc2.pv = deepcopy(result_jpc.pv)
        jpc2.version = "2"
        
        # 2. 创建当前迭代的功率值字典，用于收敛判断
        current_power_values = Dict()
        
        # 3. 根据各个换流器的模式更新相应的功率
        if !isempty(result_jpc.converter)
            # 处理工作模式为constant δs、Us (CONV_MODE==1)
            PowerFlow.update_mode_1_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # 处理工作模式为constant Ps、Qs (CONV_MODE==2)
            # 默认模式，不需要特殊处理
            
            # 处理工作模式为constant Ps、Us (CONV_MODE==3)
            PowerFlow.update_mode_3_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # 处理工作模式为constant Udc、Qs (CONV_MODE==4)
            PowerFlow.update_mode_4_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # 处理工作模式为constant Udc、Us (CONV_MODE==5)
            PowerFlow.update_mode_5_converters!(result_jpc, jpc1, jpc2, current_power_values)

            # 处理工作模式为Droop Udc、Constant Qs (CONV_MODE==6)
            PowerFlow.update_mode_6_converters!(result_jpc, jpc1, jpc2, current_power_values)

            # 处理工作模式为Droop Udc、Constant Us (CONV_MODE==7)
            PowerFlow.update_mode_7_converters!(result_jpc, jpc1, jpc2, current_power_values)

        end
        
        # 1. 独立计算交流和直流潮流（在更新负荷后进行）
        if !isempty(jpc1.busAC)
            jpc1 = runpf(jpc1, opt)
        end
        
        if !isempty(jpc2.busDC)
            jpc2 = rundcpf(jpc2, opt)
        end
        
        # 检查潮流计算是否成功
        if !isempty(jpc2.busDC)
            if !(jpc1.success && jpc2.success)
                @warn "潮流计算在第 $iter 次迭代中失败"
                break
            end
        elseif !jpc1.success
            @warn "交流潮流计算在第 $iter 次迭代中失败"
            break
        end
        
        # 4. 更新结果
        result_jpc.busAC = deepcopy(jpc1.busAC)
        result_jpc.genAC = deepcopy(jpc1.genAC)
        result_jpc.branchAC = deepcopy(jpc1.branchAC)
        # result_jpc.loadAC = deepcopy(jpc1.loadAC)  # 保持注释，因为我们要保留更新后的负荷
        result_jpc.busDC = deepcopy(jpc2.busDC)
        result_jpc.genDC = deepcopy(jpc2.genDC)
        result_jpc.branchDC = deepcopy(jpc2.branchDC)
        result_jpc.loadDC = deepcopy(jpc2.loadDC)
        result_jpc.pv = deepcopy(jpc2.pv)
        result_jpc.pv_acsystem = deepcopy(jpc1.pv_acsystem)
        result_jpc.iterationsAC = deepcopy(jpc1.iterationsAC)
        result_jpc.iterationsDC = deepcopy(jpc2.iterationsDC)
        
        # 5. 检查收敛性
        if iter > 1
            max_diff = 0.0
            for (key, value) in current_power_values
                if haskey(prev_power_values, key)
                    diff = abs(value - prev_power_values[key])
                    max_diff = max(max_diff, diff)
                end
            end
            
            if max_diff < convergence_tolerance
                @info "迭代收敛于第 $iter 次迭代，最大功率差值: $max_diff"
                break
            elseif iter == max_iterations
                @warn "达到最大迭代次数 $max_iterations,但未收敛。最大功率差值: $max_diff"
            end
        end
        
        prev_power_values = deepcopy(current_power_values)
    end
    
    # 设置求解状态
    if !isempty(result_jpc.busDC)
        if result_jpc.iterationsAC > 0 && result_jpc.iterationsDC > 0
            result_jpc.success = true
        else
            result_jpc.success = false
        end
    else
        result_jpc.success = result_jpc.iterationsAC > 0
    end
    
    return result_jpc
end

# 更新工作模式为constant δs、Us (CONV_MODE==1)的换流器
function update_mode_1_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_delta_us = filter(row -> row[CONV_MODE] == 1, eachrow(result_jpc.converter))
    if !isempty(converters_delta_us)
        converters_index = findall(row -> row[CONV_MODE] == 1, eachrow(result_jpc.converter))
        for i in eachindex(converters_delta_us[:, 1])
            conv = converters_delta_us[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            P_ac = -jpc1.genAC[gen_row, PG]
            Q_ac = -jpc1.genAC[gen_row, QG]
            
            if P_ac < 0
                P_dc = -P_ac/conv[CONV_EFF]
            else
                P_dc = -P_ac * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_Q_ac"] = Q_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_Q_AC] = Q_ac
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # 更新直流侧负荷
            update_dc_load!(result_jpc, jpc2, conv[CONV_DCBUS], P_dc)
        end
    end
end

# 更新工作模式为constant Ps、Us (CONV_MODE==3)的换流器
function update_mode_3_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_ps_us = filter(row -> row[CONV_MODE] == 3, eachrow(result_jpc.converter))
    if !isempty(converters_ps_us)
        converters_index = findall(row -> row[CONV_MODE] == 3, eachrow(result_jpc.converter))
        for i in eachindex(converters_ps_us[:, 1])
            conv = converters_ps_us[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            P_ac = -jpc1.genAC[gen_row, PG]
            Q_ac = -jpc1.genAC[gen_row, QG]
            
            if P_ac < 0
                P_dc = -P_ac/conv[CONV_EFF]
            else
                P_dc = -P_ac * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_Q_ac"] = Q_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_Q_AC] = Q_ac
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # 更新直流侧负荷
            update_dc_load!(result_jpc, jpc2, conv[CONV_DCBUS], P_dc)
        end
    end
end

# 更新工作模式为constant Udc、Qs (CONV_MODE==4)的换流器
function update_mode_4_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_udc_qs = filter(row -> row[CONV_MODE] == 4, eachrow(result_jpc.converter))
    if !isempty(converters_udc_qs)
        converters_index = findall(row -> row[CONV_MODE] == 4, eachrow(result_jpc.converter))
        for i in eachindex(converters_udc_qs[:, 1])
            conv = converters_udc_qs[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # 更新交流侧负荷
            update_ac_load!(result_jpc, jpc1, conv[CONV_ACBUS], P_ac)
        end
    end
end

# 更新工作模式为constant Udc、Us (CONV_MODE==5)的换流器
function update_mode_5_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_udc_us = filter(row -> row[CONV_MODE] == 5, eachrow(result_jpc.converter))
    if !isempty(converters_udc_us)
        converters_index = findall(row -> row[CONV_MODE] == 5, eachrow(result_jpc.converter))
        for i in eachindex(converters_udc_us[:, 1])
            conv = converters_udc_us[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # 更新交流侧发电机
            ac_gen_rows = findall(x -> x == Int(conv[CONV_ACBUS]), result_jpc.genAC[:, 1])
            ac_gen_rows_jpc1 = findall(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            if !isempty(ac_gen_rows)
                result_jpc.genAC[ac_gen_rows, PG] .= -P_ac
                jpc1.genAC[ac_gen_rows_jpc1, PG] .= -P_ac  # 更新jpc1中的交流发电机
            end
        end
    end
end

function update_mode_6_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_droop_udc_qs = filter(row -> row[CONV_MODE] == 6, eachrow(result_jpc.converter))
    if !isempty(converters_droop_udc_qs)
        converters_index = findall(row -> row[CONV_MODE] == 6, eachrow(result_jpc.converter))
        for i in eachindex(converters_droop_udc_qs[:, 1])
            conv = converters_droop_udc_qs[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            U_dc_limited = calculate_droop_voltage(-P_dc, conv[CONV_DROOP_KP],jpc2.genDC[gen_row,VG])
            jpc2.genDC[gen_row, VG] = U_dc_limited  # 更新直流侧电压
            jpc2.busDC[findfirst(jpc2.busDC[:,BUS_I].==jpc2.genDC[gen_row,GEN_BUS]), VM] = U_dc_limited  # 更新直流母线电压
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac


            # 更新交流侧负荷
            update_ac_load!(result_jpc, jpc1, conv[CONV_ACBUS], P_ac)
        end
    end

end

function update_mode_7_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_droop_udc_us = filter(row -> row[CONV_MODE] == 7, eachrow(result_jpc.converter))
    if !isempty(converters_droop_udc_us)
        converters_index = findall(row -> row[CONV_MODE] == 7, eachrow(result_jpc.converter))
        for i in eachindex(converters_droop_udc_us[:, 1])
            conv = converters_droop_udc_us[i]  # 获取第i行作为当前换流器
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            U_dc_limited = calculate_droop_voltage(-P_dc, conv[CONV_DROOP_KP],jpc2.genDC[gen_row,VG])
            jpc2.genDC[gen_row, VG] = U_dc_limited  # 更新直流侧电压
            jpc2.busDC[findfirst(jpc2.busDC[:,BUS_I].==jpc2.genDC[gen_row,GEN_BUS]), VM] = U_dc_limited  # 更新直流母线电压
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # 保存当前功率值用于收敛判断，使用行索引作为换流器ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # 更新换流器的功率
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # 更新交流侧发电机
            ac_gen_rows = findall(x -> x == Int(conv[CONV_ACBUS]), result_jpc.genAC[:, 1])
            ac_gen_rows_jpc1 = findall(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            if !isempty(ac_gen_rows)
                result_jpc.genAC[ac_gen_rows, PG] .= -P_ac
                jpc1.genAC[ac_gen_rows_jpc1, PG] .= -P_ac  # 更新jpc1中的交流发电机
            end
        end
    end
    
end

# 辅助函数：更新直流侧负荷
function update_dc_load!(result_jpc, jpc2, dc_bus_id, P_dc)
    dc_bus_id = Int(dc_bus_id)
    
    # 更新直流母线的功率
    dc_bus_rows = findall(x -> x == dc_bus_id, result_jpc.busDC[:, 1])
    dc_bus_rows_jpc2 = findall(x -> x == dc_bus_id, jpc2.busDC[:, 1])
    if !isempty(dc_bus_rows)
        result_jpc.busDC[dc_bus_rows, PD] .+= P_dc
    end
    if !isempty(dc_bus_rows_jpc2)
        jpc2.busDC[dc_bus_rows_jpc2, PD] .+= P_dc
    end
    
    # 查找对应的直流负荷
    dc_load_rows = findall(x -> x == dc_bus_id, result_jpc.loadDC[:, 2])
    dc_load_rows_jpc2 = findall(x -> x == dc_bus_id, jpc2.loadDC[:, 2])
    
    if !isempty(dc_load_rows)
        # 如果存在负荷，更新功率
        result_jpc.loadDC[dc_load_rows, LOAD_PD] .+= P_dc
    else
        # 如果没有直流负荷，则创建一个虚拟负荷
        new_load_id = isempty(result_jpc.loadDC) ? 1 : maximum(result_jpc.loadDC[:, 1]) + 1
        
        # 创建新的负荷行，初始化所有值为0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # 负荷ID
        new_load_row[LOAD_CND] = dc_bus_id  # 母线ID
        new_load_row[LOAD_STATUS] = 1       # 工作状态
        new_load_row[LOAD_PD] = P_dc        # 有功功率
        new_load_row[LOADP_PERCENT] = 1.0   # 有功百分比
        
        # 将新的负荷行添加到result_jpc的loadDC中
        result_jpc.loadDC = vcat(result_jpc.loadDC, reshape(new_load_row, 1, :))
    end
    if !isempty(dc_load_rows_jpc2)
        # 如果存在负荷，更新功率
        jpc2.loadDC[dc_load_rows_jpc2, LOAD_PD] .+= P_dc
    else
        # 如果没有直流负荷，则创建一个虚拟负荷
        new_load_id = isempty(jpc2.loadDC) ? 1 : maximum(jpc2.loadDC[:, 1]) + 1
        
        # 创建新的负荷行，初始化所有值为0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # 负荷ID
        new_load_row[LOAD_CND] = dc_bus_id  # 母线ID
        new_load_row[LOAD_STATUS] = 1       # 工作状态
        new_load_row[LOAD_PD] = P_dc        # 有功功率
        new_load_row[LOADP_PERCENT] = 1.0   # 有功百分比
        
        # 将新的负荷行添加到jpc2的loadDC中
        jpc2.loadDC = vcat(jpc2.loadDC, reshape(new_load_row, 1, :))
    end
end

# 辅助函数：更新交流侧负荷
function update_ac_load!(result_jpc, jpc1, ac_bus_id, P_ac)
    ac_bus_id = Int(ac_bus_id)
    
    # 更新交流母线的功率
    ac_bus_rows = findall(x -> x == ac_bus_id, result_jpc.busAC[:, 1])
    ac_bus_rows_jpc1 = findall(x -> x == ac_bus_id, jpc1.busAC[:, 1])
    if !isempty(ac_bus_rows)
        result_jpc.busAC[ac_bus_rows, PD] .+= P_ac
    end
    if !isempty(ac_bus_rows_jpc1)
        jpc1.busAC[ac_bus_rows_jpc1, PD] .+= P_ac
    end
    
    # 查找对应的交流负荷
    ac_load_rows = findall(x -> x == ac_bus_id, result_jpc.loadAC[:, 2])
    ac_load_rows_jpc1 = findall(x -> x == ac_bus_id, jpc1.loadAC[:, 2])
    
    if !isempty(ac_load_rows)
        # 如果存在负荷，更新功率
        result_jpc.loadAC[ac_load_rows, LOAD_PD] .+= P_ac
    else
        # 如果没有交流负荷，则创建一个虚拟负荷
        new_load_id = isempty(result_jpc.loadAC) ? 1 : maximum(result_jpc.loadAC[:, 1]) + 1
        
        # 创建新的负荷行，初始化所有值为0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # 负荷ID
        new_load_row[LOAD_CND] = ac_bus_id  # 母线ID
        new_load_row[LOAD_STATUS] = 1       # 工作状态
        new_load_row[LOAD_PD] = P_ac        # 有功功率
        new_load_row[LOADP_PERCENT] = 1.0   # 有功百分比
        
        # 将新的负荷行添加到result_jpc的loadAC中
        result_jpc.loadAC = vcat(result_jpc.loadAC, reshape(new_load_row, 1, :))
    end
    if !isempty(ac_load_rows_jpc1)
        # 如果存在负荷，更新功率
        jpc1.loadAC[ac_load_rows_jpc1, LOAD_PD] .+= P_ac
    else
        # 如果没有交流负荷，则创建一个虚拟负荷
        new_load_id = isempty(jpc1.loadAC) ? 1 : maximum(jpc1.loadAC[:, 1]) + 1
        
        # 创建新的负荷行，初始化所有值为0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # 负荷ID
        new_load_row[LOAD_CND] = ac_bus_id  # 母线ID
        new_load_row[LOAD_STATUS] = 1       # 工作状态
        new_load_row[LOAD_PD] = P_ac        # 有功功率
        new_load_row[LOADP_PERCENT] = 1.0   # 有功百分比
        
        # 将新的负荷行添加到jpc1的loadAC中
        jpc1.loadAC = vcat(jpc1.loadAC, reshape(new_load_row, 1, :))
    end
end


# 下垂控制辅助函数
function calculate_droop_voltage(P_dc, k_p,U_dc_ref=1.0, U_dc_min=0.95, U_dc_max=1.05)
    """
    计算下垂控制的直流电压
    U_dc = U_dc_ref - k_p * P_dc
    """
    U_dc_droop = U_dc_ref - k_p * P_dc
    # 电压限幅
    U_dc_limited = max(U_dc_min, min(U_dc_max, U_dc_droop))
    return U_dc_limited
end