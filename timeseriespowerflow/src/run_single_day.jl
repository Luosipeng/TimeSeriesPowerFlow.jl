function run_single_day(old_jpc, opt, day_load_matrix, day_price_line, day_irradiance_line)

    jpc = deepcopy(old_jpc)
    converter_ac_bus = Int.(jpc.converter[:, CONV_ACBUS])
    converter_ac_p_mw = jpc.converter[:, CONV_P_AC]
    converter_ac_q_mvar = jpc.converter[:, CONV_Q_AC]
    converter_dc_bus = Int.(jpc.converter[:, CONV_DCBUS])
    converter_dc_p_mw = jpc.converter[:, CONV_P_DC]

    # 创建母线ID到索引的映射
    ac_bus_map = Dict(jpc.busAC[i, BUS_I] => i for i in 1:size(jpc.busAC, 1))
    dc_bus_map = Dict(jpc.busDC[i, BUS_I] => i for i in 1:size(jpc.busDC, 1))

    # 交流侧处理 - 减去换流器的功率
    for i in 1:length(converter_ac_bus)
        ac_bus_id = converter_ac_bus[i]
        if haskey(ac_bus_map, ac_bus_id)
            idx = ac_bus_map[ac_bus_id]
            if jpc.converter[i, CONV_MODE] != 6 && jpc.converter[i, CONV_MODE] != 7
            jpc.busAC[idx, PD] -= converter_ac_p_mw[i]
            end
            if jpc.converter[i, CONV_MODE] != 7
                jpc.busAC[idx, QD] -= converter_ac_q_mvar[i]
            end
        else
            @warn "找不到ID为 $ac_bus_id 的交流母线"
        end
    end

    # 直流侧处理 - 减去换流器的功率
    for i in 1:length(converter_dc_bus)
        dc_bus_id = converter_dc_bus[i]
        if haskey(dc_bus_map, dc_bus_id)
            idx = dc_bus_map[dc_bus_id]
            if jpc.converter[i, CONV_MODE] != 6 && jpc.converter[i, CONV_MODE] != 7
                jpc.busDC[idx, PD] -= converter_dc_p_mw[i]
            end
        else
            @warn "找不到ID为 $dc_bus_id 的直流母线"
        end
    end

    # 删除或者减去对应的负荷对应的负荷
    for i in 1:size(jpc.loadAC, 1)
        load_bus = jpc.loadAC[i, LOAD_CND]
        for j in 1:length(converter_ac_bus)
            if load_bus == converter_ac_bus[j]
                # 只减去换流器的功率部分，而不是删除整个负荷
                if jpc.converter[j, CONV_MODE] != 6 && jpc.converter[j, CONV_MODE] != 7
                    jpc.loadAC[i, LOAD_PD] -= converter_ac_p_mw[j]
                end
                if jpc.converter[j, CONV_MODE] != 7
                    # 只减去换流器的无功功率部分，而不是删除整个负荷
                    jpc.loadAC[i, LOAD_QD] -= converter_ac_q_mvar[j]
                end
            end
        end
    end

    for i in 1:size(jpc.loadDC, 1)
        load_bus = jpc.loadDC[i, LOAD_CND]
        for j in 1:length(converter_dc_bus)
            if load_bus == converter_dc_bus[j]
                if jpc.converter[j, CONV_MODE] != 6 && jpc.converter[j, CONV_MODE] != 7
                    # 只减去换流器的功率部分，而不是删除整个负荷
                    jpc.loadDC[i, LOAD_PD] -= converter_dc_p_mw[j]
                end
            end
        end
    end

   new_jpc, ac_node_mapping, dc_node_mapping = TimeSeriesPowerFlow.renumber_hybrid_system(jpc)
    

    ac_node_reverse_mapping = Dict{Int, Float64}()
    dc_node_reverse_mapping = Dict{Int, Float64}()
    
    # 填充反向映射字典
    for (old_num, new_num) in ac_node_mapping
        ac_node_reverse_mapping[new_num] = old_num
    end
    
    for (old_num, new_num) in dc_node_mapping
        dc_node_reverse_mapping[new_num] = old_num
    end

    # 使用新的JPC对象进行后续操作
    # 提取负荷信息
    loadAC = new_jpc.loadAC
    loadDC = new_jpc.loadDC

    inserviced_line_AC = findall(loadAC[:, LOAD_STATUS] .== 1)
    inserviced_line_DC = findall(loadDC[:, LOAD_STATUS] .== 1)

    # 过滤掉未使用的负荷
    loadAC = loadAC[inserviced_line_AC, :]
    loadDC = loadDC[inserviced_line_DC, :]

    num_load_ac = size(loadAC, 1)
    num_load_dc = size(loadDC, 1)

    # 提取负荷数据
    loadAC_PD = zeros(num_load_ac, 24)
    loadAC_QD = zeros(num_load_ac, 24)
    # loadDC_PD = zeros(num_load_dc*24, 2)
    
    for i in 1:24
        # 找出当前小时的所有数据
        idx = findall(day_load_matrix[:, 1] .== i)
        # 为AC负荷赋值
        loadAC_PD[:, i] = day_load_matrix[idx, 3]
        loadAC_QD[:, i] = day_load_matrix[idx, 4]
        
    end
    
    loadDC_PD = zeros(num_load_dc, 24)
    for i in 1:24
        # 找出当前小时的所有数据
        # idx = findall(day_load_matrix[:, 1] .== i)
        
        # 为DC负荷赋值
        loadDC_PD[:, i] = loadDC[:, LOAD_PD]
    end
    # TODO: 这里的day_load_matrix需要根据实际情况调整
    # loadDC_PD = loadDC[:, LOAD_PD]

    # 提取负荷位置信息 - 使用新编号
    ac_bus_ids = loadAC[:, LOAD_CND]
    dc_bus_ids = loadDC[:, LOAD_CND]

    # 构建负荷连接矩阵 - 使用新编号的总节点数
    total_buses = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)
    
    # AC负荷连接矩阵
    Cld_ac = zeros(total_buses, size(loadAC, 1))
    for i in eachindex(ac_bus_ids)
        bus_id = Int(ac_bus_ids[i])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cld_ac[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的交流母线"
        end
    end
    
    # DC负荷连接矩阵
    Cld_dc = zeros(total_buses, size(loadDC, 1))
    for i in eachindex(dc_bus_ids)
        bus_id = Int(dc_bus_ids[i])
        if bus_id in new_jpc.busDC[:, BUS_I]
            Cld_dc[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的直流母线"
        end
    end
    
    # 提取发电机信息 - 使用新JPC
    genAC = new_jpc.genAC
    genDC = new_jpc.genDC
    genAC = genAC[genAC[:, GEN_STATUS] .== 1, :]  # 仅保留状态为1的发电机
    genDC = genDC[genDC[:, GEN_STATUS] .== 1, :]  # 仅保留状态为1的发电机

    # 找出平衡节点 - 使用新编号
    # slack_busAC = new_jpc.busAC[findfirst(new_jpc.busAC[:, BUS_TYPE] .== REF), BUS_I]
    # keep_indices = genAC[:, GEN_BUS] .!= slack_busAC
    # genAC = genAC[keep_indices, :]

    # 找出换流器节点 - 使用新编号
    if !isempty(new_jpc.converter)
        conv_ac_bus = new_jpc.converter[:, CONV_ACBUS]
        keep_indices = genAC[:, GEN_BUS] .!= conv_ac_bus
        genAC = genAC[keep_indices, :]
    end
    
    # 提取发电机功率
    genAC_PG = genAC[:, PG]
    # 构建发电机连接矩阵 - 使用新编号
    Cgen_ac = zeros(total_buses, size(genAC, 1))
    for i in eachindex(genAC[:, GEN_BUS])
        bus_id = Int(genAC[i, GEN_BUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cgen_ac[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的交流母线"
        end
    end

    # 提取换流器信息 - 使用新JPC
    converters = new_jpc.converter
    Pij_inv = zeros(size(converters, 1))
    Pij_rec = zeros(size(converters, 1))
    Qij_inv = zeros(size(converters, 1))
    η_rec = converters[:,CONV_EFF]
    η_inv = converters[:, CONV_EFF]
    
    for i in eachindex(converters[:, 1])
        Pac = converters[i, CONV_P_AC]
        Pdc = converters[i, CONV_P_DC]
        Qac = converters[i, CONV_Q_AC]

        if Pac <= 0
            Pij_inv[i] = -Pac
        else
            Pij_inv[i] = 0
        end

        if Pdc <= 0
            Pij_rec[i] = -Pdc
        else
            Pij_rec[i] = 0
        end

        if Qac < 0
            Qij_inv[i] = -Qac
        else
            Qij_inv[i] = 0
        end
    end

    # 构建换流器连接矩阵 - 使用新编号
    Cconv_ac = zeros(total_buses, size(converters, 1))
    for i in eachindex(converters[:, CONV_ACBUS])
        bus_id = Int(converters[i, CONV_ACBUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cconv_ac[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的交流母线"
        end
    end

    Cconv_dc = zeros(total_buses, size(converters, 1))
    for i in eachindex(converters[:, CONV_DCBUS])
        bus_id = Int(converters[i, CONV_DCBUS])
        if bus_id in new_jpc.busDC[:, BUS_I]

            Cconv_dc[ bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的直流母线"
        end
    end

    # 提取光伏AC系统信息 - 使用新JPC
    pv_ac = new_jpc.pv_acsystem

    # 提取光伏AC系统的基准光照强度
    pv_ac_base_irradiance = pv_ac[:, PV_AC_IRRADIANCE]
    # 提取PV光伏出力上限
    pv_ac_p_mw_ratio = zeros(size(pv_ac, 1), 24)
    pv_ac_p_mw_ratio = inv.(pv_ac_base_irradiance) * day_irradiance_line[:, 2]'
    pv_ac_p_mw = pv_ac[:, PV_AC_INVERTER_PAC] 

    pv_ac_inverter_pac = pv_ac[:, PV_AC_INVERTER_PAC]

    # 提取光伏AC连接矩阵 - 使用新编号
    Cpv_ac = zeros(total_buses, size(pv_ac, 1))
    for i in eachindex(pv_ac[:, PV_AC_BUS])
        bus_id = Int(pv_ac[i, PV_AC_BUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cpv_ac[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的交流母线"
        end
    end

    # 提取光伏信息 - 使用新JPC
    pvarray = new_jpc.pv

    # 提取光伏的基准光照强度
    pv_base_irradiance = pvarray[:, PV_IRRADIANCE]

    # 提取PV光伏出力上限
    # pv_max_p_mw_ratio = day_irradiance_line[:, 2] ./ pv_base_irradiance
    pv_max_p_mw_ratio = inv.(pv_base_irradiance)* day_irradiance_line[:, 2]'
    pv_max_p_mw = pvarray[:, PV_VMPP].* pvarray[:, PV_IMPP]./1000000

    pv_isc = pvarray[:, PV_ISC]
    pv_impp = pvarray[:, PV_IMPP]
    
    # 提取光伏连接矩阵 - 使用新编号
    Cpv_dc = zeros(total_buses, size(pvarray, 1))
    for i in eachindex(pvarray[:, PV_BUS])
        bus_id = Int(pvarray[i, PV_BUS])
        if bus_id in new_jpc.busDC[:, BUS_I]
            Cpv_dc[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的直流母线"
        end
    end

    # 提取储能信息 - 使用新JPC
    storage = new_jpc.storage
    # 提取储能出力信息
    ess_power_capacity_mw = storage[:, ESS_POWER_CAPACITY]
    ess_energy_capacity_mwh = storage[:, ESS_ENERGY_CAPACITY]
    ess_initial_soc = storage[:, ESS_SOC_INIT]

    # 提取储能SOC信息
    ess_min_soc = storage[:, ESS_SOC_MIN]
    ess_max_soc = storage[:, ESS_SOC_MAX]
    ess_efficiency = storage[:, ESS_EFFICIENCY]

    # 提取储能连接矩阵 - 使用新编号
    # 注意：原代码中这里的维度似乎有问题，应该是total_buses而不是只有busDC的行数
    Cstorage_ac = zeros(total_buses, size(storage, 1))
    for i in eachindex(storage[:, ESS_BUS])
        bus_id = Int(storage[i, ESS_BUS])
        if bus_id in new_jpc.busDC[:, BUS_I]  # 这里可能需要检查是否应该是busAC
            Cstorage_ac[bus_id, i] = 1
        else
            @warn "找不到ID为 $bus_id 的直流母线"
        end
    end

    result = run_dynamic_dispatch(new_jpc,
        Cld_ac, Cld_dc, 
        loadAC_PD, loadAC_QD,
        loadDC_PD,
        genAC_PG,
        Cgen_ac, 
        Cconv_ac, Cconv_dc, 
        η_rec, η_inv,
        Cpv_ac,
        Cpv_dc,
        pv_ac_p_mw_ratio,
        pv_ac_p_mw,
        pv_max_p_mw,
        pv_max_p_mw_ratio,
        Cstorage_ac,
        ess_initial_soc,
        ess_max_soc,
        ess_min_soc,
        ess_power_capacity_mw,
        ess_energy_capacity_mwh,
        ess_efficiency,
        day_price_line
        )

    # 全天潮流计算
    results_pf = Array{Any}(undef, 24)
    # 赋值
    for t in 1:24
        # 负荷幅值
        loadAC_PD_t = loadAC_PD[:, t]
        loadAC_QD_t = loadAC_QD[:, t]
        loadDC_PD_t = loadDC_PD[:, t]

        # 更新负荷
        for i in eachindex(jpc.loadAC[:, LOAD_CND])
            bus_id = jpc.loadAC[i, LOAD_CND]
            if bus_id in jpc.busAC[:, BUS_I]
                bus_index = findfirst(x -> x == bus_id, jpc.busAC[:, BUS_I])
                jpc.busAC[bus_index, PD] = loadAC_PD_t[i]
                jpc.busAC[bus_index, QD] = loadAC_QD_t[i]

                # 更新负荷
                idx = findfirst(jpc.loadAC[:, LOAD_CND] .== bus_id)
                if idx !== nothing
                    jpc.loadAC[idx, LOAD_PD] = loadAC_PD_t[i]
                    jpc.loadAC[idx, LOAD_QD] = loadAC_QD_t[i]
                else
                    @warn "找不到ID为 $bus_id 的交流母线负荷"
                end
            else
                @warn "找不到ID为 $bus_id 的交流母线"
            end
        end

        for i in eachindex(jpc.loadDC[:, LOAD_CND])
            bus_id = jpc.loadDC[i, LOAD_CND]
            if bus_id in jpc.busDC[:, BUS_I]
                bus_index = findfirst(x -> x == bus_id, jpc.busDC[:, BUS_I])
                jpc.busDC[bus_index, PD] = loadDC_PD_t[i]

                # 更新负荷
                idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_id)
                if idx !== nothing
                    jpc.loadDC[idx, LOAD_PD] = loadDC_PD_t[i]
                else
                    @warn "找不到ID为 $bus_id 的直流母线负荷"
                end
            else
                @warn "找不到ID为 $bus_id 的直流母线"
            end
        end

        # 将结果重新赋值到jpc.converter中
        for i in eachindex(converters[:, 1])
            if result["Pij_inv"][i,t] > 0
                jpc.converter[i, CONV_P_AC] = result["Pij_inv"][i,t]
            else
                jpc.converter[i, CONV_P_AC] = -result["Pij_rec"][i,t] * jpc.converter[i, CONV_EFF]
            end

            if result["Pij_rec"][i,t] > 0
                jpc.converter[i, CONV_P_DC] = result["Pij_rec"][i,t]
            else
                jpc.converter[i, CONV_P_DC] = -result["Pij_inv"][i,t] * jpc.converter[i, CONV_EFF]
            end
            # jpc.converter[i, CONV_Q_AC] = Qij_inv_result[i]
        end

        # 根据控制模式修改负荷或者发电机功率
        for i in eachindex(converters[:, 1])
            if jpc.converter[i, CONV_MODE] == 6
                # 如果是Droop_Udc_Qs模式，直接将交流功率设置为0
                genDC_idx = findfirst(jpc.genDC[:, GEN_BUS] .== jpc.converter[i, CONV_DCBUS])
                jpc.genDC[genDC_idx, PG] = -jpc.converter[i, CONV_P_DC]
                # 第一步：找到对应的索引
                idx = findfirst(jpc.busAC[:, BUS_I] .== jpc.converter[i, CONV_ACBUS])
                
                # 第二步：如果找到了索引，更新该母线的有功功率
                # if idx !== nothing
                #     jpc.busAC[idx, PD] += jpc.converter[i, CONV_P_AC]
                # end
                # # 第一步：找到对应的索引
                # idx = findfirst(jpc.loadAC[:, LOAD_CND] .== jpc.converter[i, CONV_ACBUS])

                # # 第二步：如果找到了索引，更新该负荷的有功功率
                # if idx !== nothing
                #     jpc.loadAC[idx, LOAD_PD] += jpc.converter[i, CONV_P_AC]
                # end

            elseif jpc.converter[i, CONV_MODE] == 7
                # 如果是Droop_Udc_Us模式，直接将无功功率设置为0
                genAC_idx = findfirst(jpc.genAC[:, GEN_BUS] .== jpc.converter[i, CONV_ACBUS])
                jpc.genAC[genAC_idx, PG] = -jpc.converter[i, CONV_P_AC]

                genDC_idx = findfirst(jpc.genDC[:, GEN_BUS] .== jpc.converter[i, CONV_DCBUS])
                jpc.genDC[genDC_idx, PG] = -jpc.converter[i, CONV_P_DC]

            end
        end

        #为储能负荷进行赋值
        for i in eachindex(storage[:, ESS_BUS])
            bus_id = storage[i, ESS_BUS]
            if bus_id in new_jpc.busDC[:, BUS_I]
                bus_id = dc_node_reverse_mapping[bus_id]  # 使用新编号的直流母线ID
                # ac_bus_count = size(new_jpc.busAC, 1)
                line = findfirst(jpc.busDC[:,BUS_I] .== bus_id)
                jpc.busDC[line, PD] = result["ess_charge"][i,t]- result["ess_discharge"][i,t]

                # 更新储能负荷
                idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_id)
                if idx !== nothing
                    jpc.loadDC[idx, LOAD_PD] = result["ess_charge"][i,t]- result["ess_discharge"][i,t]
                else
                    @warn "找不到ID为 $bus_id 的直流母线负荷"
                end
            else
                @warn "找不到ID为 $bus_id 的直流母线"
            end
        end

        # 根据光照强度对光伏进行赋值
        jpc.pv[:, PV_ISC] = pv_isc.*pv_max_p_mw_ratio[:,t]
        jpc.pv[:, PV_IMPP] = pv_impp.*pv_max_p_mw_ratio[:,t]

        # 更新光伏AC系统的功率
        jpc.pv_acsystem[:, PV_AC_INVERTER_PAC] = pv_ac_inverter_pac.*pv_ac_p_mw_ratio[:, t]

        results_pf[t] = PowerFlow.runhpf(jpc, opt)
    end

    return results_pf
end