# function runautopf(jpc::JPC ,opt, mode::String="lindistflow")
#     # 回收处理为纯功率负荷的换流器
#     converter_ac_bus = Int.(jpc.converter[:,1])
#     converter_dc_bus = Int.(jpc.converter[:,2])

#     converter_ac_p_mw = jpc.converter[:,4]
#     converter_ac_q_mvar = jpc.converter[:,5]
#     converter_dc_p_mw = jpc.converter[:,6]

#     # 创建母线ID到索引的映射
#     ac_bus_map = Dict(jpc.busAC[i, 1] => i for i in 1:size(jpc.busAC, 1))
#     dc_bus_map = Dict(jpc.busDC[i, 1] => i for i in 1:size(jpc.busDC, 1))

#     # 交流侧处理
#     for i in 1:length(converter_ac_bus)
#         ac_bus_id = converter_ac_bus[i]
#         if haskey(ac_bus_map, ac_bus_id)
#             idx = ac_bus_map[ac_bus_id]
#             jpc.busAC[idx, PD] -= converter_ac_p_mw[i]
#             jpc.busAC[idx, QD] -= converter_ac_q_mvar[i]
#         else
#             @warn "找不到ID为 $ac_bus_id 的交流母线"
#         end
#     end

#     # 直流侧处理
#     for i in 1:length(converter_dc_bus)
#         dc_bus_id = converter_dc_bus[i]
#         if haskey(dc_bus_map, dc_bus_id)
#             idx = dc_bus_map[dc_bus_id]
#             jpc.busDC[idx, PD] -= converter_dc_p_mw[i]
#         else
#             @warn "找不到ID为 $dc_bus_id 的直流母线"
#         end
#     end

#     # 删除对应的负荷
#     load_ac_index = findall(x -> x in converter_ac_bus, jpc.loadAC[:,2])
#     load_dc_index = findall(x -> x in converter_dc_bus, jpc.loadDC[:,2])

#     if !isempty(load_ac_index)
#         keep_ac_rows = ones(Bool, size(jpc.loadAC, 1))
#         keep_ac_rows[load_ac_index] .= false
#         jpc.loadAC = jpc.loadAC[keep_ac_rows, :]
#     end

#     if !isempty(load_dc_index)
#         keep_dc_rows = ones(Bool, size(jpc.loadDC, 1))
#         keep_dc_rows[load_dc_index] .= false
#         jpc.loadDC = jpc.loadDC[keep_dc_rows, :]
#     end
#     # 使用最大功率传输的方式进行功率自动分配
#     # 1. 计算每个换流器的功率传输能力
#     if mode == "mppt"
#         pdc_mw = PowerFlow.runmppt(jpc)
#     elseif mode == "lindistflow"
#         pac_mw, qac_mvar, pdc_mw = PowerFlow.runlindistflow(jpc)
#     else
#         error("不支持的模式: $mode")
#     end
#     if mode == "mppt"
#         # 2. 计算换流器交流测功率
#         efficiency = jpc.converter[:, 7]
#         for i in 1:length(pdc_mw)
#             if pdc_mw[i] > 0
#                 pac_mw = - pdc_mw[i] * efficiency[i]
#             else
#                 pac_mw = -pdc_mw[i] / efficiency[i]
#             end
#             jpc.converter[i, 4] = pac_mw
#             jpc.converter[i, 6] = pdc_mw[i]
#         end
#     else
#         # 2. 更新换流器的功率
#         for i in 1:length(pdc_mw)
#             jpc.converter[i, 4] = pac_mw[i]  # 更新交流功率
#             jpc.converter[i, 5] = qac_mvar[i]  # 更新无功功率
#             jpc.converter[i, 6] = pdc_mw[i]  # 更新直流功率
#         end
#     end

#     # 3. 重新将换流器处理为负荷
#     # 交流侧处理 - 将换流器的功率重新添加为负荷
#     for i in 1:length(converter_ac_bus)
#         ac_bus_id = converter_ac_bus[i]
#         if haskey(ac_bus_map, ac_bus_id)
#             idx = ac_bus_map[ac_bus_id]
#             jpc.busAC[idx, PD] += jpc.converter[i, 4]  # 添加有功功率
#             jpc.busAC[idx, QD] += jpc.converter[i, 5]  # 添加无功功率
#         else
#             @warn "找不到ID为 $ac_bus_id 的交流母线"
#         end
#     end

#     # 直流侧处理 - 将换流器的功率重新添加为负荷
#     for i in 1:length(converter_dc_bus)
#         dc_bus_id = converter_dc_bus[i]
#         if haskey(dc_bus_map, dc_bus_id)
#             idx = dc_bus_map[dc_bus_id]
#             jpc.busDC[idx, PD] += jpc.converter[i, 6]  # 添加直流功率
#         else
#             @warn "找不到ID为 $dc_bus_id 的直流母线"
#         end
#     end

#     # 可选：重新创建负荷记录（如果需要在loadAC和loadDC中也体现换流器负荷）
#     for i in 1:length(converter_ac_bus)
#         # 创建新的交流负荷记录
#         new_ac_load = zeros(size(jpc.loadAC, 2))
#         new_ac_load[1] = size(jpc.loadAC, 1) + i  # 负荷ID
#         new_ac_load[2] = converter_ac_bus[i]      # 母线ID
#         new_ac_load[3] = 1                        # 状态
#         new_ac_load[4] = jpc.converter[i, 4]      # 有功功率
#         new_ac_load[5] = jpc.converter[i, 5]      # 无功功率
#         new_ac_load[8] = 1                        # 纯功率负荷
        
#         # 将新负荷添加到loadAC
#         jpc.loadAC = vcat(jpc.loadAC, reshape(new_ac_load, 1, :))
#     end

#     for i in 1:length(converter_dc_bus)
#         # 创建新的直流负荷记录
#         new_dc_load = zeros(size(jpc.loadDC, 2))
#         new_dc_load[1] = size(jpc.loadDC, 1) + i  # 负荷ID
#         new_dc_load[2] = converter_dc_bus[i]      # 母线ID
#         new_dc_load[3] = 1                        # 状态
#         new_dc_load[4] = jpc.converter[i, 6]      # 直流功率
#         new_dc_load[8] = 1                        # 纯功率负荷
        
#         # 将新负荷添加到loadDC
#         jpc.loadDC = vcat(jpc.loadDC, reshape(new_dc_load, 1, :))
#     end


#     results = PowerFlow.runhpf(jpc, opt)
    
#     return results
# end

function runautopf(jpc::Utils.JPC, opt)
    # 剥离换流器的功率负荷
    # 交流侧处理
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

   new_jpc = PowerFlow.renumber_hybrid_system(jpc)
    
    # 使用新的JPC对象进行后续操作
    # 提取负荷信息
    loadAC = new_jpc.loadAC
    loadDC = new_jpc.loadDC
    loadAC = loadAC[loadAC[:, LOAD_STATUS] .== 1, :]  # 仅保留状态为1的负荷
    loadDC = loadDC[loadDC[:, LOAD_STATUS] .== 1, :]  # 仅保留状态为1的负荷

    # 提取负荷数据
    loadAC_PD = loadAC[:, LOAD_PD]
    loadAC_QD = loadAC[:, LOAD_QD]
    loadDC_PD = loadDC[:, LOAD_PD]

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
    slack_busAC = new_jpc.busAC[findfirst(new_jpc.busAC[:, BUS_TYPE] .== REF), BUS_I]
    keep_indices = genAC[:, GEN_BUS] .!= slack_busAC
    genAC = genAC[keep_indices, :]

    # 找出换流器节点 - 使用新编号
    conv_ac_bus = new_jpc.converter[:, CONV_ACBUS]
    keep_indices = genAC[:, GEN_BUS] .!= conv_ac_bus
    genAC = genAC[keep_indices, :]
    
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

    # 提取光伏信息 - 使用新JPC
    pvarray = new_jpc.pv
    # 提取PV光伏出力上限
    pv_max_p_mw = pvarray[:, PV_VMPP].* pvarray[:, PV_IMPP]./1000000
    
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
    storage_p_mw = storage[:, ESS_P_MW]
    storage_max_p_mw = storage[:, ESS_MAX_P_MW]
    storage_min_p_mw = storage[:, ESS_MIN_P_MW]

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

   Pij_inv_result, Qij_inv_result, Pij_rec_result, Pess_mw_result = PowerFlow.runlindistflow(
    new_jpc, opt, 
    Cld_ac, Cld_dc, 
    loadAC_PD, loadAC_QD,
    loadDC_PD,
    genAC_PG,
    Cgen_ac, 
    Cconv_ac, Cconv_dc, 
    Pij_inv, Pij_rec, Qij_inv, 
    η_rec, η_inv,
    Cpv_dc,
    pv_max_p_mw,
    Cstorage_ac,
    storage_p_mw,
    storage_max_p_mw,
    storage_min_p_mw
)

    # 将结果重新赋值到jpc.converter中
    for i in eachindex(converters[:, 1])
        if Pij_inv_result[i] > 0
            jpc.converter[i, CONV_P_AC] = Pij_inv_result[i]
        else
            jpc.converter[i, CONV_P_AC] = -Pij_rec_result[i] * jpc.converter[i, CONV_EFF]
        end

        if Pij_rec_result[i] > 0
            jpc.converter[i, CONV_P_DC] = Pij_rec_result[i]
        else
            jpc.converter[i, CONV_P_DC] = -Pij_inv_result[i] * jpc.converter[i, CONV_EFF]
        end
        jpc.converter[i, CONV_Q_AC] = Qij_inv_result[i]
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
            if idx !== nothing
                jpc.busAC[idx, QD] += jpc.converter[i, CONV_Q_AC]
            end
            # 第一步：找到对应的索引
            idx = findfirst(jpc.loadAC[:, LOAD_CND] .== jpc.converter[i, CONV_ACBUS])

            # 第二步：如果找到了索引，更新该负荷的有功功率
            if idx !== nothing
                jpc.loadAC[idx, LOAD_QD] += jpc.converter[i, CONV_Q_AC]
            end

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
            bus_index = findfirst(x -> x == bus_id, new_jpc.busDC[:, BUS_I])
            # ac_bus_count = size(new_jpc.busAC, 1)
            jpc.busDC[bus_index, PD] += Pess_mw_result[i]

            # 更新储能负荷
            idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_index)
            if idx !== nothing
                jpc.loadDC[idx, LOAD_PD] += Pess_mw_result[i]
            else
                @warn "找不到ID为 $bus_id 的直流母线负荷"
            end
        else
            @warn "找不到ID为 $bus_id 的直流母线"
        end
    end

    # 运行潮流计算
    results = PowerFlow.runhpf_iteration(jpc, opt)
    return results
end

