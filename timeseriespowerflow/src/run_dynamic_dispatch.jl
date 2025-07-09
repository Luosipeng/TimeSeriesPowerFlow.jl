using  JuMP
using  Ipopt
using  Gurobi

function run_dynamic_dispatch(new_jpc,
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
        day_price_line,
        num_hours=24)

    # 提取连接信息
    branchAC = new_jpc.branchAC
    branchDC = new_jpc.branchDC
    converter = new_jpc.converter

    rac = branchAC[:,BR_R]
    rdc = branchDC[:,BR_R]
    rconv = zeros(size(converter, 1))
    # 提取母线数量
    n_nodes = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)
    A , branch_data = TimeSeriesPowerFlow.build_incidence_matrix_td(n_nodes, branchAC, branchDC, converter)

     # 初始化
    nbr = size(A, 1)  # 支路数
    nc = size(converter, 1)  # 变流器数
    ng = size(genAC_PG, 1)  # 交流发电机数
    npv = size(new_jpc.pv, 1)  # 光伏设备数
    ns = size(new_jpc.storage, 1)  # 储能设备数

    # 搜寻支路索引
    # dcbranch = filter(x -> x[3] == 2, branch_data)
    dcbranch_indices = findall(x -> x[3] == 2, branch_data)
    converter_indices = findall(x -> x[3] == 3, branch_data)
    acbranch_indices = findall(x -> x[3] == 1, branch_data)

    r = zeros(nbr)  # 初始化支路电阻
    r[acbranch_indices] = rac  # 交流支路电阻
    r[dcbranch_indices] = rdc  # 直流支路电阻
    r[converter_indices] = rconv  # 变流器电阻

    # 优化模型配置
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)  # 输出详细程度
    set_optimizer_attribute(model, "max_iter", 3000)  # 最大迭代次数
    set_optimizer_attribute(model, "tol", 1e-6)       # 收敛容差

    # model = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(model, "OutputFlag", 0)  # 显示求解过程
    # set_optimizer_attribute(model, "TimeLimit", 3600)  # 最大求解时间(秒)
    # set_optimizer_attribute(model, "MIPGap", 1e-4)    # MIP相对间隙


    # 定义变量
    @variable(model, Pij[1:nbr,1:num_hours])  # 支路功率流
    @variable(model, Pgen[1:ng,1:num_hours])  # 逆变器功率流
    @variable(model, Pij_inv[1:nc,1:num_hours] >= 0)  # 变流器功率流
    @variable(model, Pij_rec[1:nc,1:num_hours] >= 0)  # 逆变器功率流
    @variable(model, P_pv_mw[1:npv,1:num_hours] >= 0)  # 光伏功率流
    @variable(model, soc[1:ns,1:num_hours])  # 储能功率流
    @variable(model, ess_charge[1:ns,1:num_hours])  # 储能功率流
    @variable(model, ess_discharge[1:ns,1:num_hours])  # 储能功率流
    @variable(model, 0 <= ess_mode[1:ns, 1:num_hours] <= 1)  # Relaxed binary
    # @variable(model, ess_mode[1:ns, 1:num_hours], Bin)  # 严格二进制变量

    # 光伏功率约束
    for t in 1:num_hours
        @constraint(model, 0 .<= P_pv_mw[:,t] .<= pv_max_p_mw[:].*pv_max_p_mw_ratio[:,t])  # 光伏功率注入下限
    end

    # 交流发电机功率约束
    # for i in 1:ng
    #     @constraint(model, Pgen[i, :] >= 0)  # 发电机功率注入下限
    # end

    # 换流器功率约束
    # for t in 1:num_hours
    #     @constraint(model, Pij_inv[:,t].*Pij_rec[:,t].==0)  # 逆变器功率注入下限
    # end

    # 换流器有功功率约束
    for t in 1:num_hours
         @constraint(model, Pij[converter_indices,t] .== 0)
    end
   

    # 储能功率约束
    for e in 1:ns
        for t in 1:num_hours
            # 充放电功率限制
            @constraint(model, ess_charge[e, t] >= 0)
            @constraint(model, ess_discharge[e, t] >= 0)
            @constraint(model, ess_charge[e, t] <= ess_power_capacity_mw[e])
            @constraint(model, ess_discharge[e, t] <= ess_power_capacity_mw[e])
            
            # 充放电互斥约束
            @constraint(model, ess_discharge[e, t] <= ess_power_capacity_mw[e] * (1 - ess_mode[e, t]))
            @constraint(model, ess_charge[e, t] <= ess_power_capacity_mw[e] * ess_mode[e, t])
        end
        
        # 初始SOC约束
        @constraint(model, soc[e, 1] == ess_initial_soc[e] - 
                        ess_discharge[e, 1] / ess_efficiency[e] / ess_energy_capacity_mwh[e] + 
                        ess_charge[e, 1] * ess_efficiency[e] / ess_energy_capacity_mwh[e])
        
        # SOC随时间的演变
        for t in 2:num_hours
            @constraint(model, soc[e, t] == soc[e, t-1] - 
                            ess_discharge[e, t] / ess_efficiency[e] / ess_energy_capacity_mwh[e] + 
                            ess_charge[e, t] * ess_efficiency[e] / ess_energy_capacity_mwh[e])
        end
        
        # SOC上下限约束
        for t in 1:num_hours
            @constraint(model, soc[e, t] >= ess_min_soc[e])
            @constraint(model, soc[e, t] <= ess_max_soc[e])
        end
        
        # 循环约束：最终SOC等于初始SOC
        @constraint(model, soc[e, num_hours] == ess_initial_soc[e])
    end

   # 功率平衡约束
    for t in 1:num_hours
        @constraint(model,A' * Pij[:,t] + Cld_ac*loadAC_PD[:,t] + Cld_dc*loadDC_PD[:,t] -Cgen_ac*Pgen[:,t] -Cpv_ac*(pv_ac_p_mw.*pv_ac_p_mw_ratio[:,t]) + Cconv_ac*(Pij_inv[:,t]-η_rec.*Pij_rec[:,t]) + Cconv_dc*(Pij_rec[:,t]-η_inv.*Pij_inv[:,t]) - Cpv_dc * P_pv_mw[:,t] + Cstorage_ac*(ess_charge[:,t]-ess_discharge[:,t]) .== 0)
    end

    original_objective = sum(day_price_line[t,2] * Pgen[g, t] for g in 1:ng for t in 1:num_hours)
    penalty_weight = 1000  # 惩罚权重，需要调整
    # 目标函数：最小化总功率流
    # @objective(model, Min, sum(day_price_line[t,2] * Pgen[g, t] for g in 1:ng for t in 1:num_hours))
    if nc > 0 && num_hours > 0
        complementarity_penalty = sum(Pij_inv[i,t] * Pij_rec[i,t] for i in 1:nc for t in 1:num_hours)
        @objective(model, Min, original_objective + penalty_weight * complementarity_penalty)
    else
        @objective(model, Min, original_objective)
        println("警告: nc 或 num_hours 为零，跳过互补性惩罚项")
    end
    # 目标函数：最小化交流支路功率流的平方和
    # @objective(model, Min, sum(Pij[i,t]^2 *r[i]  for i in acbranch_indices for t in 1:num_hours)  )
    

    # 求解模型
    optimize!(model)

    # 检查求解状态
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED
        if status == MOI.OPTIMAL
            println("优化成功，找到最优解！")
        elseif status == MOI.LOCALLY_SOLVED
            println("优化成功，找到局部最优解！")
        elseif status == MOI.ALMOST_LOCALLY_SOLVED
            println("优化成功，找到近似最优解！")
        end
        
        # # 输出目标函数值
        obj_value = objective_value(model)
        # println("目标函数值: ", obj_value)
        
        # # 输出关键变量值
        # println("\n发电机出力 (Pgen):")
        Pgen_values = value.(Pgen)
        # display(Pgen_values)

        # println("\n换流器功率 (Pij_inv):")
        Pij_inv_values = value.(Pij_inv)
        # display(Pij_inv_values)

        # println("\n换流器功率 (Pij_rec):")
        Pij_rec_values = value.(Pij_rec)
        # display(Pij_rec_values)
        
        # println("\n光伏出力 (P_pv_mw):")
        P_pv_values = value.(P_pv_mw)
        # display(P_pv_values)
        
        # println("\n储能充电功率 (ess_charge):")
        ess_charge_values = value.(ess_charge)
        # display(ess_charge_values)
        
        # println("\n储能放电功率 (ess_discharge):")
        ess_discharge_values = value.(ess_discharge)
        # display(ess_discharge_values)
        
        # println("\n储能SOC (soc):")
        soc_values = value.(soc)
        # display(soc_values)
        
        # 返回结果
        return Dict(
            "status" => status,
            "objective" => obj_value,
            "Pgen" => Pgen_values,
            "Pij_inv" => Pij_inv_values,
            "Pij_rec" => Pij_rec_values,
            "P_pv_mw" => P_pv_values,
            "ess_charge" => ess_charge_values,
            "ess_discharge" => ess_discharge_values,
            "soc" => soc_values
        )
    else
        println("优化未成功完成。状态: ", status)
        println("求解器消息: ", raw_status(model))
        
        # 返回错误状态
        return Dict("status" => status, "error" => raw_status(model))
    end
end

function build_incidence_matrix_td(n_nodes, branchAC, branchDC, converter)
    # 计算总支路数
    n_branches = size(branchAC, 1) + size(branchDC, 1) + size(converter, 1)
    
    # 提取所有支路的起始和终止节点，并标记支路类型
    # 类型标记: 1=AC, 2=DC, 3=Converter
    branch_data = Vector{Tuple{Int, Int, Int, Int}}(undef, n_branches)
    
    # 添加交流支路
    for i in 1:size(branchAC, 1)
        # 确保小的节点号作为起始节点
        node1 = Int(branchAC[i, 1])
        node2 = Int(branchAC[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[i] = (from, to, 1, i)  # (起始节点, 终止节点, 类型, 原始索引)
    end
    
    # 添加直流支路
    offset = size(branchAC, 1)
    for i in 1:size(branchDC, 1)
        # 确保小的节点号作为起始节点
        node1 = Int(branchDC[i, 1])
        node2 = Int(branchDC[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[offset + i] = (from, to, 2, i)
    end
    
    # 添加变流器
    offset = size(branchAC, 1) + size(branchDC, 1)
    for i in 1:size(converter, 1)
        # 确保小的节点号作为起始节点
        node1 = Int(converter[i, 1])
        node2 = Int(converter[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[offset + i] = (from, to, 3, i)
    end
    
    # 按起始节点排序
    sort!(branch_data, by = x -> (x[1], x[2]))
    
    # 创建关联矩阵
    A = zeros(n_branches, n_nodes)
    
    # 根据排序后的支路构建关联矩阵
    for (idx, (from, to, type, _)) in enumerate(branch_data)
        A[idx, from] = 1   # 从节点流出为正
        A[idx, to] = -1    # 到节点流入为负
    end
    
    return A, branch_data
end