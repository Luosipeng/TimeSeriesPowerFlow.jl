using  JuMP
using  Ipopt
using  LinearAlgebra
using  Gurobi
function runlindistflow(new_jpc, opt, 
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
        storage_min_p_mw)

    branchAC = new_jpc.branchAC
    branchDC = new_jpc.branchDC
    converter = new_jpc.converter
    # loadAC = new_jpc.loadAC
    # loadDC = new_jpc.loadDC
    # pvarray = new_jpc.pv
    Cld_ac = Cld_ac[2:end, :]  # 去除第一行
    Cld_dc = Cld_dc[2:end, :]  # 去除第一行
    Cgen_ac = Cgen_ac[2:end, :]  # 去除第一行
    Cconv_ac = Cconv_ac[2:end, :]  # 去除第一行
    Cconv_dc = Cconv_dc[2:end, :]  # 去除第一行
    Cpv_dc = Cpv_dc[2:end, :]  # 去除第一行
    Cstorage_ac = Cstorage_ac[2:end, :]  # 去除第一行

    n_nodes = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)

    A , branch_data = build_incidence_matrix(n_nodes, branchAC, branchDC, converter)
    
    
    # 初始化
    nbr = size(A, 1)  # 支路数
    nc = size(converter, 1)  # 变流器数
    npv = size(new_jpc.pv, 1)  # 光伏设备数
    ns = size(new_jpc.storage, 1)  # 储能设备数

    # 搜寻支路索引
    # dcbranch = filter(x -> x[3] == 2, branch_data)
    dcbranch_indices = findall(x -> x[3] == 2, branch_data)
    converter_indices = findall(x -> x[3] == 3, branch_data)
    acbranch_indices = findall(x -> x[3] == 1, branch_data)


#     model = Model(Gurobi.Optimizer)
#     set_optimizer_attribute(model, "OutputFlag", 1)
#     # model = Model(Ipopt.Optimizer)
#     # IPOPT特定参数设置
#     # set_optimizer_attribute(model, "print_level", 5)  # 输出详细程度
#     # set_optimizer_attribute(model, "max_iter", 3000)  # 最大迭代次数
#     # set_optimizer_attribute(model, "tol", 1e-6)       # 收敛容差

#     @variable(model, Pij[1:nbr])  # 支路有功功率
#     @variable(model, Qij[1:nbr])  # 支路无功功率
#     @variable(model, Pij_inv[1:nc])  # 换流器功率
#     @variable(model, Qij_inv[1:nc])  # 换流器无功功率
#     @variable(model, Pij_rec[1:nc])  # 直流侧功率
#     @variable(model, P_pv_mw[1:npv]) 
#     @variable(model, Pess_mw[1:ns])  # 储能设备有功功率

#     # 添加绝对值变量
#     @variable(model, Pij_abs[1:nbr] >= 0)  # 支路有功功率绝对值
#     @variable(model, Qij_abs[1:nbr] >= 0)  # 支路无功功率绝对值

#     # 添加绝对值约束
#     for i in 1:nbr
#         @constraint(model, Pij_abs[i] >= Pij[i])
#         @constraint(model, Pij_abs[i] >= -Pij[i])
#         @constraint(model, Qij_abs[i] >= Qij[i])
#         @constraint(model, Qij_abs[i] >= -Qij[i])
#     end

#     # 使用for循环对每个PV单独处理功率约束
#     for i in 1:npv
#         @constraint(model, 0 <= P_pv_mw[i] <= pv_max_p_mw[i])  # 光伏功率注入下限
#     end

#     # 储能出力约束
#     for i in 1:ns
#         @constraint(model, storage_min_p_mw[i] <= Pess_mw[i] <= storage_max_p_mw[i])  # 储能有功功率约束
#     end

#     #换流器约束
#     for i in 1:nc
#         @constraint(model, Pij_inv[i] >= 0)  # 换流器有功功率约束
#         @constraint(model, Qij_inv[i] >= 0)  # 换流器无功功率约束
#         @constraint(model, Pij_rec[i] >= 0)  # 直流侧功率约束

#         @constraint(model, Pij_inv[i] * Pij_rec[i] == 0)  # 换流器有功功率和直流侧功率的乘积为负
#     end

#     # 直流无功约束
#     @constraint(model, Qij[dcbranch_indices] .== 0.0)

#     # 变流器有功功率约束
#     @constraint(model, Pij[converter_indices] .== 0)  # 交流侧有功功率
#     @constraint(model, Qij[converter_indices] .== 0)  # 交流侧无功功率

#     # 有功功率平衡约束
#     @constraint(model, A' * Pij + Cld_ac*loadAC_PD + Cld_dc*loadDC_PD + Cconv_ac*(Pij_inv-η_rec.*Pij_rec) + Cconv_dc*(Pij_rec-η_inv.*Pij_inv) - Cpv_dc * P_pv_mw - Cstorage_ac*Pess_mw .== 0)
#     # 无功功率平衡约束
#     @constraint(model, A' * Qij + Cld_ac*loadAC_QD + Cconv_ac * Qij_inv .== 0)

#     # 目标函数
#     # @objective(model, Min, sum(Pij_abs) )
#     @objective(model, Min, 
#         sum(Pij_abs[acbranch_indices])+0.1*sum(Pij_abs[dcbranch_indices]) - 0.5 * sum(P_pv_mw) -0.3 * sum(Pess_mw)
#     )

#     # 求解模型
#     optimize!(model)

#     # 检查模型是否成功求解
#    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
#         println("模型成功求解！")
        
#         println("\n===== 约束满足情况 =====")
        
#         # 检查绝对值约束
#         println("\n绝对值约束：")
#         for i in 1:nbr
#             p_val = value(Pij[i])
#             p_abs_val = value(Pij_abs[i])
#             q_val = value(Qij[i])
#             q_abs_val = value(Qij_abs[i])
#             println("支路 $i: |Pij| = $p_abs_val (Pij = $p_val), |Qij| = $q_abs_val (Qij = $q_val)")
#         end
        
#         # 检查PV功率约束
#         println("\nPV功率约束：")
#         for i in 1:npv
#             p_pv = value(P_pv_mw[i])
#             println("PV $i: 0 <= $p_pv <= $(pv_max_p_mw[i])")
#         end
        
#         # 检查储能功率约束
#         println("\n储能功率约束：")
#         for i in 1:ns
#             p_ess = value(Pess_mw[i])
#             println("储能 $i: $(storage_min_p_mw[i]) <= $p_ess <= $(storage_max_p_mw[i])")
#         end
        
#         # 检查换流器约束
#         println("\n换流器约束：")
#         for i in 1:nc
#             p_inv = value(Pij_inv[i])
#             q_inv = value(Qij_inv[i])
#             p_rec = value(Pij_rec[i])
#             product = p_inv * p_rec
#             println("换流器 $i: Pij_inv = $p_inv (>= 0), Qij_inv = $q_inv (>= 0), Pij_rec = $p_rec (>= 0)")
#             println("   Pij_inv * Pij_rec = $product (应为0)")
#         end
        
#         # 检查直流无功约束
#         println("\n直流无功约束：")
#         for i in dcbranch_indices
#             q_val = value(Qij[i])
#             println("直流支路 $i: Qij = $q_val (应为0)")
#         end
        
#         # 检查变流器功率约束
#         println("\n变流器功率约束：")
#         for i in converter_indices
#             p_val = value(Pij[i])
#             q_val = value(Qij[i])
#             println("变流器支路 $i: Pij = $p_val (应为0), Qij = $q_val (应为0)")
#         end
        
#         # 检查功率平衡约束
#         println("\n功率平衡约束：")
#         p_balance = value.(A' * Pij + Cld_ac*loadAC_PD + Cld_dc*loadDC_PD + Cconv_ac*(Pij_inv-η_rec.*Pij_rec) + Cconv_dc*(Pij_rec-η_inv.*Pij_inv) - Cpv_dc * P_pv_mw - Cstorage_ac*Pess_mw)
#         q_balance = value.(A' * Qij + Cld_ac*loadAC_QD + Cconv_ac * Qij_inv)
#         println("有功功率平衡: $(norm(p_balance, Inf))")
#         println("无功功率平衡: $(norm(q_balance, Inf))")
        
#         # 输出详细结果
#         println("\n===== 详细结果 =====")
#         println("\n支路有功功率 Pij:")
#         for i in 1:nbr
#             println("支路 $i: $(value(Pij[i]))")
#         end
        
#         println("\n支路无功功率 Qij:")
#         for i in 1:nbr
#             println("支路 $i: $(value(Qij[i]))")
#         end
        
#         println("\n换流器功率:")
#         for i in 1:nc
#             println("换流器 $i: Pij_inv = $(value(Pij_inv[i])), Qij_inv = $(value(Qij_inv[i])), Pij_rec = $(value(Pij_rec[i]))")
#         end
        
#         println("\nPV功率:")
#         for i in 1:npv
#             println("PV $i: $(value(P_pv_mw[i]))")
#         end
        
#         println("\n储能功率:")
#         for i in 1:ns
#             println("储能 $i: $(value(Pess_mw[i]))")
#         end
        
#         # 输出目标函数值
#         println("\n目标函数值: $(objective_value(model))")
        
#     else
#         @error "模型求解失败" status=termination_status(model)
#     end


    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "FeasibilityTol", 1e-5)
    # model = Model(Ipopt.Optimizer)
    # IPOPT特定参数设置
    # set_optimizer_attribute(model, "print_level", 5)  # 输出详细程度
    # set_optimizer_attribute(model, "max_iter", 3000)  # 最大迭代次数
    # set_optimizer_attribute(model, "tol", 1e-6)       # 收敛容差

    @variable(model, Pij[1:nbr])  # 支路有功功率
    @variable(model, Qij[1:nbr])  # 支路无功功率
    @variable(model, Pij_inv[1:nc] >= 0)  # 换流器功率
    @variable(model, Qij_inv[1:nc] >= 0)  # 换流器无功功率
    @variable(model, Pij_rec[1:nc] >= 0)  # 直流侧功率
    @variable(model, P_pv_mw[1:npv]) 
    @variable(model, Pess_mw[1:ns])  # 储能设备有功功率

    # 添加绝对值变量
    @variable(model, Pij_abs[1:nbr] >= 0)  # 支路有功功率绝对值
    @variable(model, Qij_abs[1:nbr] >= 0)  # 支路无功功率绝对值

    # 添加绝对值约束
    for i in 1:nbr
        @constraint(model, Pij_abs[i] >= Pij[i])
        @constraint(model, Pij_abs[i] >= -Pij[i])
        @constraint(model, Qij_abs[i] >= Qij[i])
        @constraint(model, Qij_abs[i] >= -Qij[i])
    end

    # 使用for循环对每个PV单独处理功率约束
    for i in 1:npv
        @constraint(model, 0 <= P_pv_mw[i] <= pv_max_p_mw[i])  # 光伏功率注入下限
    end

    # 储能出力约束
    for i in 1:ns
        @constraint(model, storage_min_p_mw[i] <= Pess_mw[i] <= storage_max_p_mw[i])  # 储能有功功率约束
    end

    # 使用大M法替代双线性约束
    # 定义一个足够大的常数M
    M = 1000.0  # 根据您的问题规模调整这个值

    # 添加二进制变量表示换流器工作模式
    @variable(model, z[1:nc], Bin)

    # 换流器约束 - 使用大M法
    for i in 1:nc
        # 如果z=1，允许Pij_inv取值，否则Pij_inv必须为0
        @constraint(model, Pij_inv[i] <= M * z[i])
        
        # 如果z=0，允许Pij_rec取值，否则Pij_rec必须为0
        @constraint(model, Pij_rec[i] <= M * (1 - z[i]))
        
        # 不再需要原来的双线性约束 Pij_inv[i] * Pij_rec[i] == 0
    end

    # 直流无功约束
    @constraint(model, Qij[dcbranch_indices] .== 0.0)

    # 变流器有功功率约束
    @constraint(model, Pij[converter_indices] .== 0)  # 交流侧有功功率
    @constraint(model, Qij[converter_indices] .== 0)  # 交流侧无功功率

    # 有功功率平衡约束
    @constraint(model, A' * Pij + Cld_ac*loadAC_PD + Cld_dc*loadDC_PD -Cgen_ac*genAC_PG + Cconv_ac*(Pij_inv-η_rec.*Pij_rec) + Cconv_dc*(Pij_rec-η_inv.*Pij_inv) - Cpv_dc * P_pv_mw - Cstorage_ac*Pess_mw .== 0)
    # 无功功率平衡约束
    @constraint(model, A' * Qij + Cld_ac*loadAC_QD + Cconv_ac * Qij_inv .== 0)

    # 目标函数
    @objective(model, Min, 
        sum(Pij_abs[acbranch_indices])+0.1*sum(Pij_abs[dcbranch_indices]) - 0.5 * sum(P_pv_mw) -0.3 * sum(Pess_mw)
    )

    # 求解模型
    optimize!(model)

    # 检查模型是否成功求解
    if termination_status(model) == MOI.OPTIMAL
        println("模型成功求解！")
        
        # println("\n===== 约束满足情况 =====")
        
        # # 检查绝对值约束
        # println("\n绝对值约束：")
        # for i in 1:nbr
        #     p_val = value(Pij[i])
        #     p_abs_val = value(Pij_abs[i])
        #     q_val = value(Qij[i])
        #     q_abs_val = value(Qij_abs[i])
        #     println("支路 $i: |Pij| = $p_abs_val (Pij = $p_val), |Qij| = $q_abs_val (Qij = $q_val)")
        # end
        
        # # 检查PV功率约束
        # println("\nPV功率约束：")
        # for i in 1:npv
        #     p_pv = value(P_pv_mw[i])
        #     println("PV $i: 0 <= $p_pv <= $(pv_max_p_mw[i])")
        # end
        
        # # 检查储能功率约束
        # println("\n储能功率约束：")
        # for i in 1:ns
        #     p_ess = value(Pess_mw[i])
        #     println("储能 $i: $(storage_min_p_mw[i]) <= $p_ess <= $(storage_max_p_mw[i])")
        # end
        
        # # 检查换流器约束（大M法）
        # println("\n换流器约束（大M法）：")
        # for i in 1:nc
        #     p_inv = value(Pij_inv[i])
        #     q_inv = value(Qij_inv[i])
        #     p_rec = value(Pij_rec[i])
        #     z_val = value(z[i])
        #     println("换流器 $i: Pij_inv = $p_inv, Qij_inv = $q_inv, Pij_rec = $p_rec, z = $z_val")
        #     println("   互斥约束: z=$z_val -> " * (z_val > 0.5 ? "逆变模式" : "整流模式"))
        # end
        
        # # 检查直流无功约束
        # println("\n直流无功约束：")
        # for i in dcbranch_indices
        #     q_val = value(Qij[i])
        #     println("直流支路 $i: Qij = $q_val (应为0)")
        # end
        
        # # 检查变流器功率约束
        # println("\n变流器功率约束：")
        # for i in converter_indices
        #     p_val = value(Pij[i])
        #     q_val = value(Qij[i])
        #     println("变流器支路 $i: Pij = $p_val (应为0), Qij = $q_val (应为0)")
        # end
        
        # # 检查功率平衡约束
        # println("\n功率平衡约束：")
        # p_balance = value.(A' * Pij + Cld_ac*loadAC_PD + Cld_dc*loadDC_PD + Cconv_ac*(Pij_inv-η_rec.*Pij_rec) + Cconv_dc*(Pij_rec-η_inv.*Pij_inv) - Cpv_dc * P_pv_mw - Cstorage_ac*Pess_mw)
        # q_balance = value.(A' * Qij + Cld_ac*loadAC_QD + Cconv_ac * Qij_inv)
        # println("有功功率平衡: $(norm(p_balance, Inf))")
        # println("无功功率平衡: $(norm(q_balance, Inf))")
        
        # # 输出详细结果
        # println("\n===== 详细结果 =====")
        # println("\n支路有功功率 Pij:")
        # for i in 1:nbr
        #     println("支路 $i: $(value(Pij[i]))")
        # end
        
        # println("\n支路无功功率 Qij:")
        # for i in 1:nbr
        #     println("支路 $i: $(value(Qij[i]))")
        # end
        
        println("\n换流器功率:")
        for i in 1:nc
            println("换流器 $i: Pij_inv = $(value(Pij_inv[i])), Qij_inv = $(value(Qij_inv[i])), Pij_rec = $(value(Pij_rec[i]))")
        end
        
        println("\nPV功率:")
        for i in 1:npv
            println("PV $i: $(value(P_pv_mw[i]))")
        end
        
        println("\n储能功率:")
        for i in 1:ns
            println("储能 $i: $(value(Pess_mw[i]))")
        end
        
        # 输出目标函数值
        println("\n目标函数值: $(objective_value(model))")
        
    else
        @error "模型求解失败" status=termination_status(model)
    end

    return value.(Pij_inv), value.(Qij_inv), value.(Pij_rec), value.(Pess_mw)

end

function build_incidence_matrix(n_nodes, branchAC, branchDC, converter, ref_node=1)
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
    
    # 去掉参考节点对应的列
    return A[:, setdiff(1:n_nodes, ref_node)], branch_data
end

# 打印所有约束的详细版本
function print_all_constraints_detailed(model)
    println("===== 模型约束详情 =====")
    
    # 打印直流无功约束
    println("\n1. 直流无功约束:")
    for (i, idx) in enumerate(dcbranch_indices)
        println("   Qn[$idx] == 0.0")
    end
    
    # 打印变流器有功功率约束
    println("\n2. 变流器有功功率约束:")
    for (i, idx) in enumerate(converter_indices)
        println("   Pn[$idx] == 0")
    end
    
    # 打印变流器无功功率约束
    println("\n3. 变流器无功功率约束:")
    for (i, idx) in enumerate(converter_indices)
        println("   Qn[$idx] == 0")
    end
    
    # 打印变流器条件约束
    println("\n4. 变流器条件约束:")
    for i in 1:nc
        println("   变流器 $i:")
        println("     - 如果 z[$i] = 1, 则 Pac[$i] >= 0")
        println("     - 如果 z[$i] = 0, 则 Pac[$i] <= -$ε")
        println("     - 如果 z[$i] = 1, 则 Pac[$i] = -$(efficiency[i]) * Pdc[$i]")
        println("     - 如果 z[$i] = 0, 则 Pac[$i] = -Pdc[$i] / $(efficiency[i])")
    end
    
    # 打印有功功率平衡约束 - 详细展开每个节点的约束
    println("\n5. 有功功率平衡约束:")
    A_transpose = A'
    for i in 1:size(A_transpose, 1)
        constraint_str = "   节点 $i: "
        
        # 添加支路功率项
        branch_terms = []
        for j in 1:nbr
            if A_transpose[i, j] != 0
                push!(branch_terms, "$(A_transpose[i, j]) * Pn[$j]")
            end
        end
        
        if !isempty(branch_terms)
            constraint_str *= join(branch_terms, " + ")
        else
            constraint_str *= "0"
        end
        
        # 添加负荷功率
        constraint_str *= " + Pd[$i]"
        
        # 添加变流器交流侧功率
        inv_terms = []
        for j in 1:nc
            if Cinv[i, j] != 0
                push!(inv_terms, "$(Cinv[i, j]) * Pac[$j]")
            end
        end
        
        if !isempty(inv_terms)
            constraint_str *= " + " * join(inv_terms, " + ")
        end
        
        # 添加变流器直流侧功率
        rec_terms = []
        for j in 1:nc
            if Crec[i, j] != 0
                push!(rec_terms, "$(Crec[i, j]) * Pdc[$j]")
            end
        end
        
        if !isempty(rec_terms)
            constraint_str *= " + " * join(rec_terms, " + ")
        end
        
        constraint_str *= " == 0"
        println(constraint_str)
    end
    
    # 打印无功功率平衡约束 - 详细展开每个节点的约束
    println("\n6. 无功功率平衡约束:")
    for i in 1:size(A_transpose, 1)
        constraint_str = "   节点 $i: "
        
        # 添加支路功率项
        branch_terms = []
        for j in 1:nbr
            if A_transpose[i, j] != 0
                push!(branch_terms, "$(A_transpose[i, j]) * Qn[$j]")
            end
        end
        
        if !isempty(branch_terms)
            constraint_str *= join(branch_terms, " + ")
        else
            constraint_str *= "0"
        end
        
        # 添加负荷功率
        constraint_str *= " + Qd[$i]"
        
        # 添加变流器交流侧无功功率
        inv_terms = []
        for j in 1:nc
            if Cinv[i, j] != 0
                push!(inv_terms, "$(Cinv[i, j]) * Qac[$j]")
            end
        end
        
        if !isempty(inv_terms)
            constraint_str *= " + " * join(inv_terms, " + ")
        end
        
        constraint_str *= " == 0"
        println(constraint_str)
    end
end

# 调用函数打印所有约束
# print_all_constraints_detailed(model)

"""
    process_pv_devices(pvarray, new_jpc, indices; a=14.5, base_mva=100)

处理光伏设备并计算其功率注入。

# 参数
- `pvarray`: 光伏设备数组，包含所有光伏设备的参数
- `new_jpc`: 电力系统数据结构
- `indices`: 包含各列索引的命名元组或字典，例如 (PV_IN_SERVICE=1, PV_BUS=2, ...)
- `a`: 经验公式参数，默认为14.5
- `base_mva`: 基准功率，单位为MVA，默认为100

# 返回值
- `Spv`: 复数功率注入向量，表示每个节点的光伏功率注入
"""
function process_pv_devices(pvarray, new_jpc, n_nodes; a=14.5, base_mva=100)
    if size(pvarray, 1) < 1
        Ppv_max = zeros(Float64, size(pvarray,1))
        Ppv_min = zeros(Float64, size(pvarray,1))
        Cpv = zeros(n_nodes, size(pvarray, 1))
        return Ppv_max,Ppv_min , Cpv
    end
    
    # 处理光伏设备
    active_pv = pvarray[pvarray[:, PV_IN_SERVICE] .> 0, :]
    
    # 如果没有活跃的光伏设备，返回零向量
    if isempty(active_pv)
        return zeros(Complex{Float64}, size(new_jpc.busDC, 1))
    end
    
    pvbus = Int.(active_pv[:, PV_BUS])  # 光伏连接的母线编号
    bus_indice = findall(x -> x in pvbus, new_jpc.busDC[:,BUS_I])  # 光伏连接的母线索引
    
    # 如果没有找到匹配的母线，返回零向量
    if isempty(bus_indice)
        return zeros(Complex{Float64}, size(new_jpc.busDC, 1))
    end
    
    base_kvs = new_jpc.busDC[bus_indice, BASE_KV]
    V_max = 1.05 .*ones(length(bus_indice))  # 假设光伏设备的电压为1.0 pu
    V_min = 0.95 .*ones(length(bus_indice))  # 假设光伏设备的电压为1.0 pu

    # 提取光伏参数
    Vocs = active_pv[:, PV_VOC]      # 开路电压 (V)
    Iscs = active_pv[:, PV_ISC]      # 短路电流 (A)
    areas = active_pv[:, PV_AREA]    # 面积或其他缩放因子

    # 将标幺值电压转换为实际电压 (V)
    voltage_ratios = 1
    V_arrays_max = V_max .* base_kvs .* 1000 .* voltage_ratios  # 转换为伏特
    V_arrays_min = V_min .* base_kvs .* 1000 .* voltage_ratios  # 转换为伏特
    
    # 初始化电流
    I_arrays_max = zeros(size(active_pv, 1))
    I_arrays_min = zeros(size(active_pv, 1))

    # 创建掩码进行向量化操作
    valid_v_mask_max = (V_arrays_max .>= 0) .& (V_arrays_max .<= Vocs)
    valid_v_mask_min = (V_arrays_min .>= 0) .& (V_arrays_min .<= Vocs)
    
    # 计算比例
    ratiosmax = zeros(size(V_arrays_max))
    ratiosmax[valid_v_mask_max] = V_arrays_max[valid_v_mask_max] ./ Vocs[valid_v_mask_max]
    ratiosmin = zeros(size(V_arrays_min))
    ratiosmin[valid_v_mask_min] = V_arrays_min[valid_v_mask_min] ./ Vocs[valid_v_mask_min]
    
    # 使用经验公式计算电流（向量化操作）
    I_arrays_max[valid_v_mask_max] = Iscs[valid_v_mask_max] .* areas[valid_v_mask_max] .* (1 .- ratiosmax[valid_v_mask_max].^a)
    I_arrays_min[valid_v_mask_min] = Iscs[valid_v_mask_min] .* areas[valid_v_mask_min] .* (1 .- ratiosmin[valid_v_mask_min].^a)

    P_arrays_max = V_arrays_max .* I_arrays_max  # 最大功率
    P_arrays_min = V_arrays_min .* I_arrays_min  # 最小功率
    P_pus_max = P_arrays_max ./ (1e6)  # 从W转换为标幺值
    P_pus_min = P_arrays_min ./ (1e6)  # 从W转换为标幺值
    Ppv_max = -P_pus_max  # 负值表示功率注入
    Ppv_min = -P_pus_min  # 负值表示功率注入
    # 获取节点数量
    
    # 初始化功率注入向量
    Cpv = zeros(n_nodes, size(pvarray, 1))
    # 更新功率注入
    for i in 1:length(pvbus)
        Cpv[pvbus[i],i] = 1.0
    end
    
    return Ppv_max, Ppv_min, Cpv
end




