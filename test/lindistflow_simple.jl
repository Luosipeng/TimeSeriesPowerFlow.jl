"""
This file is used to test the linear AC/DC hybrid distribution flow model.
"""

using JuMP
using Gurobi  # 替换Ipopt为Gurobi
using LinearAlgebra

# 定义配电网络参数
# 改进后的连续换流器模型部分
function create_distflow_model_no_slack(
    Nac::Int,                  # 交流节点数量
    Ndc::Int,                  # 直流节点数量
    C_inv::Vector{Int},        # 逆变器位置
    C_rec::Vector{Int},        # 整流器位置
    P_inv0::Vector{Float64},   # 逆变器有功功率初始值
    Q_inv0::Vector{Float64},   # 逆变器无功功率初始值
    P_rec0::Vector{Float64},   # 整流器有功功率初始值
    acP_load::Vector{Float64}, # 交流负载有功功率
    acQ_load::Vector{Float64}, # 交流负载无功功率
    dcP_load::Vector{Float64}, # 直流负载有功功率
    v_min::Float64,            # 最小电压限制
    v_max::Float64,            # 最大电压限制
    acbranch_idx,
    dcbranch_idx,
    acbranches,
    dcbranches,
    Aac,
    Adc,
    Rac,
    Xac,
    Rdc,
)
    # 创建优化模型，使用Gurobi求解器
    model = Model(Gurobi.Optimizer)
    
    # 设置Gurobi参数 - 提高精度和稳定性
    set_optimizer_attribute(model, "OutputFlag", 1)       # 显示求解过程
    set_optimizer_attribute(model, "MIPGap", 1e-6)        # MIP相对间隙，提高精度
    set_optimizer_attribute(model, "TimeLimit", 600)      # 求解时间限制（秒），增加时间
    set_optimizer_attribute(model, "Threads", 0)          # 使用所有可用线程
    set_optimizer_attribute(model, "NumericFocus", 3)     # 最高数值稳定性
    set_optimizer_attribute(model, "FeasibilityTol", 1e-9) # 可行性容差，提高精度
    set_optimizer_attribute(model, "OptimalityTol", 1e-9) # 最优性容差，提高精度
    set_optimizer_attribute(model, "BarConvTol", 1e-10)   # 障碍法收敛容差
    set_optimizer_attribute(model, "ScaleFlag", 2)        # 积极缩放
     
    # 决策变量
    @variable(model, Pac[1:length(acbranches)]) # 支路有功功率
    @variable(model, Qac[1:length(acbranches)]) # 支路无功功率
    @variable(model, Pdc[1:length(dcbranches)]) # 直流支路有功功率
    
    # 节点电压变量 - 注意根节点电压固定为1.0
    @variable(model, vac[1:Nac])  # 交流节点电压幅值
    @variable(model, vdc[1:Ndc])  # 直流节点电压幅值

    # 换流器功率变量 - 使用与初始值相同的维度
    @variable(model, P_inv[1:length(P_inv0)])
    @variable(model, Q_inv[1:length(Q_inv0)])
    @variable(model, P_rec[1:length(P_rec0)])
    
    # 为非根节点设置电压约束
    for i in 2:Nac
        @constraint(model, v_min^2 <= vac[i] <= v_max^2)
    end
    
    for i in 2:Ndc
        @constraint(model, v_min^2 <= vdc[i] <= v_max^2)
    end
    
    # 固定根节点电压为1.0的平方
    @constraint(model, vac[1] == 1.0)
    @constraint(model, vdc[1] == 1.0)

    # 目标函数：最小化与初始值的最大偏差（无穷范数）
    @variable(model, max_dev >= 0)  # 最大偏差变量
    
    # 对每个换流器功率设置偏差约束
    for i in 1:length(P_inv)
        @constraint(model, P_inv[i] - P_inv0[i] <= max_dev)
        @constraint(model, P_inv0[i] - P_inv[i] <= max_dev)
    end
    
    for i in 1:length(Q_inv)
        @constraint(model, Q_inv[i] - Q_inv0[i] <= max_dev)
        @constraint(model, Q_inv0[i] - Q_inv[i] <= max_dev)
    end
    
    for i in 1:length(P_rec)
        @constraint(model, P_rec[i] - P_rec0[i] <= max_dev)
        @constraint(model, P_rec0[i] - P_rec[i] <= max_dev)
    end
    
    # 最小化最大偏差
    @objective(model, Min, max_dev)
    
    # 创建逆变器功率向量（使用JuMP表达式）
    P_inverter = zeros(AffExpr, Nac-1)
    Q_inverter = zeros(AffExpr, Nac-1)
    
    # 将逆变器功率映射到对应节点
    for i in 1:length(C_inv)
        if C_inv[i] > 1  # 确保不是根节点
            node_idx = C_inv[i] - 1  # 调整索引（去掉根节点）
            P_inverter[node_idx] = P_inv[i]
            Q_inverter[node_idx] = Q_inv[i]
        end
    end

    # 创建整流器功率向量（使用JuMP表达式）
    P_rectifier = zeros(AffExpr, Ndc-1)
    
    # 将整流器功率映射到对应节点
    for i in 1:length(C_rec)
        if C_rec[i] > 1  # 确保不是根节点
            node_idx = C_rec[i] - 1  # 调整索引（去掉根节点）
            P_rectifier[node_idx] = P_rec[i]
        end
    end
    
    # 交流节点功率平衡方程 - 完全消除松弛变量，强制功率平衡
    # 确保维度匹配 - 只取非根节点的负载
    acP_load_non_root = acP_load[2:end]  # 去掉根节点
    acQ_load_non_root = acQ_load[2:end]  # 去掉根节点
    
    # 确保acP_n和acQ_n的维度正确
    # 首先检查分支索引是否正确
    branch_indices = [get(acbranch_idx, k, 0) for k in 2:Nac]
    if 0 in branch_indices
        error("无法找到某些节点的对应支路索引")
    end
    
    acP_n = Pac[branch_indices]
    acQ_n = Qac[branch_indices]
    
    # 确保维度匹配后再添加约束
    if length(acP_n) == length(acP_load_non_root) && length(acP_n) == length(P_inverter)
        @constraint(model, Aac' * acP_n + acP_load_non_root - P_inverter .== 0)
    else
        error("交流有功功率平衡约束维度不匹配: acP_n=$(length(acP_n)), acP_load_non_root=$(length(acP_load_non_root)), P_inverter=$(length(P_inverter))")
    end
    
    if length(acQ_n) == length(acQ_load_non_root) && length(acQ_n) == length(Q_inverter)
        @constraint(model, Aac' * acQ_n + acQ_load_non_root - Q_inverter .== 0)
    else
        error("交流无功功率平衡约束维度不匹配: acQ_n=$(length(acQ_n)), acQ_load_non_root=$(length(acQ_load_non_root)), Q_inverter=$(length(Q_inverter))")
    end

    # 直流节点功率平衡方程 - 完全消除松弛变量，强制功率平衡
    dcP_load_non_root = dcP_load[2:end]  # 去掉根节点
    
    # 确保dcP_n的维度正确
    branch_indices_dc = [get(dcbranch_idx, k, 0) for k in 2:Ndc]
    if 0 in branch_indices_dc
        error("无法找到某些直流节点的对应支路索引")
    end
    
    dcP_n = Pdc[branch_indices_dc]
    
    # 确保维度匹配后再添加约束
    if length(dcP_n) == length(dcP_load_non_root) && length(dcP_n) == length(P_rectifier)
        @constraint(model, Adc' * dcP_n + dcP_load_non_root + P_rectifier .== 0)
    else
        error("直流功率平衡约束维度不匹配: dcP_n=$(length(dcP_n)), dcP_load_non_root=$(length(dcP_load_non_root)), P_rectifier=$(length(P_rectifier))")
    end

    # 交流节点电压平衡方程
    @constraint(model, -2.0 * Rac * acP_n - 2.0 * Xac * acQ_n + vac[1:end-1] .== vac[2:end])

    # 直流节点电压平衡方程
    @constraint(model, -2.0 * Rdc * dcP_n + vdc[1:end-1] .== vdc[2:end])
    
    # 使用连续换流器效率模型替代二元变量
    # 定义换流器效率参数
    η_inv = 0.9  # 逆变器效率
    η_rec = 0.9  # 整流器效率
    
    # 使用辅助变量和平滑近似替代二元变量
    @variable(model, P_inv_pos[1:length(P_inv)] >= 0)  # 逆变器正向功率
    @variable(model, P_inv_neg[1:length(P_inv)] >= 0)  # 逆变器负向功率
    
    # 平滑近似参数
    smooth_factor = 1e-4  # 平滑因子，调整为更小的值以提高精度
    
    for i in 1:length(P_inv)
        # 分解P_inv为正负两部分
        @constraint(model, P_inv[i] == P_inv_pos[i] - P_inv_neg[i])
        
        # 平滑过渡的连续换流器模型
        @constraint(model, P_rec[i] == P_inv_pos[i] / η_inv + P_inv_neg[i] * η_rec)
    end
    
    # 增加一个额外的目标项，鼓励P_inv_pos和P_inv_neg不同时为大值
    @variable(model, complementarity_slack[1:length(P_inv)] >= 0)
    
    for i in 1:length(P_inv)
        @constraint(model, P_inv_pos[i] * P_inv_neg[i] <= complementarity_slack[i] + smooth_factor)
    end
    
    # 修改目标函数，加入对互补松弛的惩罚
    complementarity_weight = 1e-3  # 互补松弛项的权重，增加权重以更强调互补性
    @objective(model, Min, max_dev + complementarity_weight * sum(complementarity_slack))
    
    return model, Pac, Qac, Pdc, vac, vdc, P_inv, Q_inv, P_rec
end



# 构建网络关联矩阵
function build_incidence_matrix(n_nodes, branches)
    # 创建关联矩阵
    A = zeros(length(branches), n_nodes)
    
    for (i, (from, to)) in enumerate(branches)
        A[i, from] = 1  # 从节点流出为正
        A[i, to] = -1   # 到节点流入为负
    end
    
    return A[:, 2:end]  # 去掉根节点的列
end

function reindex_network_by_path(bus, branch, root)
    # 构建无向图
    n_nodes = size(bus, 1)
    graph = Dict{Int, Vector{Int}}()
    
    # 从支路信息构建图
    for i in eachindex(branch[:,1])
        from = Int(branch[i, 1])
        to = Int(branch[i, 2])
        
        if !haskey(graph, from)
            graph[from] = Int[]
        end
        if !haskey(graph, to)
            graph[to] = Int[]
        end
        
        push!(graph[from], to)
        push!(graph[to], from)  # 无向图
    end
    
    # 使用BFS进行遍历和重新编号
    old_to_new = Dict{Int, Int}()
    visited = Set{Int}()
    queue = [root]
    push!(visited, root)
    
    # 根节点编号为1
    old_to_new[root] = 1
    new_idx = 2
    
    while !isempty(queue)
        node = popfirst!(queue)
        
        # 遍历相邻节点
        for neighbor in get(graph, node, Int[])
            if !(neighbor in visited)
                push!(visited, neighbor)
                old_to_new[neighbor] = new_idx
                new_idx += 1
                push!(queue, neighbor)
            end
        end
    end
    
    # 检查是否所有节点都被访问到
    if length(old_to_new) != n_nodes
        unvisited = [i for i in 1:n_nodes if !haskey(old_to_new, Int(bus[i, 1]))]
        error("无法访问所有节点。未访问节点: $unvisited")
    end
    
    # 创建新的bus矩阵
    new_bus = copy(bus)
    for i in 1:n_nodes
        old_idx = Int(bus[i, 1])
        new_bus[i, 1] = old_to_new[old_idx]
    end
    
    # 创建新的branch矩阵，并确保小编号为from，大编号为to
    new_branch = copy(branch)
    for i in eachindex(branch[:,1])
        from_old = Int(branch[i, 1])
        to_old = Int(branch[i, 2])
        
        from_new = old_to_new[from_old]
        to_new = old_to_new[to_old]
        
        # 确保小编号为from，大编号为to
        if from_new < to_new
            new_branch[i, 1] = from_new
            new_branch[i, 2] = to_new
        else
            new_branch[i, 1] = to_new
            new_branch[i, 2] = from_new
        end
    end
    
    # 返回新的矩阵和节点映射
    return new_bus, new_branch, old_to_new
end

function identify_slack_nodes(busAC, busDC)
    # 找出交流侧的松弛节点（类型为3的节点）
    ac_slack_indices = findall(busAC[:, 2] .== 3)
    if isempty(ac_slack_indices)
        error("交流侧没有找到松弛节点（类型为3的节点）")
    elseif length(ac_slack_indices) > 1
        println("警告：交流侧有多个松弛节点，选择第一个")
    end
    ac_root = Int(busAC[ac_slack_indices[1], 1])
    
    # 找出直流侧的松弛节点（类型为2的节点）
    dc_slack_indices = findall(busDC[:, 2] .== 2)
    if isempty(dc_slack_indices)
        error("直流侧没有找到松弛节点（类型为2的节点）")
    elseif length(dc_slack_indices) > 1
        println("警告：直流侧有多个松弛节点，选择第一个")
    end
    dc_root = Int(busDC[dc_slack_indices[1], 1])
    
    return ac_root, dc_root
end

function convert_branch_matrix_to_tuples(branch_matrix)
    branches = []
    for i in eachindex(branch_matrix[:,1])
        from = Int(branch_matrix[i, 1])
        to = Int(branch_matrix[i, 2])
        push!(branches, (from, to))
    end
    return branches
end

function process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    # 识别松弛节点（根节点）
    ac_root, dc_root = identify_slack_nodes(busAC, busDC)
    
    # 对交流网络进行基于路径的重新编号，使根节点为1
    new_busAC, new_branchAC, ac_node_map = reindex_network_by_path(busAC, branchAC, ac_root)
    
    # 对直流网络进行基于路径的重新编号，使根节点为1
    new_busDC, new_branchDC, dc_node_map = reindex_network_by_path(busDC, branchDC, dc_root)
    
    # 对节点矩阵按节点编号排序
    ac_sort_idx = sortperm(new_busAC[:, 1])
    new_busAC = new_busAC[ac_sort_idx, :]
    
    dc_sort_idx = sortperm(new_busDC[:, 1])
    new_busDC = new_busDC[dc_sort_idx, :]
    
    # 处理节点编号
    Nac = size(new_busAC, 1)
    Ndc = size(new_busDC, 1)
    
    # 对支路矩阵按from节点编号从小到大排序
    ac_branch_sort_idx = sortperm(new_branchAC[:, 1])
    new_branchAC = new_branchAC[ac_branch_sort_idx, :]
    
    dc_branch_sort_idx = sortperm(new_branchDC[:, 1])
    new_branchDC = new_branchDC[dc_branch_sort_idx, :]
    
    # 转换排序后的支路数据为元组列表
    acbranches = convert_branch_matrix_to_tuples(new_branchAC)
    dcbranches = convert_branch_matrix_to_tuples(new_branchDC)
    
    # 创建支路索引映射
    acbranch_idx = Dict{Int, Int}()
    for (i, (from, to)) in enumerate(acbranches)
        acbranch_idx[to] = i
    end
    
    dcbranch_idx = Dict{Int, Int}()
    for (i, (from, to)) in enumerate(dcbranches)
        dcbranch_idx[to] = i
    end
    
    # 提取负载数据并转换为p.u.
    acP_load = new_busAC[:, 3] ./ baseMVA
    acQ_load = new_busAC[:, 4] ./ baseMVA
    dcP_load = new_busDC[:, 3] ./ baseMVA

    # 创建对角阻抗矩阵
    Rac_values = new_branchAC[:, 3]
    Xac_values = new_branchAC[:, 4]
    Rdc_values = new_branchDC[:, 3]
    
    # 创建对角矩阵
    Rac = Diagonal(Rac_values)
    Xac = Diagonal(Xac_values)
    Rdc = Diagonal(Rdc_values)
    
    # 返回处理后的数据
    return Dict(
        "Nac" => Nac,
        "Ndc" => Ndc,
        "acbranches" => acbranches,
        "dcbranches" => dcbranches,
        "acbranch_idx" => acbranch_idx,
        "dcbranch_idx" => dcbranch_idx,
        "acP_load" => acP_load,
        "acQ_load" => acQ_load,
        "dcP_load" => dcP_load,
        "busAC" => new_busAC,
        "busDC" => new_busDC,
        "branchAC" => new_branchAC,
        "branchDC" => new_branchDC,
        "ac_node_map" => ac_node_map,
        "dc_node_map" => dc_node_map,
        "Rac" => Rac,
        "Xac" => Xac,
        "Rdc" => Rdc
    )
end

function run_lindistflow_no_slack(mpc)
    # 创建mpc的深拷贝，避免修改原始数据
    mpc_result = deepcopy(mpc)
    
    baseMVA = mpc_result["baseMVA"]  # 基准功率
    # 网络参数
    busAC = mpc_result["busAC"]
    busDC = mpc_result["busDC"]
    branchAC = mpc_result["branchAC"]
    branchDC = mpc_result["branchDC"]
    C_inv = mpc_result["Ci"]
    C_rec = mpc_result["Cr"]
    v_min = 0.9
    v_max = 1.05
    P_inv0 = mpc_result["P_inv"] ./baseMVA
    Q_inv0 = mpc_result["Q_inv"]./baseMVA
    P_rec0 = mpc_result["P_inv_dc"] ./baseMVA
    
    # 添加调试信息
    println("交流节点数: ", size(busAC, 1))
    println("直流节点数: ", size(busDC, 1))
    println("交流支路数: ", size(branchAC, 1))
    println("直流支路数: ", size(branchDC, 1))
    println("逆变器数量: ", length(C_inv))
    println("整流器数量: ", length(C_rec))
    
    network_data = process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    
    # 处理数据
    Nac = network_data["Nac"]
    Ndc = network_data["Ndc"]
    acbranches = network_data["acbranches"]
    dcbranches = network_data["dcbranches"]
    acP_load = network_data["acP_load"]
    acQ_load = network_data["acQ_load"]
    dcP_load = network_data["dcP_load"]
    acbranch_idx = network_data["acbranch_idx"]
    dcbranch_idx = network_data["dcbranch_idx"]
    Rac = network_data["Rac"]
    Xac = network_data["Xac"]
    Rdc = network_data["Rdc"]
    ac_node_map = network_data["ac_node_map"]
    dc_node_map = network_data["dc_node_map"]
    
    # 添加调试信息
    println("\n重新编号后:")
    println("交流节点数 (Nac): ", Nac)
    println("直流节点数 (Ndc): ", Ndc)
    println("交流支路数: ", length(acbranches))
    println("直流支路数: ", length(dcbranches))
    println("交流负载数: ", length(acP_load))
    println("直流负载数: ", length(dcP_load))
    
    # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
    C_inv_new = [ac_node_map[node] for node in C_inv]
    C_rec_new = [dc_node_map[node] for node in C_rec]
    
    println("\n更新后的逆变器位置: ", C_inv_new)
    println("更新后的整流器位置: ", C_rec_new)
    
    # 创建交流支路节点关联矩阵
    Aac = build_incidence_matrix(Nac, acbranches)
    println("\nAac维度: ", size(Aac))

    # 创建直流支路节点关联矩阵
    Adc = build_incidence_matrix(Ndc, dcbranches)
    println("Adc维度: ", size(Adc))
    
    # 检查acbranch_idx映射
    println("\nacbranch_idx映射:")
    for (k, v) in acbranch_idx
        println("节点 $k -> 支路索引 $v")
    end
    
    # 检查dcbranch_idx映射
    println("\ndcbranch_idx映射:")
    for (k, v) in dcbranch_idx
        println("节点 $k -> 支路索引 $v")
    end
    
    # 检查是否所有非根节点都有对应的支路索引
    missing_ac_nodes = [i for i in 2:Nac if !haskey(acbranch_idx, i)]
    if !isempty(missing_ac_nodes)
        println("\n警告: 以下交流节点没有对应的支路索引: ", missing_ac_nodes)
    end
    
    missing_dc_nodes = [i for i in 2:Ndc if !haskey(dcbranch_idx, i)]
    if !isempty(missing_dc_nodes)
        println("\n警告: 以下直流节点没有对应的支路索引: ", missing_dc_nodes)
    end

    # 使用无松弛变量的模型
    model, Pac, Qac, Pdc, vac, vdc, P_inv, Q_inv, P_rec = create_distflow_model_no_slack(
        Nac, Ndc, C_inv_new, C_rec_new, P_inv0, Q_inv0, P_rec0, acP_load, acQ_load, 
        dcP_load, v_min, v_max, acbranch_idx, dcbranch_idx, acbranches, dcbranches,
        Aac, Adc, Rac, Xac, Rdc
    )
   
    # 求解模型
    optimize!(model)
    
    # 创建新旧节点编号的反向映射
    ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
    dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
    
    # 输出结果并更新mpc_result
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
        println("\n===== 优化成功! =====")
        
        # 输出目标函数值
        obj_value = objective_value(model)
        println("\n目标函数值 (换流器功率最大偏差): ", obj_value, " p.u. (", obj_value * baseMVA * 1000, " kW)")
        
        # 设置更高的精度显示
        digits_pu = 12  # 标幺值显示的小数位数
        digits_kw = 8   # kW/kvar显示的小数位数
        
        # 验证功率平衡
        println("\n========== 功率平衡验证 ==========")
        
        # 交流有功功率平衡验证
        println("\n----- 交流有功功率平衡验证 -----")
        println("节点ID    不平衡量(p.u.)                不平衡量(kW)")
        
        # 获取变量值
        branch_indices_ac = [get(acbranch_idx, k, 0) for k in 2:Nac]
        acP_n_val = value.(Pac[branch_indices_ac])
        P_inverter_val = zeros(Nac-1)
        
        for i in 1:length(C_inv)
            if C_inv_new[i] > 1  # 确保不是根节点
                node_idx = C_inv_new[i] - 1  # 调整索引（去掉根节点）
                P_inverter_val[node_idx] = value(P_inv[i])
            end
        end
        
        # 计算功率不平衡量
        acP_mismatch = Aac' * acP_n_val + acP_load[2:end] - P_inverter_val
        
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]  # 跳过根节点
            p_mismatch = acP_mismatch[i]
            p_mismatch_kw = p_mismatch * baseMVA * 1000  # 转换为kW
            
            println(lpad(string(old_node), 6), "    ", 
                    lpad(string(round(p_mismatch, digits=digits_pu)), 24), "    ", 
                    lpad(string(round(p_mismatch_kw, digits=digits_kw)), 16))
        end
        
        # 更新mpc中的换流器功率
        mpc_result["P_inv"] = value.(P_inv) .* baseMVA
        mpc_result["Q_inv"] = value.(Q_inv) .* baseMVA
        mpc_result["P_inv_dc"] = value.(P_rec) .* baseMVA
        
        # 更新mpc_result中的busAC和busDC
        for i in 1:Nac
            old_node = ac_new_to_old[i]
            voltage_value = sqrt(value(vac[i]))
            
            # 查找节点在busAC中的位置并更新
            idx = findfirst(x -> x == old_node, busAC[:, 1])
            if idx !== nothing
                busAC[idx, 8] = voltage_value
            end
        end
        
        for i in 1:Ndc
            old_node = dc_new_to_old[i]
            voltage_value = sqrt(value(vdc[i]))
            
            # 查找节点在busDC中的位置并更新
            idx = findfirst(x -> x == old_node, busDC[:, 1])
            if idx !== nothing
                busDC[idx, 8] = voltage_value
            end
        end
        
        mpc_result["busAC"] = busAC
        mpc_result["busDC"] = busDC
    else
        println("\n===== 优化失败! =====")
        println("终止状态: ", termination_status(model))
        println("终止原因: ", raw_status(model))
    end
    
    # 添加优化状态标志
    mpc_result["success"] = (termination_status(model) == MOI.OPTIMAL || 
                            termination_status(model) == MOI.LOCALLY_SOLVED)
    
    return mpc_result
end
lindistflow_result=run_lindistflow_no_slack(mpc)