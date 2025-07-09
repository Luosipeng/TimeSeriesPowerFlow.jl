"""
This file is used to test the linear AC/DC hybrid distribution flow model.
"""

using JuMP
using Gurobi  # 替换Ipopt为Gurobi
using LinearAlgebra

# 定义配电网络参数
function create_distflow_model(
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
    
    # 设置Gurobi参数
    set_optimizer_attribute(model, "OutputFlag", 1)       # 显示求解过程
    set_optimizer_attribute(model, "MIPGap", 1e-4)        # MIP相对间隙
    set_optimizer_attribute(model, "TimeLimit", 300)      # 求解时间限制（秒）
    set_optimizer_attribute(model, "Threads", 0)          # 使用所有可用线程
    set_optimizer_attribute(model, "NumericFocus", 3)     # 提高数值稳定性
     
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
    
    # 添加松弛变量，用于功率平衡约束
    slack_bound = 1e-6  # 松弛变量的界限，适当放宽
    @variable(model, -slack_bound <= p_slack[1:Nac-1] <= slack_bound)
    @variable(model, -slack_bound <= q_slack[1:Nac-1] <= slack_bound)
    @variable(model, -slack_bound <= dc_slack[1:Ndc-1] <= slack_bound)
    
    # 为非根节点设置电压约束
    for i in 2:Nac
        @constraint(model, v_min^2 <= vac[i] <= v_max^2)
    end
    
    for i in 2:Ndc
        @constraint(model, v_min^2 <= vdc[i] <= v_max^2)
    end
    
    # 固定根节点电压为1.0的平方
    # @constraint(model, vac[1] == 1.0)
    # @constraint(model, vdc[1] == 1.0)

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
    
    @constraint(model, max_dev >= vac[1] - 1.0)  # 确保最大偏差非负
    @constraint(model, max_dev >= 1.0 - vac[1] )  # 确保最大偏差非负

    @constraint(model, max_dev >= vdc[1] - 1.0)  # 确保最大偏差非负
    @constraint(model, max_dev >= 1.0 - vdc[1] )  # 确保最大偏差非负

    # 最小化最大偏差
    @objective(model, Min, max_dev)
    
    # 创建逆变器功率向量（使用JuMP表达式，而不是Float64向量）
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

    # 创建整流器功率向量（使用JuMP表达式，而不是Float64向量）
    P_rectifier = zeros(AffExpr, Ndc-1)
    
    # 将整流器功率映射到对应节点
    for i in 1:length(C_rec)
        if C_rec[i] > 1  # 确保不是根节点
            node_idx = C_rec[i] - 1  # 调整索引（去掉根节点）
            P_rectifier[node_idx] = P_rec[i]
        end
    end
    
    # 交流节点功率平衡方程（添加松弛变量）
    acP_load = acP_load[2:end]  # 去掉根节点
    acP_n = Pac[map(k -> acbranch_idx[k], 2:Nac)]
   
    acQ_load = acQ_load[2:end]  # 去掉根节点
    acQ_n = Qac[map(k -> acbranch_idx[k], 2:Nac)]

    @constraint(model, Aac' * acP_n + acP_load - P_inverter + p_slack .== 0)
    @constraint(model, Aac' * acQ_n + acQ_load - Q_inverter + q_slack .== 0)

    # 直流节点功率平衡方程（添加松弛变量）
    dcP_load = dcP_load[2:end]  # 去掉根节点
    dcP_n = Pdc[map(k -> dcbranch_idx[k], 2:Ndc)]

    @constraint(model, Adc' * dcP_n + dcP_load + P_rectifier + dc_slack .== 0)

    # 交流节点电压平衡方程
    @constraint(model, -2.0 * Rac * acP_n - 2.0 * Xac * acQ_n + vac[1:end-1] .== vac[2:end])

    # 直流节点电压平衡方程
    @constraint(model, -2.0 * Rdc * dcP_n + vdc[1:end-1] .== vdc[2:end])
    
    # 改进换流器效率模型 - 使用二元变量
    @variable(model, z[1:length(P_inv)], Bin)  # 使用二元变量而非连续变量
    
    # 估计P_inv的上下界
    P_inv_bound = 10.0  # 根据实际问题调整这个值
    
    for i in 1:length(P_inv)
        # 限制P_inv的范围
        @constraint(model, -P_inv_bound <= P_inv[i] <= P_inv_bound)
        
        # 约束z的值基于P_inv的符号
        @constraint(model, P_inv[i] <= P_inv_bound * z[i])           # 当P_inv > 0时，z必须为1
        @constraint(model, P_inv[i] >= -P_inv_bound * (1 - z[i]))    # 当P_inv < 0时，z必须为0
        
        # 根据z的值确定P_rec的计算方式
        @constraint(model, P_rec[i] <= P_inv[i]/0.9 + P_inv_bound*(1-z[i]))
        @constraint(model, P_rec[i] >= P_inv[i]/0.9 - P_inv_bound*(1-z[i]))
        
        @constraint(model, P_rec[i] <= P_inv[i]*0.9 + P_inv_bound*z[i])
        @constraint(model, P_rec[i] >= P_inv[i]*0.9 - P_inv_bound*z[i])
    end
        
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

# function run_fixed_gen_test_case()
#     baseMVA = 100.0  # 基准功率
#     # 网络参数
#     busAC = [1.0 3.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 10.0 1.0 1.05 0.8;
#             2.0 1.0 0.1275 0.07902 0.0 0.0 0.0 1.0 0.0 10.0 1.0 1.05 0.8]
#     busDC = [1.0 1.0 0.150 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8;
#             2.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8;
#             3.0 2.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8]
#     branchAC = [1.0  2.0  5.0   5.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0]
#     branchDC = [1.0  2.0  5.8912  0.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0;
#                 2.0  3.0  3.0193  0.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0]
#     C_inv = [2]
#     C_rec = [1]
#     v_min = 0.8
#     v_max = 1.05
#     P_inv0 = [0.1]./baseMVA
#     Q_inv0 = [0.05]./baseMVA
#     P_rec0 = [0.1111]./baseMVA
    
#     network_data = process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    
#     # 处理数据
#     Nac = network_data["Nac"]
#     Ndc = network_data["Ndc"]
#     acbranches = network_data["acbranches"]
#     dcbranches = network_data["dcbranches"]
#     acP_load = network_data["acP_load"]
#     acQ_load = network_data["acQ_load"]
#     dcP_load = network_data["dcP_load"]
#     acbranch_idx = network_data["acbranch_idx"]
#     dcbranch_idx = network_data["dcbranch_idx"]
#     Rac = network_data["Rac"]
#     Xac = network_data["Xac"]
#     Rdc = network_data["Rdc"]
#     ac_node_map = network_data["ac_node_map"]
#     dc_node_map = network_data["dc_node_map"]
    
#     # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
#     C_inv_new = [ac_node_map[node] for node in C_inv]
#     C_rec_new = [dc_node_map[node] for node in C_rec]
    
#     # 创建交流支路节点关联矩阵
#     Aac = build_incidence_matrix(Nac, acbranches)

#     # 创建直流支路节点关联矩阵
#     Adc = build_incidence_matrix(Ndc, dcbranches)

#     model, Pac, Qac, Pdc, vac, vdc, P_inv, Q_inv, P_rec = create_distflow_model(
#         Nac, Ndc, C_inv_new, C_rec_new, P_inv0, Q_inv0, P_rec0, acP_load, acQ_load, 
#         dcP_load, v_min, v_max, acbranch_idx, dcbranch_idx, acbranches, dcbranches,
#         Aac, Adc, Rac, Xac, Rdc
#     )
   
#     # 求解模型
#     optimize!(model)
    
#     # 输出结果
#     if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
#         println("\n优化成功!")
#         println("目标函数值 (最大功率失配): ", objective_value(model))
        
#         # 创建新旧节点编号的反向映射
#         ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
#         dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
        
#         # 输出交流支路功率（使用原始编号）
#         println("\n交流支路功率:")
#         for i in eachindex(acbranches)
#             from_new, to_new = acbranches[i]
#             from_old = ac_new_to_old[from_new]
#             to_old = ac_new_to_old[to_new]
#             println("支路 ($from_old, $to_old): P = $(value(Pac[i])), Q = $(value(Qac[i]))")
#         end
        
#         # 输出直流支路功率（使用原始编号）
#         println("\n直流支路功率:")
#         for i in eachindex(dcbranches)
#             from_new, to_new = dcbranches[i]
#             from_old = dc_new_to_old[from_new]
#             to_old = dc_new_to_old[to_new]
#             println("支路 ($from_old, $to_old): P = $(value(Pdc[i]))")
#         end
        
#         # 输出交流节点电压（使用原始编号）
#         println("\n交流节点电压:")
#         for i in 1:Nac
#             old_node = ac_new_to_old[i]
#             println("节点 $old_node: v = $(sqrt(value(vac[i])))")
#         end
        
#         # 输出直流节点电压（使用原始编号）
#         println("\n直流节点电压:")
#         for i in 1:Ndc
#             old_node = dc_new_to_old[i]
#             println("节点 $old_node: v = $(sqrt(value(vdc[i])))")
#         end
        
#         # 输出逆变器和整流器功率
#         println("\n逆变器功率:")
#         for i in eachindex(C_inv)
#             old_node = C_inv[i]
#             println("节点 $old_node: P = $(value(P_inv[i])), Q = $(value(Q_inv[i]))")
#         end
        
#         println("\n整流器功率:")
#         for i in eachindex(C_rec)
#             old_node = C_rec[i]
#             println("节点 $old_node: P = $(value(P_rec[i]))")
#         end
#     else
#         println("优化失败: ", termination_status(model))
#     end
    
#     return model, Pac, Qac, Pdc, vac, vdc, P_inv, Q_inv, P_rec
# end

function run_lindistflow(mpc)
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
    
    network_data = PowerFlow.process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    
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
    
    # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
    C_inv_new = [ac_node_map[node] for node in C_inv]
    C_rec_new = [dc_node_map[node] for node in C_rec]
    
    # 创建交流支路节点关联矩阵
    Aac = PowerFlow.build_incidence_matrix(Nac, acbranches)

    # 创建直流支路节点关联矩阵
    Adc = PowerFlow.build_incidence_matrix(Ndc, dcbranches)

    model, Pac, Qac, Pdc, vac, vdc, P_inv, Q_inv, P_rec = PowerFlow.create_distflow_model(
        Nac, Ndc, C_inv_new, C_rec_new, P_inv0, Q_inv0, P_rec0, acP_load, acQ_load, 
        dcP_load, v_min, v_max, acbranch_idx, dcbranch_idx, acbranches, dcbranches,
        Aac, Adc, Rac, Xac, Rdc
    )
   
    # 获取模型中的松弛变量
    p_slack_var = model[:p_slack]
    q_slack_var = model[:q_slack]
    dc_slack_var = model[:dc_slack]
    
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
        
        # 输出松弛变量（功率失配）
        println("\n========== 功率失配信息 (松弛变量) ==========")
        
        # 交流有功功率失配
        println("\n----- 交流有功功率失配 (p_slack) -----")
        println("节点ID    失配值(p.u.)                失配值(kW)")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]  # 跳过根节点
            p_slack_val = value(p_slack_var[i])
            p_slack_kw = p_slack_val * baseMVA * 1000  # 转换为kW
            
            println(lpad(string(old_node), 6), "    ", 
                    lpad(string(round(p_slack_val, digits=digits_pu)), 24), "    ", 
                    lpad(string(round(p_slack_kw, digits=digits_kw)), 16))
        end
        
        # 交流无功功率失配
        println("\n----- 交流无功功率失配 (q_slack) -----")
        println("节点ID    失配值(p.u.)                失配值(kvar)")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]  # 跳过根节点
            q_slack_val = value(q_slack_var[i])
            q_slack_kvar = q_slack_val * baseMVA * 1000  # 转换为kvar
            
            println(lpad(string(old_node), 6), "    ", 
                    lpad(string(round(q_slack_val, digits=digits_pu)), 24), "    ", 
                    lpad(string(round(q_slack_kvar, digits=digits_kw)), 16))
        end
        
        # 直流功率失配
        println("\n----- 直流功率失配 (dc_slack) -----")
        println("节点ID    失配值(p.u.)                失配值(kW)")
        for i in 1:Ndc-1
            old_node = dc_new_to_old[i+1]  # 跳过根节点
            dc_slack_val = value(dc_slack_var[i])
            dc_slack_kw = dc_slack_val * baseMVA * 1000  # 转换为kW
            
            println(lpad(string(old_node), 6), "    ", 
                    lpad(string(round(dc_slack_val, digits=digits_pu)), 24), "    ", 
                    lpad(string(round(dc_slack_kw, digits=digits_kw)), 16))
        end
        
        # 计算总功率失配
        total_p_slack = sum(abs.(value.(p_slack_var))) * baseMVA * 1000  # kW
        total_q_slack = sum(abs.(value.(q_slack_var))) * baseMVA * 1000  # kvar
        total_dc_slack = sum(abs.(value.(dc_slack_var))) * baseMVA * 1000  # kW
        
        println("\n----- 总功率失配 -----")
        println("交流有功总失配: ", round(total_p_slack, digits=digits_kw), " kW")
        println("交流无功总失配: ", round(total_q_slack, digits=digits_kw), " kvar")
        println("直流有功总失配: ", round(total_dc_slack, digits=digits_kw), " kW")
        
        # 输出最大失配值
        max_p_slack = maximum(abs.(value.(p_slack_var))) * baseMVA * 1000  # kW
        max_q_slack = maximum(abs.(value.(q_slack_var))) * baseMVA * 1000  # kvar
        max_dc_slack = maximum(abs.(value.(dc_slack_var))) * baseMVA * 1000  # kW
        
        println("\n----- 最大功率失配 -----")
        println("交流有功最大失配: ", round(max_p_slack, digits=digits_kw), " kW")
        println("交流无功最大失配: ", round(max_q_slack, digits=digits_kw), " kvar")
        println("直流有功最大失配: ", round(max_dc_slack, digits=digits_kw), " kW")
        
        # 输出原始松弛变量值（不四舍五入）
        println("\n========== 原始松弛变量值 (完整精度) ==========")
        
        println("\n----- 交流有功功率失配 (p_slack) -----")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]
            p_slack_val = value(p_slack_var[i])
            p_slack_kw = p_slack_val * baseMVA * 1000
            println("节点 $old_node: ", p_slack_val, " p.u. (", p_slack_kw, " kW)")
        end
        
        println("\n----- 交流无功功率失配 (q_slack) -----")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]
            q_slack_val = value(q_slack_var[i])
            q_slack_kvar = q_slack_val * baseMVA * 1000
            println("节点 $old_node: ", q_slack_val, " p.u. (", q_slack_kvar, " kvar)")
        end
        
        println("\n----- 直流功率失配 (dc_slack) -----")
        for i in 1:Ndc-1
            old_node = dc_new_to_old[i+1]
            dc_slack_val = value(dc_slack_var[i])
            dc_slack_kw = dc_slack_val * baseMVA * 1000
            println("节点 $old_node: ", dc_slack_val, " p.u. (", dc_slack_kw, " kW)")
        end
        
        # 输出节点电压
        println("\n========== 节点电压 ==========")
        
        # 交流节点电压
        println("\n----- 交流节点电压 -----")
        println("节点ID   电压(p.u.)")
        for i in 1:Nac
            old_node = ac_new_to_old[i]
            voltage_value = sqrt(value(vac[i]))
            
            println(lpad(string(old_node), 6), "   ", 
                    lpad(string(round(voltage_value, digits=6)), 12))
            
            # 查找节点在busAC中的位置并更新
            idx = findfirst(x -> x == old_node, busAC[:, 1])
            if idx !== nothing
                busAC[idx, 8] = voltage_value
            end
        end
        
        # 直流节点电压
        println("\n----- 直流节点电压 -----")
        println("节点ID   电压(p.u.)")
        for i in 1:Ndc
            old_node = dc_new_to_old[i]
            voltage_value = sqrt(value(vdc[i]))
            
            println(lpad(string(old_node), 6), "   ", 
                    lpad(string(round(voltage_value, digits=6)), 12))
            
            # 查找节点在busDC中的位置并更新
            idx = findfirst(x -> x == old_node, busDC[:, 1])
            if idx !== nothing
                busDC[idx, 8] = voltage_value
            end
        end
        
        
    else
        println("\n===== 优化失败! =====")
        println("终止状态: ", termination_status(model))
        println("终止原因: ", raw_status(model))
        
        # 检查是否有解决方案
        if has_values(model)
            println("\n虽然未达到最优，但有解决方案。输出松弛变量:")
            
            # 输出目标函数值
            obj_value = objective_value(model)
            println("\n目标函数值 (换流器功率偏差): ", obj_value, " p.u. (", obj_value * baseMVA * 1000, " kW)")
            
            # 输出松弛变量（完整精度）
            println("\n----- 交流有功功率失配 (p_slack) -----")
            for i in 1:Nac-1
                old_node = ac_new_to_old[i+1]
                p_slack_val = value(p_slack_var[i])
                p_slack_kw = p_slack_val * baseMVA * 1000
                println("节点 $old_node: ", p_slack_val, " p.u. (", p_slack_kw, " kW)")
            end
            
            println("\n----- 交流无功功率失配 (q_slack) -----")
            for i in 1:Nac-1
                old_node = ac_new_to_old[i+1]
                q_slack_val = value(q_slack_var[i])
                q_slack_kvar = q_slack_val * baseMVA * 1000
                println("节点 $old_node: ", q_slack_val, " p.u. (", q_slack_kvar, " kvar)")
            end
            
            println("\n----- 直流功率失配 (dc_slack) -----")
            for i in 1:Ndc-1
                old_node = dc_new_to_old[i+1]
                dc_slack_val = value(dc_slack_var[i])
                dc_slack_kw = dc_slack_val * baseMVA * 1000
                println("节点 $old_node: ", dc_slack_val, " p.u. (", dc_slack_kw, " kW)")
            end
            
            # 保存松弛变量
            mpc_result["p_slack"] = value.(p_slack_var)
            mpc_result["q_slack"] = value.(q_slack_var)
            mpc_result["dc_slack"] = value.(dc_slack_var)
            mpc_result["p_slack_kw"] = value.(p_slack_var) .* baseMVA .* 1000
            mpc_result["q_slack_kvar"] = value.(q_slack_var) .* baseMVA .* 1000
            mpc_result["dc_slack_kw"] = value.(dc_slack_var) .* baseMVA .* 1000
        else
            println("没有可行解")
        end
    end
    # 更新mpc中的换流器功率
    mpc_result["P_inv"] = value.(P_inv) .* baseMVA
    mpc_result["Q_inv"] = value.(Q_inv) .* baseMVA
    mpc_result["P_inv_dc"] = value.(P_rec) .* baseMVA
    # 更新mpc_result中的busAC和busDC
    mpc_result["busAC"] = busAC
    mpc_result["busDC"] = busDC
    
    # 添加优化状态标志
    mpc_result["success"] = (termination_status(model) == MOI.OPTIMAL || 
                            termination_status(model) == MOI.LOCALLY_SOLVED)
    
    return mpc_result
end
