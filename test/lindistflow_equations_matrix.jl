"""
This file is used to solve the linear AC/DC hybrid distribution flow model (LinDistFlow).
"""

using LinearAlgebra

# 定义线性化交直流混合配电网络的潮流计算函数
function solve_linear_acdc_power_flow(
    Nac::Int,                  # 交流节点数量
    Ndc::Int,                  # 直流节点数量
    C_inv::Vector{Int},        # 逆变器位置
    C_rec::Vector{Int},        # 整流器位置
    P_inv::Vector{Float64},    # 逆变器有功功率
    Q_inv::Vector{Float64},    # 逆变器无功功率
    P_rec::Vector{Float64},    # 整流器有功功率
    acP_load::Vector{Float64}, # 交流负载有功功率
    acQ_load::Vector{Float64}, # 交流负载无功功率
    dcP_load::Vector{Float64}, # 直流负载有功功率
    acbranches,                # 交流支路
    dcbranches,                # 直流支路
    Rac,                       # 交流支路电阻
    Xac,                       # 交流支路电抗
    Rdc                        # 直流支路电阻
)
    # 创建节点到支路的映射
    acbranch_idx = Dict{Int, Int}()
    for (i, (from, to)) in enumerate(acbranches)
        acbranch_idx[to] = i
    end
    
    dcbranch_idx = Dict{Int, Int}()
    for (i, (from, to)) in enumerate(dcbranches)
        dcbranch_idx[to] = i
    end
    
    # 创建支路节点关联矩阵
    Aac = build_incidence_matrix(Nac, acbranches)
    Adc = build_incidence_matrix(Ndc, dcbranches)
    
    # 提取非根节点的负载
    acP_load_nonroot = acP_load[2:end]
    acQ_load_nonroot = acQ_load[2:end]
    dcP_load_nonroot = dcP_load[2:end]
    
    # 创建逆变器和整流器功率向量
    P_inverter = zeros(Nac-1)
    Q_inverter = zeros(Nac-1)
    P_rectifier = zeros(Ndc-1)
    
    for i in eachindex(C_inv)
        if C_inv[i] > 1  # 非根节点
            P_inverter[C_inv[i]-1] = P_inv[i]
            Q_inverter[C_inv[i]-1] = Q_inv[i]
        end
    end
    
    for i in eachindex(C_rec)
        if C_rec[i] > 1  # 非根节点
            P_rectifier[C_rec[i]-1] = P_rec[i]
        end
    end
    
    # 构建线性方程组
    # 1. 交流支路功率方程
    n_acbranches = length(acbranches)
    n_dcbranches = length(dcbranches)
    
    # 方程总数: 交流支路功率(P,Q) + 直流支路功率 + 交流节点电压 + 直流节点电压
    n_eqs = 2*n_acbranches + n_dcbranches + (Nac-1) + (Ndc-1)
    
    # 变量总数: 交流支路功率(P,Q) + 直流支路功率 + 交流节点电压 + 直流节点电压
    n_vars = 2*n_acbranches + n_dcbranches + (Nac-1) + (Ndc-1)
    
    # 构建系数矩阵和常数向量
    A = zeros(n_eqs, n_vars)
    b = zeros(n_eqs)
    
    # 变量索引
    idx_Pac = 1:n_acbranches
    idx_Qac = (n_acbranches+1):(2*n_acbranches)
    idx_Pdc = (2*n_acbranches+1):(2*n_acbranches+n_dcbranches)
    idx_vac = (2*n_acbranches+n_dcbranches+1):(2*n_acbranches+n_dcbranches+Nac-1)
    idx_vdc = (2*n_acbranches+n_dcbranches+Nac):(2*n_acbranches+n_dcbranches+Nac-1+Ndc-1)
    
    # 方程索引
    eq_idx = 1
    
    # 1. 交流节点功率平衡方程 (P)
    for i in 1:(Nac-1)
        # 支路功率平衡
        A[eq_idx, idx_Pac] = Aac[i, :]
        
        # 右侧常数: 负载 - 逆变器功率
        b[eq_idx] = -acP_load_nonroot[i] + P_inverter[i]
        
        eq_idx += 1
    end
    
    # 2. 交流节点功率平衡方程 (Q)
    for i in 1:(Nac-1)
        # 支路功率平衡
        A[eq_idx, idx_Qac] = Aac[i, :]
        
        # 右侧常数: 负载 - 逆变器功率
        b[eq_idx] = -acQ_load_nonroot[i] + Q_inverter[i]
        
        eq_idx += 1
    end
    
    # 3. 直流节点功率平衡方程
    for i in 1:(Ndc-1)
        # 支路功率平衡
        A[eq_idx, idx_Pdc] = Adc[i, :]
        
        # 右侧常数: 负载 + 整流器功率
        b[eq_idx] = -dcP_load_nonroot[i] - P_rectifier[i]
        
        eq_idx += 1
    end
    
    # 4. 交流节点电压方程
    for i in 2:Nac
        branch_idx = acbranch_idx[i]
        
        # 电压降方程: v_to = v_from - 2(r*P + x*Q)
        # 对应的支路功率
        A[eq_idx, idx_Pac[branch_idx]] = 2.0 * Rac[branch_idx, branch_idx]
        A[eq_idx, idx_Qac[branch_idx]] = 2.0 * Xac[branch_idx, branch_idx]
        
        # 对应的节点电压
        from_node = acbranches[branch_idx][1]
        if from_node == 1
            # 根节点电压为1.0
            b[eq_idx] += 1.0
        else
            # 非根节点电压为变量
            A[eq_idx, idx_vac[from_node-1]] = -1.0
        end
        
        # 目标节点电压
        A[eq_idx, idx_vac[i-1]] = 1.0
        
        eq_idx += 1
    end
    
    # 5. 直流节点电压方程
    for i in 2:Ndc
        branch_idx = dcbranch_idx[i]
        
        # 电压降方程: v_to = v_from - 2*r*P
        # 对应的支路功率
        A[eq_idx, idx_Pdc[branch_idx]] = 2.0 * Rdc[branch_idx, branch_idx]
        
        # 对应的节点电压
        from_node = dcbranches[branch_idx][1]
        if from_node == 1
            # 根节点电压为1.0
            b[eq_idx] += 1.0
        else
            # 非根节点电压为变量
            A[eq_idx, idx_vdc[from_node-1]] = -1.0
        end
        
        # 目标节点电压
        A[eq_idx, idx_vdc[i-1]] = 1.0
        
        eq_idx += 1
    end
    
    # 求解线性方程组
    x = A \ b
    
    # 提取结果
    Pac = x[idx_Pac]
    Qac = x[idx_Qac]
    Pdc = x[idx_Pdc]
    vac = ones(Nac)  # 根节点电压为1.0
    vac[2:end] = x[idx_vac]
    vdc = ones(Ndc)  # 根节点电压为1.0
    vdc[2:end] = x[idx_vdc]
    
    return vac, vdc, Pac, Qac, Pdc
end

# 构建网络关联矩阵
function build_incidence_matrix(n_nodes, branches)
    # 创建关联矩阵
    A = zeros(n_nodes-1, length(branches))
    
    for (i, (from, to)) in enumerate(branches)
        if from != 1  # 非根节点
            A[from-1, i] = 1  # 从节点流出为正
        end
        if to != 1    # 非根节点
            A[to-1, i] = -1   # 到节点流入为负
        end
    end
    
    return A
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

function run_linear_power_flow()
    baseMVA = 100.0  # 基准功率
    # 网络参数
    busAC = [1.0 3.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 10.0 1.0 1.05 0.8;
            2.0 1.0 0.1275 0.07902 0.0 0.0 0.0 1.0 0.0 10.0 1.0 1.05 0.8]
    busDC = [1.0 1.0 0.150 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8;
            2.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8;
            3.0 2.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.824 1.0 1.05 0.8]
    branchAC = [1.0  2.0  5.0   5.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0]
    branchDC = [1.0  2.0  5.8912  0.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0;
                2.0  3.0  3.0193  0.0  0.0  100.0  100.0  100.0  0.0  0.0  1.0  -180.0  180.0  0.0]
    C_inv = [2]
    C_rec = [1]
    P_inv = [0.1]./baseMVA
    Q_inv = [0.05]./baseMVA
    P_rec = [0.1111]./baseMVA
    
    network_data = process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    
    # 处理数据
    Nac = network_data["Nac"]
    Ndc = network_data["Ndc"]
    acbranches = network_data["acbranches"]
    dcbranches = network_data["dcbranches"]
    acP_load = network_data["acP_load"]
    acQ_load = network_data["acQ_load"]
    dcP_load = network_data["dcP_load"]
    Rac = network_data["Rac"]
    Xac = network_data["Xac"]
    Rdc = network_data["Rdc"]
    ac_node_map = network_data["ac_node_map"]
    dc_node_map = network_data["dc_node_map"]
    
    # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
    C_inv_new = [ac_node_map[node] for node in C_inv]
    C_rec_new = [dc_node_map[node] for node in C_rec]
    
    # 使用线性方法求解潮流
    vac, vdc, Pac, Qac, Pdc = solve_linear_acdc_power_flow(
        Nac, Ndc, C_inv_new, C_rec_new, P_inv, Q_inv, P_rec,
        acP_load, acQ_load, dcP_load,
        acbranches, dcbranches,
        Rac, Xac, Rdc
    )
    
    # 创建新旧节点编号的反向映射
    ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
    dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
    
    # 输出结果
    println("\n线性潮流计算结果:")
    
    # 输出交流节点电压（使用原始编号）
    println("\n交流节点电压:")
    for i in 1:Nac
        old_node = ac_new_to_old[i]
        println("节点 $old_node: v = $(sqrt(vac[i]))")
    end
    
    # 输出直流节点电压（使用原始编号）
    println("\n直流节点电压:")
    for i in 1:Ndc
        old_node = dc_new_to_old[i]
        println("节点 $old_node: v = $(sqrt(vdc[i]))")
    end
    
    # 输出交流支路功率（使用原始编号）
    println("\n交流支路功率:")
    for i in eachindex(acbranches)
        from_new, to_new = acbranches[i]
        from_old = ac_new_to_old[from_new]
        to_old = ac_new_to_old[to_new]
        println("支路 ($from_old, $to_old): P = $(Pac[i]), Q = $(Qac[i])")
    end
    
    # 输出直流支路功率（使用原始编号）
    println("\n直流支路功率:")
    for i in eachindex(dcbranches)
        from_new, to_new = dcbranches[i]
        from_old = dc_new_to_old[from_new]
        to_old = dc_new_to_old[to_new]
        println("支路 ($from_old, $to_old): P = $(Pdc[i])")
    end
    
    return vac, vdc, Pac, Qac, Pdc
end

# 运行线性潮流计算
run_linear_power_flow()
