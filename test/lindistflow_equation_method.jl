"""
This file implements a direct equation-solving approach for AC/DC hybrid distribution flow model.
"""

using LinearAlgebra
using Plots

# 定义配电网络参数
function solve_distflow_equations(
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
    v_min::Float64,            # 最小电压限制
    v_max::Float64,            # 最大电压限制
    p_gen_fixed::Vector{Float64}, # 节点固定有功发电量
    q_gen_fixed::Vector{Float64}, # 节点固定无功发电量
    childrenAC,
    childrenDC,
    parentsAC,
    parentsDC,
    acbranch_idx,
    dcbranch_idx,
    acbranches,
    dcbranches,
    branchAC_new,  # 交流支路参数矩阵
    branchDC_new   # 直流支路参数矩阵
)
    # 初始化变量
    # 支路功率
    Pac = zeros(length(acbranches))
    Qac = zeros(length(acbranches))
    Pdc = zeros(length(dcbranches))
    
    # 节点电压
    vac = ones(Nac)  # 初始化所有交流节点电压为1.0 p.u.
    vdc = ones(Ndc)  # 初始化所有直流节点电压为1.0 p.u.
    
    # 正向扫描：计算支路功率
    # 交流网络正向扫描
    for n in Nac:-1:2  # 从叶节点向根节点
        parent = parentsAC[n]
        branch_idx = 0
        
        # 查找连接到父节点的支路索引
        if haskey(acbranch_idx, (parent, n))
            branch_idx = acbranch_idx[(parent, n)]
        else
            error("找不到连接节点 $parent 和 $n 的支路")
        end
        
        # 计算流入节点n的功率
        P_in = acP_load[n]
        Q_in = acQ_load[n]
        
        # 减去节点n的发电量
        if n in C_inv
            inv_index = findfirst(x -> x == n, C_inv)
            P_in -= P_inv[inv_index]
            Q_in -= Q_inv[inv_index]
        end
        P_in -= p_gen_fixed[n]
        Q_in -= q_gen_fixed[n]
        
        # 加上流向所有子节点的功率
        if haskey(childrenAC, n) && !isempty(childrenAC[n])
            for child in childrenAC[n]
                child_branch_idx = acbranch_idx[(n, child)]
                P_in += Pac[child_branch_idx]
                Q_in += Qac[child_branch_idx]
            end
        end
        
        # 更新支路功率
        Pac[branch_idx] = P_in
        Qac[branch_idx] = Q_in
    end
    
    # 直流网络正向扫描
    for n in Ndc:-1:2  # 从叶节点向根节点
        parent = parentsDC[n]
        branch_idx = 0
        
        # 查找连接到父节点的支路索引
        if haskey(dcbranch_idx, (parent, n))
            branch_idx = dcbranch_idx[(parent, n)]
        else
            error("找不到连接节点 $parent 和 $n 的支路")
        end
        
        # 计算流入节点n的功率
        P_in = dcP_load[n]
        
        # 如果节点有整流器，考虑整流器功率
        if n in C_rec
            rec_index = findfirst(x -> x == n, C_rec)
            P_in += P_rec[rec_index]  # 整流器输出为负，所以这里是加
        end
        
        # 加上流向所有子节点的功率
        if haskey(childrenDC, n) && !isempty(childrenDC[n])
            for child in childrenDC[n]
                child_branch_idx = dcbranch_idx[(n, child)]
                P_in += Pdc[child_branch_idx]
            end
        end
        
        # 更新支路功率
        Pdc[branch_idx] = P_in
    end
    
    # 反向扫描：计算节点电压
    # 交流网络反向扫描
    for n in 2:Nac
        parent = parentsAC[n]
        branch_idx = acbranch_idx[(parent, n)]
        
        r = branchAC_new[branch_idx, 3]  # 支路电阻
        x = branchAC_new[branch_idx, 4]  # 支路电抗
        
        # 计算电压降
        vac[n] = vac[parent] - 2 * (r * Pac[branch_idx] + x * Qac[branch_idx])
        
        # 检查电压是否在限制范围内
        if vac[n] < v_min
            println("警告：交流节点 $n 电压低于下限 $v_min")
            vac[n] = v_min
        elseif vac[n] > v_max
            println("警告：交流节点 $n 电压高于上限 $v_max")
            vac[n] = v_max
        end
    end
    
    # 直流网络反向扫描
    for n in 2:Ndc
        parent = parentsDC[n]
        branch_idx = dcbranch_idx[(parent, n)]
        
        r = branchDC_new[branch_idx, 3]  # 支路电阻
        
        # 计算电压降
        vdc[n] = vdc[parent] - 2 * r * Pdc[branch_idx]
        
        # 检查电压是否在限制范围内
        if vdc[n] < v_min
            println("警告：直流节点 $n 电压低于下限 $v_min")
            vdc[n] = v_min
        elseif vdc[n] > v_max
            println("警告：直流节点 $n 电压高于上限 $v_max")
            vdc[n] = v_max
        end
    end
    
    # 计算功率失配
    slack_pac = zeros(Nac-1)
    slack_qac = zeros(Nac-1)
    slack_pdc = zeros(Ndc-1)
    
    # 交流节点功率失配
    for n in 2:Nac
        # 计算节点n的净注入功率
        P_net = p_gen_fixed[n] - acP_load[n]
        Q_net = q_gen_fixed[n] - acQ_load[n]
        
        # 如果节点有逆变器，加上逆变器功率
        if n in C_inv
            inv_index = findfirst(x -> x == n, C_inv)
            P_net += P_inv[inv_index]
            Q_net += Q_inv[inv_index]
        end
        
        # 从父节点流入的功率
        parent = parentsAC[n]
        parent_branch_idx = acbranch_idx[(parent, n)]
        P_net += Pac[parent_branch_idx]
        Q_net += Qac[parent_branch_idx]
        
        # 流向子节点的功率
        if haskey(childrenAC, n) && !isempty(childrenAC[n])
            for child in childrenAC[n]
                child_branch_idx = acbranch_idx[(n, child)]
                P_net -= Pac[child_branch_idx]
                Q_net -= Qac[child_branch_idx]
            end
        end
        
        # 记录功率失配
        slack_pac[n-1] = P_net
        slack_qac[n-1] = Q_net
    end
    
    # 直流节点功率失配
    for n in 2:Ndc
        # 计算节点n的净注入功率
        P_net = -dcP_load[n]
        
        # 如果节点有整流器，减去整流器功率
        if n in C_rec
            rec_index = findfirst(x -> x == n, C_rec)
            P_net -= P_rec[rec_index]
        end
        
        # 从父节点流入的功率
        parent = parentsDC[n]
        parent_branch_idx = dcbranch_idx[(parent, n)]
        P_net += Pdc[parent_branch_idx]
        
        # 流向子节点的功率
        if haskey(childrenDC, n) && !isempty(childrenDC[n])
            for child in childrenDC[n]
                child_branch_idx = dcbranch_idx[(n, child)]
                P_net -= Pdc[child_branch_idx]
            end
        end
        
        # 记录功率失配
        slack_pdc[n-1] = P_net
    end
    
    # 计算最大功率失配
    max_slack = maximum([maximum(abs.(slack_pac)), maximum(abs.(slack_qac)), maximum(abs.(slack_pdc))])
    
    return Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack
end

function run_fixed_gen_test_case()
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
    v_min = 0.8
    v_max = 1.05
    P_inv = [0.1]./baseMVA
    Q_inv = [0.05]./baseMVA
    P_rec = [0.1111]./baseMVA
    η_rec = [0.9]
    η_inv = [0.9]
    p_gen_fixed = zeros(2)  # 初始化为正确的大小
    q_gen_fixed = zeros(2)  # 初始化为正确的大小
    
    network_data = process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    
    # 打印支路索引映射，以便调试
    println("\n交流支路索引映射:")
    for (branch, idx) in network_data["acbranch_idx"]
        println("支路 $branch 的索引: $idx")
    end
    
    println("\n直流支路索引映射:")
    for (branch, idx) in network_data["dcbranch_idx"]
        println("支路 $branch 的索引: $idx")
    end
    
    # 打印父节点和子节点映射，以便调试
    println("\n交流节点的父节点映射:")
    for (node, parent) in network_data["parentsAC"]
        println("节点 $node 的父节点: $parent")
    end
    
    println("\n交流节点的子节点映射:")
    for (node, children) in network_data["childrenAC"]
        println("节点 $node 的子节点: $children")
    end
    
    println("\n直流节点的父节点映射:")
    for (node, parent) in network_data["parentsDC"]
        println("节点 $node 的父节点: $parent")
    end
    
    println("\n直流节点的子节点映射:")
    for (node, children) in network_data["childrenDC"]
        println("节点 $node 的子节点: $children")
    end
    
    # Process the data
    Nac = network_data["Nac"]
    Ndc = network_data["Ndc"]
    busAC = network_data["busAC"]
    busDC = network_data["busDC"]
    acbranches = network_data["acbranches"]
    dcbranches = network_data["dcbranches"]
    acP_load = network_data["acP_load"]
    acQ_load = network_data["acQ_load"]
    dcP_load = network_data["dcP_load"]
    childrenAC = network_data["childrenAC"]
    childrenDC = network_data["childrenDC"]
    parentsAC = network_data["parentsAC"]
    parentsDC = network_data["parentsDC"]
    acbranch_idx = network_data["acbranch_idx"]
    dcbranch_idx = network_data["dcbranch_idx"]
    branchAC_new = network_data["branchAC"]
    branchDC_new = network_data["branchDC"]
    
    # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
    ac_node_map = network_data["ac_node_map"]
    dc_node_map = network_data["dc_node_map"]
    
    # 将原始的逆变器和整流器位置转换为重新编号后的位置
    C_inv_new = [ac_node_map[node] for node in C_inv]
    C_rec_new = [dc_node_map[node] for node in C_rec]
    
    println("原始逆变器位置: $C_inv, 重新编号后: $C_inv_new")
    println("原始整流器位置: $C_rec, 重新编号后: $C_rec_new")
    
    println("======= 混合交直流配电网络潮流计算 =======")
    println("\n网络参数:")
    println("交流节点数量: $Nac (包括根节点1)")
    println("直流节点数量: $Ndc (包括根节点1)")
    println("交流支路连接: $acbranches")
    println("直流支路连接: $dcbranches")
    println("逆变器位置: 节点 $C_inv_new")
    println("整流器位置: 节点 $C_rec_new")
    println("逆变器有功功率: $P_inv")
    println("逆变器无功功率: $Q_inv")
    println("整流器有功功率: $P_rec")
    println("交流负载有功功率: $acP_load")
    println("交流负载无功功率: $acQ_load")
    println("直流负载有功功率: $dcP_load")
    println("电压限制: [$v_min, $v_max]")
       
    # 求解潮流方程
    Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack = solve_distflow_equations(
        Nac, Ndc, C_inv_new, C_rec_new, P_inv, Q_inv, P_rec, acP_load, acQ_load, dcP_load, 
        v_min, v_max, p_gen_fixed, q_gen_fixed, childrenAC, childrenDC, parentsAC, parentsDC,
        acbranch_idx, dcbranch_idx, acbranches, dcbranches, branchAC_new, branchDC_new
    )
   
    # 输出结果
    println("\n求解完成!")
    println("最大功率失配: ", max_slack)
    
    # 创建新旧节点编号的反向映射
    ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
    dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
    
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
    
    # 输出功率失配（使用原始编号）
    println("\n功率失配:")
    println("交流有功功率失配:")
    for i in 1:Nac-1
        old_node = ac_new_to_old[i+1]  # i+1是因为失配变量是从第2个节点开始的
        println("节点 $old_node: $(slack_pac[i])")
    end
    
    println("交流无功功率失配:")
    for i in 1:Nac-1
        old_node = ac_new_to_old[i+1]
        println("节点 $old_node: $(slack_qac[i])")
    end
    
    println("直流有功功率失配:")
    for i in 1:Ndc-1
        old_node = dc_new_to_old[i+1]
        println("节点 $old_node: $(slack_pdc[i])")
    end
    
    return Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack
end

# 保留原有的网络处理函数
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

function create_branch_index_mapping(branches)
    # 创建支路索引映射：(from, to) => branch_index
    branch_idx = Dict{Tuple{Int, Int}, Int}()
    for (i, (from, to)) in enumerate(branches)
        branch_idx[(from, to)] = i
    end
    return branch_idx
end

function build_undirected_graph(branches)
    # 创建无向图（邻接表表示）
    graph = Dict{Int, Vector{Int}}()
    
    for (from, to) in branches
        if !haskey(graph, from)
            graph[from] = Int[]
        end
        if !haskey(graph, to)
            graph[to] = Int[]
        end
        push!(graph[from], to)
        push!(graph[to], from)  # 无向图，两个方向都添加
    end
    
    return graph
end

function create_tree_mappings(branches, n_nodes, root)
    # 首先构建无向图
    undirected_graph = build_undirected_graph(branches)
    
    # 创建子节点映射和父节点映射
    children = Dict(i => Int[] for i in 1:n_nodes)
    parents = Dict(i => 0 for i in 1:n_nodes)
    
    # 使用BFS从根节点开始构建树
    visited = falses(n_nodes)
    queue = [root]
    visited[root] = true
    
    while !isempty(queue)
        node = popfirst!(queue)
        
        for neighbor in get(undirected_graph, node, Int[])
            if !visited[neighbor]
                push!(children[node], neighbor)  # node是neighbor的父节点
                parents[neighbor] = node         # neighbor的父节点是node
                visited[neighbor] = true
                push!(queue, neighbor)
            end
        end
    end
    
    # 验证父节点映射和子节点映射的一致性
    for node in 1:n_nodes
        if node != root  # 根节点没有父节点
            parent = parents[node]
            if !(node in children[parent])
                error("映射不一致：节点 $node 的父节点是 $parent,但 $node 不在 $parent 的子节点列表中")
            end
        end
    end
    
    return children, parents
end

function process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
    # 识别松弛节点（根节点）
    ac_root, dc_root = identify_slack_nodes(busAC, busDC)
    println("原始交流侧松弛节点（根节点）: $ac_root")
    println("原始直流侧松弛节点（根节点）: $dc_root")
    
    # 对交流网络进行基于路径的重新编号，使根节点为1
    new_busAC, new_branchAC, ac_node_map = reindex_network_by_path(busAC, branchAC, ac_root)
    println("交流节点重新编号映射: $ac_node_map")
    
        # 对直流网络进行基于路径的重新编号，使根节点为1
        new_busDC, new_branchDC, dc_node_map = reindex_network_by_path(busDC, branchDC, dc_root)
        println("直流节点重新编号映射: $dc_node_map")
        
        # 对节点矩阵按节点编号排序
        ac_sort_idx = sortperm(new_busAC[:, 1])
        new_busAC = new_busAC[ac_sort_idx, :]
        
        dc_sort_idx = sortperm(new_busDC[:, 1])
        new_busDC = new_busDC[dc_sort_idx, :]
        
        # 更新根节点编号（两侧都应该是1）
        ac_root_new = ac_node_map[ac_root]  # 应该是1
        dc_root_new = dc_node_map[dc_root]  # 应该是1
        
        # 处理节点编号
        Nac = size(new_busAC, 1)
        Ndc = size(new_busDC, 1)
        
        # 转换支路数据为元组列表
        acbranches = convert_branch_matrix_to_tuples(new_branchAC)
        dcbranches = convert_branch_matrix_to_tuples(new_branchDC)
        
        # 基于重新编号后的网络构建树结构
        childrenAC, parentsAC = create_tree_mappings(acbranches, Nac, ac_root_new)
        childrenDC, parentsDC = create_tree_mappings(dcbranches, Ndc, dc_root_new)
        
        # 创建支路索引映射
        acbranch_idx = Dict{Tuple{Int, Int}, Int}()
        for (i, (from, to)) in enumerate(acbranches)
            acbranch_idx[(from, to)] = i
        end
        
        dcbranch_idx = Dict{Tuple{Int, Int}, Int}()
        for (i, (from, to)) in enumerate(dcbranches)
            dcbranch_idx[(from, to)] = i
        end
        
        # 提取负载数据并转换为p.u.
        # 现在节点矩阵已经按编号排序，所以行索引等于节点编号
        acP_load = new_busAC[:, 3] ./ baseMVA
        acQ_load = new_busAC[:, 4] ./ baseMVA
        dcP_load = new_busDC[:, 3] ./ baseMVA
        
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
            "childrenAC" => childrenAC,
            "parentsAC" => parentsAC,
            "childrenDC" => childrenDC,
            "parentsDC" => parentsDC,
            "ac_root" => ac_root_new,
            "dc_root" => dc_root_new,
            "busAC" => new_busAC,
            "busDC" => new_busDC,
            "branchAC" => new_branchAC,
            "branchDC" => new_branchDC,
            "ac_node_map" => ac_node_map,
            "dc_node_map" => dc_node_map
        )
    end
    
    # 添加迭代求解功能，实现牛顿-拉夫森法求解交直流混合潮流
    function solve_distflow_iteratively(
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
        v_min::Float64,            # 最小电压限制
        v_max::Float64,            # 最大电压限制
        p_gen_fixed::Vector{Float64}, # 节点固定有功发电量
        q_gen_fixed::Vector{Float64}, # 节点固定无功发电量
        childrenAC,
        childrenDC,
        parentsAC,
        parentsDC,
        acbranch_idx,
        dcbranch_idx,
        acbranches,
        dcbranches,
        branchAC_new,  # 交流支路参数矩阵
        branchDC_new   # 直流支路参数矩阵
    )
        # 初始化变量
        # 支路功率
        Pac = zeros(length(acbranches))
        Qac = zeros(length(acbranches))
        Pdc = zeros(length(dcbranches))
        
        # 节点电压
        vac = ones(Nac)  # 初始化所有交流节点电压为1.0 p.u.
        vdc = ones(Ndc)  # 初始化所有直流节点电压为1.0 p.u.
        
        # 迭代参数
        max_iter = 100
        tolerance = 1e-6
        converged = false
        iter = 0
        
        # 主迭代循环
        while !converged && iter < max_iter
            iter += 1
            
            # 保存上一次迭代的值，用于检查收敛性
            Pac_old = copy(Pac)
            Qac_old = copy(Qac)
            Pdc_old = copy(Pdc)
            vac_old = copy(vac)
            vdc_old = copy(vdc)
            
            # 正向扫描：计算支路功率
            # 交流网络正向扫描
            for n in Nac:-1:2  # 从叶节点向根节点
                parent = parentsAC[n]
                branch_idx = 0
                
                # 查找连接到父节点的支路索引
                if haskey(acbranch_idx, (parent, n))
                    branch_idx = acbranch_idx[(parent, n)]
                else
                    error("找不到连接节点 $parent 和 $n 的支路")
                end
                
                # 计算流入节点n的功率
                P_in = acP_load[n]
                Q_in = acQ_load[n]
                
                # 减去节点n的发电量
                if n in C_inv
                    inv_index = findfirst(x -> x == n, C_inv)
                    P_in -= P_inv[inv_index]
                    Q_in -= Q_inv[inv_index]
                end
                P_in -= p_gen_fixed[n]
                Q_in -= q_gen_fixed[n]
                
                # 加上流向所有子节点的功率
                if haskey(childrenAC, n) && !isempty(childrenAC[n])
                    for child in childrenAC[n]
                        child_branch_idx = acbranch_idx[(n, child)]
                        P_in += Pac[child_branch_idx]
                        Q_in += Qac[child_branch_idx]
                    end
                end
                
                # 更新支路功率
                Pac[branch_idx] = P_in
                Qac[branch_idx] = Q_in
            end
            
            # 直流网络正向扫描
            for n in Ndc:-1:2  # 从叶节点向根节点
                parent = parentsDC[n]
                branch_idx = 0
                
                # 查找连接到父节点的支路索引
                if haskey(dcbranch_idx, (parent, n))
                    branch_idx = dcbranch_idx[(parent, n)]
                else
                    error("找不到连接节点 $parent 和 $n 的支路")
                end
                
                # 计算流入节点n的功率
                P_in = dcP_load[n]
                
                # 如果节点有整流器，考虑整流器功率
                if n in C_rec
                    rec_index = findfirst(x -> x == n, C_rec)
                    P_in += P_rec[rec_index]  # 整流器输出为负，所以这里是加
                end
                
                # 加上流向所有子节点的功率
                if haskey(childrenDC, n) && !isempty(childrenDC[n])
                    for child in childrenDC[n]
                        child_branch_idx = dcbranch_idx[(n, child)]
                        P_in += Pdc[child_branch_idx]
                    end
                end
                
                # 更新支路功率
                Pdc[branch_idx] = P_in
            end
            
            # 反向扫描：计算节点电压
            # 交流网络反向扫描
            for n in 2:Nac
                parent = parentsAC[n]
                branch_idx = acbranch_idx[(parent, n)]
                
                r = branchAC_new[branch_idx, 3]  # 支路电阻
                x = branchAC_new[branch_idx, 4]  # 支路电抗
                
                # 计算电压降
                vac[n] = vac[parent] - 2 * (r * Pac[branch_idx] + x * Qac[branch_idx])
                
                # 检查电压是否在限制范围内
                if vac[n] < v_min
                    vac[n] = v_min
                elseif vac[n] > v_max
                    vac[n] = v_max
                end
            end
            
            # 直流网络反向扫描
            for n in 2:Ndc
                parent = parentsDC[n]
                branch_idx = dcbranch_idx[(parent, n)]
                
                r = branchDC_new[branch_idx, 3]  # 支路电阻
                
                # 计算电压降
                vdc[n] = vdc[parent] - 2 * r * Pdc[branch_idx]
                
                # 检查电压是否在限制范围内
                if vdc[n] < v_min
                    vdc[n] = v_min
                elseif vdc[n] > v_max
                    vdc[n] = v_max
                end
            end
            
            # 检查收敛性
            max_diff_P = maximum(abs.(Pac - Pac_old))
            max_diff_Q = maximum(abs.(Qac - Qac_old))
            max_diff_Pdc = maximum(abs.(Pdc - Pdc_old))
            max_diff_vac = maximum(abs.(vac - vac_old))
            max_diff_vdc = maximum(abs.(vdc - vdc_old))
            
            max_diff = maximum([max_diff_P, max_diff_Q, max_diff_Pdc, max_diff_vac, max_diff_vdc])
            
            if max_diff < tolerance
                converged = true
                println("迭代收敛于第 $iter 次迭代，最大差值: $max_diff")
            elseif iter % 10 == 0
                println("迭代 $iter: 最大差值 = $max_diff")
            end
        end
        
        if !converged
            println("警告：达到最大迭代次数 $max_iter,未收敛。最大差值: $max_diff")
        end
        
        # 计算功率失配
        slack_pac = zeros(Nac-1)
        slack_qac = zeros(Nac-1)
        slack_pdc = zeros(Ndc-1)
        
        # 交流节点功率失配
        for n in 2:Nac
            # 计算节点n的净注入功率
            P_net = p_gen_fixed[n] - acP_load[n]
            Q_net = q_gen_fixed[n] - acQ_load[n]
            
            # 如果节点有逆变器，加上逆变器功率
            if n in C_inv
                inv_index = findfirst(x -> x == n, C_inv)
                P_net += P_inv[inv_index]
                Q_net += Q_inv[inv_index]
            end
            
            # 从父节点流入的功率
            parent = parentsAC[n]
            parent_branch_idx = acbranch_idx[(parent, n)]
            P_net += Pac[parent_branch_idx]
            Q_net += Qac[parent_branch_idx]
            
            # 流向子节点的功率
            if haskey(childrenAC, n) && !isempty(childrenAC[n])
                for child in childrenAC[n]
                    child_branch_idx = acbranch_idx[(n, child)]
                    P_net -= Pac[child_branch_idx]
                    Q_net -= Qac[child_branch_idx]
                end
            end
            
            # 记录功率失配
            slack_pac[n-1] = P_net
            slack_qac[n-1] = Q_net
        end
        
        # 直流节点功率失配
        for n in 2:Ndc
            # 计算节点n的净注入功率
            P_net = -dcP_load[n]
            
            # 如果节点有整流器，减去整流器功率
            if n in C_rec
                rec_index = findfirst(x -> x == n, C_rec)
                P_net -= P_rec[rec_index]
            end
            
            # 从父节点流入的功率
            parent = parentsDC[n]
            parent_branch_idx = dcbranch_idx[(parent, n)]
            P_net += Pdc[parent_branch_idx]
            
            # 流向子节点的功率
            if haskey(childrenDC, n) && !isempty(childrenDC[n])
                for child in childrenDC[n]
                    child_branch_idx = dcbranch_idx[(n, child)]
                    P_net -= Pdc[child_branch_idx]
                end
            end
            
            # 记录功率失配
            slack_pdc[n-1] = P_net
        end
        
        # 计算最大功率失配
        max_slack = maximum([maximum(abs.(slack_pac)), maximum(abs.(slack_qac)), maximum(abs.(slack_pdc))])
        
        return Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack, iter, converged
    end
    
    # 修改主函数以使用迭代求解
    function run_fixed_gen_test_case_iterative()
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
        v_min = 0.8
        v_max = 1.05
        P_inv = [0.1]./baseMVA
        Q_inv = [0.05]./baseMVA
        P_rec = [0.1111]./baseMVA
        η_rec = [0.9]
        η_inv = [0.9]
        p_gen_fixed = zeros(2)  # 初始化为正确的大小
        q_gen_fixed = zeros(2)  # 初始化为正确的大小
        
        network_data = process_network_data(busAC, busDC, branchAC, branchDC, baseMVA)
        
        # 打印支路索引映射，以便调试
        println("\n交流支路索引映射:")
        for (branch, idx) in network_data["acbranch_idx"]
            println("支路 $branch 的索引: $idx")
        end
        
        println("\n直流支路索引映射:")
        for (branch, idx) in network_data["dcbranch_idx"]
            println("支路 $branch 的索引: $idx")
        end
        
        # 打印父节点和子节点映射，以便调试
        println("\n交流节点的父节点映射:")
        for (node, parent) in network_data["parentsAC"]
            println("节点 $node 的父节点: $parent")
        end
        
        println("\n交流节点的子节点映射:")
        for (node, children) in network_data["childrenAC"]
            println("节点 $node 的子节点: $children")
        end
        
        println("\n直流节点的父节点映射:")
        for (node, parent) in network_data["parentsDC"]
            println("节点 $node 的父节点: $parent")
        end
        
        println("\n直流节点的子节点映射:")
        for (node, children) in network_data["childrenDC"]
            println("节点 $node 的子节点: $children")
        end
        
        # Process the data
        Nac = network_data["Nac"]
        Ndc = network_data["Ndc"]
        busAC = network_data["busAC"]
        busDC = network_data["busDC"]
        acbranches = network_data["acbranches"]
        dcbranches = network_data["dcbranches"]
        acP_load = network_data["acP_load"]
        acQ_load = network_data["acQ_load"]
        dcP_load = network_data["dcP_load"]
        childrenAC = network_data["childrenAC"]
        childrenDC = network_data["childrenDC"]
        parentsAC = network_data["parentsAC"]
        parentsDC = network_data["parentsDC"]
        acbranch_idx = network_data["acbranch_idx"]
        dcbranch_idx = network_data["dcbranch_idx"]
        branchAC_new = network_data["branchAC"]
        branchDC_new = network_data["branchDC"]
        
        # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
        ac_node_map = network_data["ac_node_map"]
        dc_node_map = network_data["dc_node_map"]
        
        # 将原始的逆变器和整流器位置转换为重新编号后的位置
        C_inv_new = [ac_node_map[node] for node in C_inv]
        C_rec_new = [dc_node_map[node] for node in C_rec]
        
        println("原始逆变器位置: $C_inv, 重新编号后: $C_inv_new")
        println("原始整流器位置: $C_rec, 重新编号后: $C_rec_new")
        
        println("======= 混合交直流配电网络潮流计算 (迭代求解) =======")
        println("\n网络参数:")
        println("交流节点数量: $Nac (包括根节点1)")
        println("直流节点数量: $Ndc (包括根节点1)")
        println("交流支路连接: $acbranches")
        println("直流支路连接: $dcbranches")
        println("逆变器位置: 节点 $C_inv_new")
        println("整流器位置: 节点 $C_rec_new")
        println("逆变器有功功率: $P_inv")
        println("逆变器无功功率: $Q_inv")
        println("整流器有功功率: $P_rec")
        println("交流负载有功功率: $acP_load")
        println("交流负载无功功率: $acQ_load")
        println("直流负载有功功率: $dcP_load")
        println("电压限制: [$v_min, $v_max]")
           
        # 求解潮流方程（迭代方法）
        Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack, iter, converged = solve_distflow_iteratively(
            Nac, Ndc, C_inv_new, C_rec_new, P_inv, Q_inv, P_rec, acP_load, acQ_load, dcP_load, 
            v_min, v_max, p_gen_fixed, q_gen_fixed, childrenAC, childrenDC, parentsAC, parentsDC,
            acbranch_idx, dcbranch_idx, acbranches, dcbranches, branchAC_new, branchDC_new
        )
       
        # 输出结果
        println("\n求解完成!")
        println("迭代次数: $iter")
        println("收敛状态: $(converged ? "已收敛" : "未收敛")")
        println("最大功率失配: ", max_slack)
        
        # 创建新旧节点编号的反向映射
        ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
        dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
        
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
        
        # 输出功率失配（使用原始编号）
        println("\n功率失配:")
        println("交流有功功率失配:")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]  # i+1是因为失配变量是从第2个节点开始的
            println("节点 $old_node: $(slack_pac[i])")
        end
        
        println("交流无功功率失配:")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]
            println("节点 $old_node: $(slack_qac[i])")
        end
        
        println("直流有功功率失配:")
        for i in 1:Ndc-1
            old_node = dc_new_to_old[i+1]
            println("节点 $old_node: $(slack_pdc[i])")
        end
        
        return Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc, max_slack, iter, converged
    end
    
    # 运行主程序
    run_fixed_gen_test_case_iterative()
    
