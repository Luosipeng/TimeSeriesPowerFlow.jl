"""
This file is used to test the linear AC/DC hybrid distribution flow model.
"""

using JuMP
using Ipopt
using LinearAlgebra

# 定义配电网络参数
function create_distflow_model(
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
    # 创建优化模型
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
     
    # 决策变量
    @variable(model, Pac[1:length(acbranches)]) # 支路有功功率
    @variable(model, Qac[1:length(acbranches)]) # 支路无功功率
    @variable(model, Pdc[1:length(dcbranches)]) # 直流支路有功功率
    
    # 节点电压变量 - 注意根节点电压固定为1.0
    @variable(model, vac[1:Nac])  # 交流节点电压幅值的平方
    @variable(model, vdc[1:Ndc])  # 直流节点电压幅值的平方
    
    # 为非根节点设置电压约束
    for i in 2:Nac
        @constraint(model, v_min <= vac[i] <= v_max)
    end
    
    for i in 2:Ndc
        @constraint(model, v_min <= vdc[i] <= v_max)
    end
    
    # 固定根节点电压为1.0
    @constraint(model, vac[1] == 1.0)
    @constraint(model, vdc[1] == 1.0)
    
    @variable(model, slack_pac[1:Nac-1])  # 有功功率失配松弛变量，无约束
    @variable(model, slack_qac[1:Nac-1])  # 无功功率失配松弛变量，无约束
    @variable(model, slack_pdc[1:Ndc-1])  # 有功功率失配松弛变量，无约束

    # 无穷范数目标函数变量
    @variable(model, max_slack >= 0)  # 无穷范数目标函数变量

    # 目标函数：最小化最大失配量的绝对值（无穷范数）
    @objective(model, Min, max_slack)

    # 无穷范数约束（考虑正负值）
    for i in 1:Nac-1
        @constraint(model, slack_pac[i] <= max_slack)
        @constraint(model, -slack_pac[i] <= max_slack)
        @constraint(model, slack_qac[i] <= max_slack)
        @constraint(model, -slack_qac[i] <= max_slack)
    end
    
    for i in 1:Ndc-1
        @constraint(model, slack_pdc[i] <= max_slack)
        @constraint(model, -slack_pdc[i] <= max_slack)
    end
    
    # 逆变器功率
    P_inverter = zeros(Nac-1,1)
    P_inverter[C_inv.-1] = P_inv

    Q_inverter = zeros(Nac-1,1)
    Q_inverter[C_inv.-1] = Q_inv

    # 整流器功率
    P_rectifier = zeros(Ndc-1,1)
    P_rectifier[C_rec.-1] = P_rec
    
    # 交流节点功率平衡方程
    acP_load = acP_load[2:end]  # 去掉根节点
    acP_n = Pac[map(k -> acbranch_idx[k], 2:Nac)]
   
    acQ_load = acQ_load[2:end]  # 去掉根节点
    acQ_n = Qac[map(k -> acbranch_idx[k], 2:Nac)]

    @constraint(model,Aac' * acP_n +acP_load - P_inverter .== slack_pac)
    @constraint(model,Aac' * acQ_n +acQ_load - Q_inverter .== slack_qac)

    # 直流节点功率平衡方程
    dcP_load = dcP_load[2:end]  # 去掉根节点
    dcP_n = Pdc[map(k -> dcbranch_idx[k], 2:Ndc)]

    @constraint(model,Adc' * dcP_n +dcP_load + P_rectifier .== slack_pdc)

    # 交流节点电压平衡方程
    @constraint(model, -2.0 * Rac * acP_n - 2.0 * Xac * acQ_n .+vac[1:end-1] .== vac[2:end])

    # 直流节点电压平衡方程
    @constraint(model, -2.0 * Rdc * dcP_n .+vdc[1:end-1] .== vdc[2:end])

    return model, Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc
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

function run_fixed_gen_test_case(mpc)
    baseMVA = mpc["baseMVA"]  # 基准功率
    # 网络参数
    busAC = mpc["busAC"]
    busDC = mpc["busDC"]
    branchAC = mpc["branchAC"]
    branchDC = mpc["branchDC"]
    C_inv = mpc["Ci"]
    C_rec = mpc["Cr"]
    v_min = 0.8
    v_max = 1.05
    P_inv = mpc["P_inv"] ./baseMVA
    Q_inv = mpc["Q_inv"] ./baseMVA
    P_rec = mpc["P_inv_dc"] ./baseMVA
    
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
    
    # 更新 C_inv 和 C_rec 以匹配重新编号后的节点
    C_inv_new = [ac_node_map[node] for node in C_inv]
    C_rec_new = [dc_node_map[node] for node in C_rec]
    
    # 创建交流支路节点关联矩阵
    Aac = build_incidence_matrix(Nac, acbranches)

    # 创建直流支路节点关联矩阵
    Adc = build_incidence_matrix(Ndc, dcbranches)

    model, Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc = create_distflow_model(
        Nac, Ndc, C_inv_new, C_rec_new, P_inv, Q_inv, P_rec, acP_load, acQ_load, 
        dcP_load, v_min, v_max, acbranch_idx, dcbranch_idx, acbranches, dcbranches,
        Aac, Adc, Rac, Xac, Rdc
    )
   
    # 求解模型
    optimize!(model)
    
    # 输出结果
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
        println("\n优化成功!")
        println("目标函数值 (最大功率失配): ", objective_value(model))
        
        # 创建新旧节点编号的反向映射
        ac_new_to_old = Dict(v => k for (k, v) in ac_node_map)
        dc_new_to_old = Dict(v => k for (k, v) in dc_node_map)
        
        # 输出交流支路功率（使用原始编号）
        println("\n交流支路功率:")
        for i in eachindex(acbranches)
            from_new, to_new = acbranches[i]
            from_old = ac_new_to_old[from_new]
            to_old = ac_new_to_old[to_new]
            println("支路 ($from_old, $to_old): P = $(value(Pac[i])), Q = $(value(Qac[i]))")
        end
        
        # 输出直流支路功率（使用原始编号）
        println("\n直流支路功率:")
        for i in eachindex(dcbranches)
            from_new, to_new = dcbranches[i]
            from_old = dc_new_to_old[from_new]
            to_old = dc_new_to_old[to_new]
            println("支路 ($from_old, $to_old): P = $(value(Pdc[i]))")
        end
        
        # 输出交流节点电压（使用原始编号）
        println("\n交流节点电压:")
        for i in 1:Nac
            old_node = ac_new_to_old[i]
            println("节点 $old_node: v = $(sqrt(value(vac[i])))")
        end
        
        # 输出直流节点电压（使用原始编号）
        println("\n直流节点电压:")
        for i in 1:Ndc
            old_node = dc_new_to_old[i]
            println("节点 $old_node: v = $(sqrt(value(vdc[i])))")
        end
        
        # 输出功率失配（使用原始编号）
        println("\n功率失配:")
        println("交流有功功率失配:")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]  # i+1是因为失配变量是从第2个节点开始的
            println("节点 $old_node: $(value(slack_pac[i]))")
        end
        
        println("交流无功功率失配:")
        for i in 1:Nac-1
            old_node = ac_new_to_old[i+1]
            println("节点 $old_node: $(value(slack_qac[i]))")
        end
        
        println("直流有功功率失配:")
        for i in 1:Ndc-1
            old_node = dc_new_to_old[i+1]
            println("节点 $old_node: $(value(slack_pdc[i]))")
        end
    else
        println("优化失败: ", termination_status(model))
    end
    
    return model, Pac, Qac, Pdc, vac, vdc, slack_pac, slack_qac, slack_pdc
end

# 运行测试用例
run_fixed_gen_test_case(mpc)
