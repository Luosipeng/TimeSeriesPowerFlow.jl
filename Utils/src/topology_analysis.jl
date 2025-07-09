"""
    create_node_mapping(case::JuliaPowerCase)

为节点创建编号映射，检测重复节点。
返回节点名称到ID的映射字典。
"""
function create_node_mapping(case::JuliaPowerCase)
    # 从case中提取所有节点名称
    nodes = [bus.name for bus in case.busesAC]
    
    count_dict = Dict{String, Int}()
    duplicates = String[]
    
    # 统计节点出现次数，检测重复节点
    for node in nodes
        node_str = string(node)
        count_dict[node_str] = get(count_dict, node_str, 0) + 1
        if count_dict[node_str] > 1 && !(node_str in duplicates)
            push!(duplicates, node_str)
        end
    end
    
    # 输出重复节点警告
    if !isempty(duplicates)
        println("\n警告：发现重复的节点名称：")
        for name in duplicates
            println(" - ", name, " (出现 ", count_dict[name], " 次)")
        end
        println()
    end
    
    # 创建节点到ID的映射
    node_dict = Dict{String, Int}()
    id = 1
    for node in unique(nodes)
        node_str = string(node)
        if !haskey(node_dict, node_str)
            node_dict[node_str] = id
            id += 1
        end
    end
    return node_dict
end

"""
    filter_active_edges(case::JuliaPowerCase)

过滤出状态值为true的活动连接。
返回活动连接的列表。
"""
function filter_active_edges(case::JuliaPowerCase)
    # 从case中提取所有边
    active_edges = []
    
    # 处理交流线路
    for line in case.branchesAC
        if line.in_service
            push!(active_edges, (line.from_bus, line.to_bus, line))
        end
    end
    
    # 处理变压器
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            push!(active_edges, (transformer.hv_bus, transformer.lv_bus, transformer))
        end
    end
    
    # 处理断路器
    for hvcb in case.hvcbs
        if hvcb.closed
            push!(active_edges, (hvcb.bus_from, hvcb.bus_to, hvcb))
        end
    end
    
    return active_edges
end

"""
    update_edge_ids(edges, node_dict)

更新连接表的始末端ID。
返回边的ID对列表。
"""
function update_edge_ids(edges, node_dict)
    edge_ids = []
    
    for (from_bus, to_bus, _) in edges
        from_id = get(node_dict, string(from_bus), 0)
        to_id = get(node_dict, string(to_bus), 0)
        push!(edge_ids, (from_id, to_id))
    end
    
    return edge_ids
end

"""
    is_switch_element(element)

判断元件是否为开关类型。
"""
function is_switch_element(element)
    # 检查元素是否为断路器类型
    if typeof(element) == DistributionSystem.HighVoltageCircuitBreaker
        return true
    end
    
    # 其他可能的开关类型检查
    if hasfield(typeof(element), :type)
        type_str = lowercase(strip(string(element.type)))
        switch_types = ["开关", "断路器", "隔离开关", "负荷开关", "刀闸", "刀开关"]
        is_switch = any(type -> type == type_str || occursin(type, type_str), switch_types)
        is_not_room = !occursin("房", type_str) && !occursin("室", type_str) && !occursin("站", type_str)
        return is_switch && is_not_room
    end
    
    return false
end

"""
    is_substation(element)

判断元件是否为变电站类型。
"""
function is_substation(element)
    if hasfield(typeof(element), :type)
        type_str = lowercase(strip(string(element.type)))
        return occursin("变电站", type_str) || occursin("变电所", type_str)
    end
    return false
end

"""
    count_switch_ports(case::JuliaPowerCase)

统计开关端口数。
返回节点ID到端口数的映射字典。
"""
function count_switch_ports(case::JuliaPowerCase)
    port_counts = Dict{Int, Int}()
    
    # 初始化所有节点的连接数为0
    for bus in case.busesAC
        port_counts[bus.bus_id] = 0
    end
    
    # 统计每个节点的连接数
    # 处理交流线路
    for line in case.branchesAC
        if line.in_service
            port_counts[line.from_bus] = get(port_counts, line.from_bus, 0) + 1
            port_counts[line.to_bus] = get(port_counts, line.to_bus, 0) + 1
        end
    end
    
    # 处理变压器
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            port_counts[transformer.hv_bus] = get(port_counts, transformer.hv_bus, 0) + 1
            port_counts[transformer.lv_bus] = get(port_counts, transformer.lv_bus, 0) + 1
        end
    end
    
    # 处理断路器
    for hvcb in case.hvcbs
        if hvcb.closed
            port_counts[hvcb.bus_from] = get(port_counts, hvcb.bus_from, 0) + 1
            port_counts[hvcb.bus_to] = get(port_counts, hvcb.bus_to, 0) + 1
        end
    end
    
    return port_counts
end

"""
    create_virtual_node(case::JuliaPowerCase, bus_id::Int, virtual_name::String)

创建虚拟节点数据。
"""
function create_virtual_node(case::JuliaPowerCase, bus_id::Int, virtual_name::String)
    # 找到原始节点
    original_bus = nothing
    for bus in case.busesAC
        if bus.bus_id == bus_id
            original_bus = bus
            break
        end
    end
    
    if original_bus === nothing
        error("找不到ID为 $bus_id 的节点")
    end
    
    # 创建虚拟节点（复制原始节点的属性）
    virtual_bus = deepcopy(original_bus)
    virtual_bus.index = length(case.busesAC) + 1  # 新ID
    virtual_bus.bus_id = length(case.busesAC) + 1  # 新ID
    virtual_bus.name = virtual_name
    
    return virtual_bus
end

"""
    create_virtual_connection(case::JuliaPowerCase, from_bus::Int, to_bus::Int)

创建虚拟连接数据。
"""
function create_virtual_connection(case::JuliaPowerCase, from_bus::Int, to_bus::Int)
    # 创建一个新的断路器作为虚拟连接
    virtual_cb = DistributionSystem.HighVoltageCircuitBreaker(
        length(case.hvcbs) + 1,  # 新ID
        "VIRTUAL_CB_$(from_bus)_$(to_bus)",  # 名称
        from_bus,  # 起始节点
        to_bus,    # 终止节点
        "l",       # 类型
        0,         # 序号
        true,     # 闭合状态
        "CB",      # 类型
        0.0,       # 额定电流
        true       # 工作状态
    )
    
    return virtual_cb
end

"""
    process_switch_nodes(case::JuliaPowerCase)

处理开关节点，添加虚拟节点。
返回更新后的case对象和虚拟连接标记。
"""
function process_switch_nodes(case::JuliaPowerCase)
    # 深拷贝case以便修改
    new_case = deepcopy(case)
    
    # 创建一个集合来标记虚拟连接
    virtual_connections = Set{Tuple{Int, Int}}()
    
    # 获取端口数
    port_counts = count_switch_ports(case)
    
    # 创建一个集合来跟踪已处理的断路器
    processed_hvcbs = Set{Int}()
    
    # 找出所有开关元件
    switch_elements = []
    
    # 检查断路器
    for hvcb in case.hvcbs
        if hvcb.closed && !(hvcb.index in processed_hvcbs)
            bus_from = hvcb.bus_from
            bus_to = hvcb.bus_to
            
            # 检查两端节点是否都满足条件
            from_eligible = get(port_counts, bus_from, 0) > 1
            to_eligible = get(port_counts, bus_to, 0) > 1
            
            if from_eligible && to_eligible
                # 如果两端都满足条件，选择端口数更多的一端
                if get(port_counts, bus_from, 0) >= get(port_counts, bus_to, 0)
                    push!(switch_elements, (bus_from, hvcb))
                else
                    push!(switch_elements, (bus_to, hvcb))
                end
            elseif from_eligible
                push!(switch_elements, (bus_from, hvcb))
            elseif to_eligible
                push!(switch_elements, (bus_to, hvcb))
            end
            
            # 标记此断路器为已处理
            push!(processed_hvcbs, hvcb.index)
        end
    end
    
    # 处理每个开关节点
    for (bus_id, element) in switch_elements
        port_count = get(port_counts, bus_id, 0)
        
        # 处理任意端口数的情况（端口数至少为2）
        if port_count >= 2
            # 获取原始节点名称
            original_name = ""
            for bus in case.busesAC
                if bus.bus_id == bus_id
                    original_name = bus.name
                    break
                end
            end
            
            # 创建一个虚拟节点
            virtual_name = string(original_name) * "_虚拟节点"
            virtual_bus = create_virtual_node(new_case, bus_id, virtual_name)
            push!(new_case.busesAC, virtual_bus)
            
            # 连接原始节点和虚拟节点
            virtual_cb = create_virtual_connection(new_case, bus_id, virtual_bus.bus_id)
            push!(new_case.hvcbs, virtual_cb)
            
            # 添加第一个虚拟连接（原始节点-虚拟节点）
            if bus_id < virtual_bus.bus_id
                push!(virtual_connections, (bus_id, virtual_bus.bus_id))
            else
                push!(virtual_connections, (virtual_bus.bus_id, bus_id))
            end
            
            # 获取开关元件的另一端节点
            other_bus_id = element.bus_from == bus_id ? element.bus_to : element.bus_from
            
            # 创建新断路器连接虚拟节点和另一端节点
            new_hvcb = deepcopy(element)
            new_hvcb.index = maximum([cb.index for cb in new_case.hvcbs]) + 1
            
            if element.bus_from == bus_id
                new_hvcb.bus_from = virtual_bus.bus_id
                new_hvcb.bus_to = other_bus_id
            else
                new_hvcb.bus_from = other_bus_id
                new_hvcb.bus_to = virtual_bus.bus_id
            end
            
            push!(new_case.hvcbs, new_hvcb)
            
            # 添加第二个虚拟连接（虚拟节点-另一端节点）
            if virtual_bus.bus_id < other_bus_id
                push!(virtual_connections, (virtual_bus.bus_id, other_bus_id))
            else
                push!(virtual_connections, (other_bus_id, virtual_bus.bus_id))
            end
            
            # 将原始断路器设为不闭合
            for j in 1:length(new_case.hvcbs)
                if new_case.hvcbs[j].index == element.index
                    new_case.hvcbs[j].closed = false
                    break
                end
            end
        end
    end
    
    # 更新节点名称到ID的映射
    new_case.bus_name_to_id = Dict{String, Int}()
    for bus in new_case.busesAC
        new_case.bus_name_to_id[bus.name] = bus.bus_id
    end
    
    return new_case, virtual_connections
end



"""
    identify_partitions(case::JuliaPowerCase)

识别网络中的分区。
返回节点分区和边分区。
"""
function identify_partitions(case::JuliaPowerCase)
    # 创建图结构
    G = SimpleGraph(length(case.busesAC))
    
    # 添加边
    for line in case.branchesAC
        if line.in_service
            add_edge!(G, line.from_bus, line.to_bus)
        end
    end
    
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            add_edge!(G, transformer.hv_bus, transformer.lv_bus)
        end
    end
    
    for hvcb in case.hvcbs
        if hvcb.closed && hvcb.in_service
            add_edge!(G, hvcb.bus_from, hvcb.bus_to)
        end
    end
    
    # 找出连通分量
    components = connected_components(G)
    
    # 创建节点分区和边分区
    node_partitions = Dict{Int, Int}()
    edge_partitions = Dict{Tuple{Int, Int}, Int}()
    
    # 为每个节点分配分区
    for (partition_id, component) in enumerate(components)
        for node_id in component
            node_partitions[node_id] = partition_id
        end
    end
    
    # 为每条边分配分区
    for line in case.branchesAC
        if line.in_service
            from_partition = get(node_partitions, line.from_bus, 0)
            to_partition = get(node_partitions, line.to_bus, 0)
            
            # 边应该属于同一个分区
            if from_partition == to_partition
                edge_partitions[(line.from_bus, line.to_bus)] = from_partition
            end
        end
    end
    
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            from_partition = get(node_partitions, transformer.hv_bus, 0)
            to_partition = get(node_partitions, transformer.lv_bus, 0)
            
            if from_partition == to_partition
                edge_partitions[(transformer.hv_bus, transformer.lv_bus)] = from_partition
            end
        end
    end
    
    for hvcb in case.hvcbs
        if hvcb.closed
            from_partition = get(node_partitions, hvcb.bus_from, 0)
            to_partition = get(node_partitions, hvcb.bus_to, 0)
            
            if from_partition == to_partition
                edge_partitions[(hvcb.bus_from, hvcb.bus_to)] = from_partition
            end
        end
    end
    
    return node_partitions, edge_partitions
end

"""
    get_edge_endpoints(e)

获取边的源节点和目标节点。
"""
function get_edge_endpoints(e)
    # 尝试不同的方法获取边的端点
    try
        return (e.src, e.dst)
    catch
        try
            # 如果边是元组，直接返回
            return e
        catch
            # 如果以上方法都失败，尝试将边转换为字符串并解析
            edge_str = string(e)
            # 假设格式为 "Edge 1 => 2" 或类似格式
            m = match(r"Edge\s+(\d+)\s*=>\s*(\d+)", edge_str)
            if m !== nothing
                return (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
            end
            # 如果所有方法都失败，抛出错误
            error("无法获取边 $e 的端点")
        end
    end
end

"""
    dfs_tree(graph, root)

使用深度优先搜索构建生成树。
"""
function dfs_tree(graph, root)
    tree_edges = Vector{Tuple{Int, Int}}()
    visited = falses(nv(graph))
    parent = zeros(Int, nv(graph))
    
    # 辅助函数，执行DFS
    function dfs_helper(u)
        visited[u] = true
        for v in neighbors(graph, u)
            if !visited[v]
                push!(tree_edges, (u, v))
                parent[v] = u
                dfs_helper(v)
            end
        end
    end
    
    # 对每个连通分量执行DFS
    for v in 1:nv(graph)
        if !visited[v]
            dfs_helper(v)
        end
    end
    
    return tree_edges, parent
end

"""
    find_path_to_lca(u, v, parent)

找到从节点u和v到它们的最近公共祖先的路径。
"""
function find_path_to_lca(u, v, parent)
    # 从u到根的路径
    path_u_to_root = Vector{Int}()
    current = u
    while current != 0
        push!(path_u_to_root, current)
        current = parent[current]
    end
    
    # 从v向上到最近公共祖先
    path_v_to_lca = Vector{Int}()
    current = v
    lca_found = false
    lca = 0
    
    while current != 0
        if current in path_u_to_root
            # 找到最近公共祖先
            lca = current
            lca_found = true
            break
        end
        
        push!(path_v_to_lca, current)
        current = parent[current]
    end
    
    if !lca_found
        # 如果没有找到公共祖先，返回空路径
        return Vector{Int}(), false
    end
    
    # 截断path_u_to_root到lca
    path_u_to_lca = Vector{Int}()
    for node in path_u_to_root
        if node == lca
            break
        end
        push!(path_u_to_lca, node)
    end
    
    # 构建从u到v的路径
    path_u_to_v = copy(path_u_to_lca)
    
    # 将path_v_to_lca反转并添加到path_u_to_v
    for i in length(path_v_to_lca):-1:1
        push!(path_u_to_v, path_v_to_lca[i])
    end
    
    return path_u_to_v, true
end

"""
    find_fundamental_cycles(G::SimpleGraph)

使用基于深度优先搜索的方法找出无向图中的基本环路。
"""
function find_fundamental_cycles(G::SimpleGraph)
    n = nv(G)  # 获取图中的节点数
    
    # 创建一个映射，将图中的实际节点ID映射到连续的索引
    node_to_index = Dict{Int, Int}()
    index_to_node = Dict{Int, Int}()
    
    # 填充映射
    index = 1
    for v in vertices(G)
        node_to_index[v] = index
        index_to_node[index] = v
        index += 1
    end
    
    # 创建一个新的图，使用连续的索引
    G_continuous = SimpleGraph(n)
    
    # 添加原图中的边，但使用连续的索引
    for e in edges(G)
        src_node, dst_node = get_edge_endpoints(e)
        add_edge!(G_continuous, node_to_index[src_node], node_to_index[dst_node])
    end
    
    # 现在在连续索引的图上查找环路
    cycles = Vector{Vector{Int}}()
    
    # 对于每个连通分量，找出基本环路
    components = connected_components(G_continuous)
    
    for component in components
        if length(component) > 0
            root = component[1]
            tree_edges, parent = dfs_tree(G_continuous, root)
            
            # 收集非树边
            non_tree_edges = Vector{Tuple{Int, Int}}()
            for e in edges(G_continuous)
                src_idx, dst_idx = get_edge_endpoints(e)
                if !((src_idx, dst_idx) in tree_edges || (dst_idx, src_idx) in tree_edges)
                    push!(non_tree_edges, (src_idx, dst_idx))
                end
            end
            
            # 对于每个非树边，找到一个环路
            for (u, v) in non_tree_edges
                path_u_to_v, success = find_path_to_lca(u, v, parent)
                
                if !success
                    # 如果没有找到公共祖先，可能是图不连通
                    continue
                end
                
                # 添加边(v,u)完成环路
                push!(path_u_to_v, u)
                
                # 将连续索引映射回原始节点ID
                original_cycle = [index_to_node[idx] for idx in path_u_to_v]
                push!(cycles, original_cycle)
            end
        end
    end
    
    return cycles
end

"""
    create_and_plot_graph_by_partition(case::JuliaPowerCase, node_partitions, edge_partitions, virtual_connections)

为每个分区创建图并找出环路。
"""
function create_and_plot_graph_by_partition(case::JuliaPowerCase, node_partitions, edge_partitions, virtual_connections)
    # 获取唯一的分区ID
    partition_ids = unique(values(node_partitions))
    
    # 为每个分区创建和绘制图形
    cycles_by_partition = Dict{Int, Vector{Vector{Int}}}()
    
    for partition_id in partition_ids
        # 创建该分区的子图
        partition_nodes = [node_id for (node_id, part_id) in node_partitions if part_id == partition_id]
        partition_edges = [(from, to) for ((from, to), part_id) in edge_partitions if part_id == partition_id]
        
        if isempty(partition_nodes)
            println("分区 $partition_id 没有节点，跳过")
            continue
        end
        
        # 创建图，使用实际节点ID
        G = SimpleGraph()
        
        # 添加节点
        for node_id in partition_nodes
            add_vertex!(G)
        end
        
        # 创建节点ID到图索引的映射
        node_to_vertex = Dict{Int, Int}()
        vertex_to_node = Dict{Int, Int}()
        
        for (i, node_id) in enumerate(partition_nodes)
            node_to_vertex[node_id] = i
            vertex_to_node[i] = node_id
        end
        
        # 添加边，使用图索引
        for (from, to) in partition_edges
            if haskey(node_to_vertex, from) && haskey(node_to_vertex, to)
                from_vertex = node_to_vertex[from]
                to_vertex = node_to_vertex[to]
                add_edge!(G, from_vertex, to_vertex)
            else
                println("警告：边 ($from, $to) 的一个或两个端点不在分区 $partition_id 中")
            end
        end
        
        # 找出环路
        cycles = Vector{Vector{Int}}()
        
        try
            # 使用自定义的环路检测函数
            graph_cycles = find_fundamental_cycles(G)
            
            # 将图索引映射回原始节点ID
            for cycle in graph_cycles
                original_cycle = [vertex_to_node[v] for v in cycle]
                push!(cycles, original_cycle)
            end
            
            println("分区 $partition_id: 发现 $(length(cycles)) 个环路")
        catch e
            println("警告：分区 $partition_id 的环路检测出错: $(typeof(e)): $(e)")
            println(stacktrace())
        end
        
        cycles_by_partition[partition_id] = cycles
    end
    
    return cycles_by_partition
end



"""
    write_results_with_partitions(output_file, case::JuliaPowerCase, cycles_by_partition, virtual_connections, node_partitions, edge_partitions)

写入结果，包括分区信息。自动覆盖已存在的文件。
"""
function write_results_with_partitions(output_file, case::JuliaPowerCase, cycles_by_partition, virtual_connections, node_partitions, edge_partitions)
    # 创建结果数据框
    nodes_df = DataFrame(
        ID = Int[],
        Name = String[],
        Type = String[],
        Partition = Int[]
    )
    
    edges_df = DataFrame(
        From_ID = Int[],
        To_ID = Int[],
        From_Name = String[],
        To_Name = String[],
        Type = String[],
        Virtual = Bool[],
        Partition = Int[]
    )
    
    cycles_df = DataFrame(
        Partition = Int[],
        Cycle_ID = Int[],
        Nodes = String[]
    )
    
    # 填充节点数据
    for bus in case.busesAC
        push!(nodes_df, [
            bus.index,
            bus.name,
            "Bus",
            get(node_partitions, bus.index, 0)
        ])
    end
    
    # 填充边数据
    # 处理交流线路
    for line in case.branchesAC
        if line.in_service
            from_name = ""
            to_name = ""
            
            # 查找节点名称
            for bus in case.busesAC
                if bus.index == line.from_bus
                    from_name = bus.name
                end
                if bus.index == line.to_bus
                    to_name = bus.name
                end
            end
            
            is_virtual = (line.from_bus, line.to_bus) in virtual_connections || 
                         (line.to_bus, line.from_bus) in virtual_connections
            
            push!(edges_df, [
                line.from_bus,
                line.to_bus,
                from_name,
                to_name,
                "Line",
                is_virtual,
                get(edge_partitions, (line.from_bus, line.to_bus), 0)
            ])
        end
    end
    
    # 处理变压器
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            from_name = ""
            to_name = ""
            
            # 查找节点名称
            for bus in case.busesAC
                if bus.index == transformer.hv_bus
                    from_name = bus.name
                end
                if bus.index == transformer.lv_bus
                    to_name = bus.name
                end
            end
            
            is_virtual = (transformer.hv_bus, transformer.lv_bus) in virtual_connections || 
                         (transformer.lv_bus, transformer.hv_bus) in virtual_connections
            
            push!(edges_df, [
                transformer.hv_bus,
                transformer.lv_bus,
                from_name,
                to_name,
                "Transformer",
                is_virtual,
                get(edge_partitions, (transformer.hv_bus, transformer.lv_bus), 0)
            ])
        end
    end
    
    # 处理断路器
    for hvcb in case.hvcbs
        if hvcb.closed
            from_name = ""
            to_name = ""
            
            # 查找节点名称
            for bus in case.busesAC
                if bus.index == hvcb.bus_from
                    from_name = bus.name
                end
                if bus.index == hvcb.bus_to
                    to_name = bus.name
                end
            end
            
            is_virtual = (hvcb.bus_from, hvcb.bus_to) in virtual_connections || 
                         (hvcb.bus_to, hvcb.bus_from) in virtual_connections
            
            push!(edges_df, [
                hvcb.bus_from,
                hvcb.bus_to,
                from_name,
                to_name,
                "CircuitBreaker",
                is_virtual,
                get(edge_partitions, (hvcb.bus_from, hvcb.bus_to), 0)
            ])
        end
    end
    
    # 填充环路数据
    for (partition_id, cycles) in cycles_by_partition
        for (cycle_id, cycle) in enumerate(cycles)
            # 将环路中的节点ID转换为名称
            node_names = []
            for node_id in cycle
                for bus in case.busesAC
                    if bus.index == node_id
                        push!(node_names, bus.name)
                        break
                    end
                end
            end
            
            push!(cycles_df, [
                partition_id,
                cycle_id,
                join(node_names, " -> ") * " -> " * node_names[1]  # 闭合环路
            ])
        end
    end
    
    # 如果文件已存在，先删除它
    if isfile(output_file)
        rm(output_file)
    end
    
    # 创建一个新的XLSX文件
    XLSX.writetable(output_file, 
        Nodes = (collect(eachcol(nodes_df)), names(nodes_df)),
        Edges = (collect(eachcol(edges_df)), names(edges_df)),
        Cycles = (collect(eachcol(cycles_df)), names(cycles_df))
    )
    
    println("结果已保存至 $output_file")
    
    return Dict(
        "nodes" => nodes_df,
        "edges" => edges_df,
        "cycles" => cycles_df
    )
end

"""
    generate_partition_report(output_file, case::JuliaPowerCase, node_partitions, edge_partitions)

生成详细的分区报告。
"""
function generate_partition_report(output_file, case::JuliaPowerCase, node_partitions, edge_partitions)
    # 获取唯一的分区ID
    partition_ids = unique(values(node_partitions))
    
    # 创建分区统计数据框
    partition_stats = DataFrame(
        Partition_ID = Int[],
        Node_Count = Int[],
        Edge_Count = Int[],
        Load_Count = Int[],
        Total_Load_MW = Float64[],
        Has_Cycles = Bool[]
    )
    
    # 为每个分区计算统计信息
    for partition_id in partition_ids
        # 计算该分区的节点数
        nodes_in_partition = [node_id for (node_id, part_id) in node_partitions if part_id == partition_id]
        node_count = length(nodes_in_partition)
        
        # 计算该分区的边数
        edges_in_partition = [(from, to) for ((from, to), part_id) in edge_partitions if part_id == partition_id]
        edge_count = length(edges_in_partition)
        
        # 计算该分区的负荷数和总负荷
        load_count = 0
        total_load = 0.0
        
        for load in case.loadsAC
            if load.in_service && load.bus in nodes_in_partition
                load_count += 1
                total_load += load.p_mw
            end
        end
        
        # 检查该分区是否有环路
        has_cycles = edge_count >= node_count
        
        push!(partition_stats, [
            partition_id,
            node_count,
            edge_count,
            load_count,
            total_load,
            has_cycles
        ])
    end
    
    # 创建分区节点详情数据框
    partition_nodes = DataFrame(
        Partition_ID = Int[],
        Node_ID = Int[],
        Node_Name = String[],
        Node_Type = String[],
        Has_Load = Bool[],
        Load_MW = Float64[]
    )
    
    # 填充分区节点详情
    for (node_id, partition_id) in node_partitions
        node_name = ""
        node_type = ""
        has_load = false
        load_mw = 0.0
        
        # 查找节点名称和类型
        for bus in case.busesAC
            if bus.index == node_id
                node_name = bus.name
                node_type = "Bus"
                break
            end
        end
        
        # 检查该节点是否有负荷
        for load in case.loadsAC
            if load.in_service && load.bus == node_id
                has_load = true
                load_mw += load.p_mw
            end
        end
        
        push!(partition_nodes, [
            partition_id,
            node_id,
            node_name,
            node_type,
            has_load,
            load_mw
        ])
    end
    
    # 创建分区边详情数据框
    partition_edges = DataFrame(
        Partition_ID = Int[],
        From_ID = Int[],
        To_ID = Int[],
        From_Name = String[],
        To_Name = String[],
        Edge_Type = String[],
        Is_Virtual = Bool[]
    )
    
    # 填充分区边详情
    for ((from_id, to_id), partition_id) in edge_partitions
        from_name = ""
        to_name = ""
        edge_type = ""
        is_virtual = false
        
        # 查找节点名称
        for bus in case.busesAC
            if bus.index == from_id
                from_name = bus.name
            end
            if bus.index == to_id
                to_name = bus.name
            end
        end
        
        # 确定边的类型
        # 检查是否为线路
        for line in case.branchesAC
            if line.in_service && ((line.from_bus == from_id && line.to_bus == to_id) || 
                                   (line.from_bus == to_id && line.to_bus == from_id))
                edge_type = "Line"
                break
            end
        end
        
        # 检查是否为变压器
        if edge_type == ""
            for transformer in case.transformers_2w_etap
                if transformer.in_service && ((transformer.hv_bus == from_id && transformer.lv_bus == to_id) || 
                                             (transformer.hv_bus == to_id && transformer.lv_bus == from_id))
                    edge_type = "Transformer"
                    break
                end
            end
        end
        
        # 检查是否为断路器
        if edge_type == ""
            for hvcb in case.hvcbs
                if hvcb.closed && ((hvcb.bus_from == from_id && hvcb.bus_to == to_id) || 
                                  (hvcb.bus_from == to_id && hvcb.bus_to == from_id))
                    edge_type = "CircuitBreaker"
                    
                    # 检查是否为虚拟连接
                    if startswith(hvcb.name, "VIRTUAL_CB_")
                        is_virtual = true
                    end
                    
                    break
                end
            end
        end
        
        push!(partition_edges, [
            partition_id,
            from_id,
            to_id,
            from_name,
            to_name,
            edge_type,
            is_virtual
        ])
    end
    
    # 生成报告文件名
    report_file = replace(output_file, ".xlsx" => "_partition_report.xlsx")
    
    # 如果文件已存在，先删除它
    if isfile(report_file)
        rm(report_file)
    end
    
    # 将报告写入Excel文件
    XLSX.writetable(report_file, 
        Partition_Stats = (collect(eachcol(partition_stats)), names(partition_stats)),
        Partition_Nodes = (collect(eachcol(partition_nodes)), names(partition_nodes)),
        Partition_Edges = (collect(eachcol(partition_edges)), names(partition_edges))
    )
    
    println("分区报告已保存至 $report_file")
    
    return report_file
end


"""
    extract_edges_from_case(case::JuliaPowerCase)

从case中提取所有边的信息。
返回边的列表，每个边是一个元组(from_bus, to_bus, edge_object)。
"""
function extract_edges_from_case(case::JuliaPowerCase)
    edges = []
    
    # 处理交流线路
    for line in case.branchesAC
        if line.in_service
            try
                src = line.from_bus
                dst = line.to_bus
                push!(edges, (src, dst, line))
            catch
                try
                    src = line.hv_bus
                    dst = line.lv_bus
                    push!(edges, (src, dst, line))
                catch
                    src = line.bus_from
                    dst = line.bus_to
                    push!(edges, (src, dst, line))
                end
            end
        end
    end
    
    # 处理变压器
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            try
                src = transformer.from_bus
                dst = transformer.to_bus
                push!(edges, (src, dst, transformer))
            catch
                try
                    src = transformer.hv_bus
                    dst = transformer.lv_bus
                    push!(edges, (src, dst, transformer))
                catch
                    src = transformer.bus_from
                    dst = transformer.bus_to
                    push!(edges, (src, dst, transformer))
                end
            end
        end
    end
    
    # 处理断路器
    for hvcb in case.hvcbs
        if hvcb.closed
            try
                src = hvcb.from_bus
                dst = hvcb.to_bus
                push!(edges, (src, dst, hvcb))
            catch
                try
                    src = hvcb.hv_bus
                    dst = hvcb.lv_bus
                    push!(edges, (src, dst, hvcb))
                catch
                    src = hvcb.bus_from
                    dst = hvcb.bus_to
                    push!(edges, (src, dst, hvcb))
                end
            end
        end
    end
    
    return edges
end

"""
    topology_analysis(case::JuliaPowerCase; output_file = "./output_result.xlsx", debug=false)

主函数，执行完整的处理流程。
"""
function topology_analysis(case::JuliaPowerCase; output_file = "./output_result.xlsx", debug=false)
    # 创建节点映射
    node_dict = create_node_mapping(case)
    
    # 提取所有边
    edges = extract_edges_from_case(case)
    
    # 统计开关端口数
    port_counts = count_switch_ports(case)
    
    # 处理开关节点，添加虚拟节点，并获取虚拟连接标记
    new_case, virtual_connections = process_switch_nodes(case)
    
    # 识别分区
    node_partitions, edge_partitions = identify_partitions(new_case)
    
    # 创建和保存图形，同时获取环路信息，按分区分别绘制
    cycles = create_and_plot_graph_by_partition(new_case, node_partitions, edge_partitions, virtual_connections)
    
    # 写入结果
    results = write_results_with_partitions(output_file, new_case, cycles, virtual_connections, node_partitions, edge_partitions)
    
    # 生成详细的分区报告
    generate_partition_report(output_file, new_case, node_partitions, edge_partitions)
    
    return results, new_case
end

# 辅助函数：合并多个数据框为一个，添加"table"列
function gather_tables_as_one(tables::Vector{Pair{String, DataFrame}})::DataFrame
    # 收集所有可能的列
    allcols = Symbol[]
    for (_, df) in tables
        for c in names(df)
            if c ∉ allcols
                push!(allcols, c)
            end
        end
    end
    
    # 对于每个表，确保它有所有列，填充缺失值，添加"table"列
    combined = DataFrame()
    for (tname, df) in tables
        df_local = copy(df)
        for c in allcols
            if c ∉ names(df_local)
                df_local[!, c] = missing
            end
        end
        df_local[!, :table] = fill(tname, nrow(df_local))
        append!(combined, df_local)
    end
    return combined
end

