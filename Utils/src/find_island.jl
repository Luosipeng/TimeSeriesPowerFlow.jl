function find_islands(jpc::JPC)
    # 获取所有母线编号
    bus_ids = jpc.busAC[:, 1]
    n_bus = length(bus_ids)
    
    # 创建母线编号到索引的映射
    bus_id_to_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_ids)
        bus_id_to_idx[Int(bus_id)] = idx
    end
    
    # 获取母线类型索引
    PQ, PV, REF = idx_bus()[1:3]
    
    # 创建邻接矩阵
    adj = zeros(Int, n_bus, n_bus)
    
    # 只考虑状态为1的支路
    for i in 1:size(jpc.branchAC, 1)
        branch = jpc.branchAC[i, :]
        if branch[11] == 1  # BR_STATUS = 11
            f_bus_id = Int(branch[1])
            t_bus_id = Int(branch[2])
            
            # 使用映射获取正确的索引
            f_idx = bus_id_to_idx[f_bus_id]
            t_idx = bus_id_to_idx[t_bus_id]
            
            adj[f_idx, t_idx] = 1
            adj[t_idx, f_idx] = 1
        end
    end
    
    # 使用DFS找到连通分量
    visited = falses(n_bus)
    groups = Vector{Vector{Int}}()
    
    # 首先找到所有连通分量
    for i in 1:n_bus
        if !visited[i]
            component = Int[]
            stack = [i]
            visited[i] = true
            
            while !isempty(stack)
                node_idx = pop!(stack)
                node_id = Int(bus_ids[node_idx])  # 转换回原始母线编号
                push!(component, node_id)
                
                # 检查相邻节点
                for neighbor_idx in 1:n_bus
                    if adj[node_idx, neighbor_idx] == 1 && !visited[neighbor_idx]
                        push!(stack, neighbor_idx)
                        visited[neighbor_idx] = true
                    end
                end
            end
            
            push!(groups, sort(component))
        end
    end
    
    # 找出孤立节点和全PQ节点组
    isolated = Int[]
    groups_to_remove = Int[]
    
    # 检查每个组
    for (idx, group) in enumerate(groups)
        if length(group) == 1
            # 处理单节点组
            node_id = group[1]
            node_idx = bus_id_to_idx[node_id]
            connections = sum(adj[node_idx, :])
            
            # 找到对应的母线在jpc.busAC中的行
            bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
            bus_type = Int(jpc.busAC[bus_row, 2])
            
            if connections == 0 && bus_type == PQ
                # 孤立的PQ节点移到isolated
                push!(isolated, node_id)
                push!(groups_to_remove, idx)
            end
        else
            # 检查多节点组是否全是PQ节点
            all_pq = true
            for node_id in group
                # 找到对应的母线在jpc.busAC中的行
                bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
                bus_type = Int(jpc.busAC[bus_row, 2])
                
                if bus_type == PV || bus_type == REF
                    all_pq = false
                    break
                end
            end
            
            if all_pq
                # 如果全是PQ节点，将整个组移到isolated
                append!(isolated, group)
                push!(groups_to_remove, idx)
            end
        end
    end
    
    # 从后向前删除组，以避免索引变化带来的问题
    sort!(groups_to_remove, rev=true)
    for idx in groups_to_remove
        deleteat!(groups, idx)
    end
      
    # 返回结果
    return groups, isolated
end


function find_islands_acdc(jpc::JPC)
    # 获取所有AC和DC母线编号
    bus_ac_ids = jpc.busAC[:, 1]
    bus_dc_ids = size(jpc.busDC, 1) > 0 ? jpc.busDC[:, 1] : Int[]
    
    n_bus_ac = length(bus_ac_ids)
    n_bus_dc = length(bus_dc_ids)
    n_bus_total = n_bus_ac + n_bus_dc
    
    # 如果没有DC母线，直接使用原始的find_islands函数
    if n_bus_dc == 0
        return find_islands(jpc)
    end
    
    # 创建统一编号系统
    unified_ids = vcat(bus_ac_ids, bus_dc_ids)
    
    # 创建母线编号到索引的映射
    bus_id_to_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(unified_ids)
        bus_id_to_idx[Int(bus_id)] = idx
    end
    
    # 创建索引到原始系统（AC/DC）的映射
    idx_to_system = Dict{Int, Symbol}()
    for i in 1:n_bus_ac
        idx_to_system[i] = :AC
    end
    for i in (n_bus_ac+1):n_bus_total
        idx_to_system[i] = :DC
    end
    
    # 创建原始编号到统一索引的映射
    ac_id_to_unified_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_ac_ids)
        ac_id_to_unified_idx[Int(bus_id)] = idx
    end
    
    dc_id_to_unified_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_dc_ids)
        dc_id_to_unified_idx[Int(bus_id)] = n_bus_ac + idx
    end
    
    # 获取母线类型索引
    PQ, PV, REF = idx_bus()[1:3]
    
    # 创建邻接矩阵
    adj = zeros(Int, n_bus_total, n_bus_total)
    
    # 处理AC支路
    for i in 1:size(jpc.branchAC, 1)
        branch = jpc.branchAC[i, :]
        if branch[11] == 1  # BR_STATUS = 11
            f_bus_id = Int(branch[1])
            t_bus_id = Int(branch[2])
            
            # 使用映射获取正确的索引
            f_idx = ac_id_to_unified_idx[f_bus_id]
            t_idx = ac_id_to_unified_idx[t_bus_id]
            
            adj[f_idx, t_idx] = 1
            adj[t_idx, f_idx] = 1
        end
    end
    
    # 处理DC支路
    if size(jpc.branchDC, 1) > 0
        for i in 1:size(jpc.branchDC, 1)
            branch = jpc.branchDC[i, :]
            if branch[11] == 1  # 假设DC支路也有状态字段，位置与AC相同
                f_bus_id = Int(branch[1])
                t_bus_id = Int(branch[2])
                
                # 使用映射获取正确的索引
                if haskey(dc_id_to_unified_idx, f_bus_id) && haskey(dc_id_to_unified_idx, t_bus_id)
                    f_idx = dc_id_to_unified_idx[f_bus_id]
                    t_idx = dc_id_to_unified_idx[t_bus_id]
                    
                    adj[f_idx, t_idx] = 1
                    adj[t_idx, f_idx] = 1
                end
            end
        end
    end
    
    # 处理转换器（连接AC和DC系统）
    if size(jpc.converter, 1) > 0
        for i in 1:size(jpc.converter, 1)
            converter = jpc.converter[i, :]
            # 假设转换器默认是投运的，除非明确指定为0
            is_active = true
            if size(converter, 1) >= 3
                is_active = converter[3] != 0
            end
            
            if is_active
                ac_bus_id = Int(converter[1])
                dc_bus_id = Int(converter[2])
                
                # 检查这些母线是否存在
                if haskey(ac_id_to_unified_idx, ac_bus_id) && haskey(dc_id_to_unified_idx, dc_bus_id)
                    ac_idx = ac_id_to_unified_idx[ac_bus_id]
                    dc_idx = dc_id_to_unified_idx[dc_bus_id]
                    
                    adj[ac_idx, dc_idx] = 1
                    adj[dc_idx, ac_idx] = 1
                end
            end
        end
    end
    
    # 使用DFS找到连通分量
    visited = falses(n_bus_total)
    unified_groups = Vector{Vector{Int}}()
    
    # 首先找到所有连通分量（包括AC和DC节点）
    for i in 1:n_bus_total
        if !visited[i]
            component = Int[]
            stack = [i]
            visited[i] = true
            
            while !isempty(stack)
                node_idx = pop!(stack)
                node_id = Int(unified_ids[node_idx])  # 获取原始母线编号
                push!(component, node_idx)  # 存储索引而不是ID
                
                # 检查相邻节点
                for neighbor_idx in 1:n_bus_total
                    if adj[node_idx, neighbor_idx] == 1 && !visited[neighbor_idx]
                        push!(stack, neighbor_idx)
                        visited[neighbor_idx] = true
                    end
                end
            end
            
            push!(unified_groups, component)
        end
    end
    
    # 将统一索引转换回原始母线编号，并分离AC和DC部分
    groups = Vector{Vector{Int}}()
    
    for group in unified_groups
        ac_component = Int[]
        
        for node_idx in group
            if idx_to_system[node_idx] == :AC
                node_id = Int(unified_ids[node_idx])
                push!(ac_component, node_id)
            end
        end
        
        if !isempty(ac_component)
            push!(groups, sort(ac_component))
        end
    end
    
    # 找出孤立节点和全PQ节点组
    isolated = Int[]
    groups_to_remove = Int[]
    
    # 检查每个组
    for (idx, group) in enumerate(groups)
        if length(group) == 1
            # 处理单节点组
            node_id = group[1]
            node_idx = ac_id_to_unified_idx[node_id]
            
            # 检查该节点是否有连接（包括到DC系统的连接）
            connections = sum(adj[node_idx, :])
            
            # 找到对应的母线在jpc.busAC中的行
            bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, 2])
                
                if connections == 0 && bus_type == PQ
                    # 孤立的PQ节点移到isolated
                    push!(isolated, node_id)
                    push!(groups_to_remove, idx)
                end
            end
        else
            # 检查多节点组是否全是PQ节点
            all_pq = true
            has_generator = false
            
            for node_id in group
                # 找到对应的母线在jpc.busAC中的行
                bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
                if bus_row !== nothing
                    bus_type = Int(jpc.busAC[bus_row, 2])
                    
                    if bus_type == PV || bus_type == REF
                        all_pq = false
                        break
                    end
                end
                
                # 检查是否有连接到该母线的发电机
                for j in 1:size(jpc.genAC, 1)
                    gen_bus = Int(jpc.genAC[j, 1])  # 假设第二列是母线编号
                    if gen_bus == node_id && jpc.genAC[j, 8] == 1  # 假设第8列是状态
                        has_generator = true
                        all_pq = false
                        break
                    end
                end
                
                if !all_pq
                    break
                end
            end
            
            # 找到与该组相连的DC母线
            connected_dc_buses = Set{Int}()
            for ac_bus_id in group
                if haskey(ac_id_to_unified_idx, ac_bus_id)
                    ac_idx = ac_id_to_unified_idx[ac_bus_id]
                    
                    # 查找与此AC节点直接相连的DC节点
                    for dc_idx in (n_bus_ac+1):n_bus_total
                        if adj[ac_idx, dc_idx] == 1
                            dc_bus_id = unified_ids[dc_idx]
                            push!(connected_dc_buses, dc_bus_id)
                        end
                    end
                end
            end
            
            # 检查这些DC母线是否有发电机
            if !isempty(connected_dc_buses) && size(jpc.sgenDC, 1) > 0
                for j in 1:size(jpc.sgenDC, 1)
                    sgen_bus = Int(jpc.sgenDC[j, 1])  # 假设第一列是母线编号
                    # 假设sgenDC的第3列是状态字段，如果不存在则假设为投运
                    is_active = size(jpc.sgenDC, 2) >= 3 ? jpc.sgenDC[j, 3] != 0 : true
                    
                    if sgen_bus in connected_dc_buses && is_active
                        has_generator = true
                        all_pq = false
                        break
                    end
                end
            end
            
            if all_pq && !has_generator
                # 如果全是PQ节点且没有发电机，将整个组移到isolated
                append!(isolated, group)
                push!(groups_to_remove, idx)
            end
        end
    end
    
    # 从后向前删除组，以避免索引变化带来的问题
    sort!(groups_to_remove, rev=true)
    for idx in groups_to_remove
        deleteat!(groups, idx)
    end
      
    # 返回结果
    return groups, isolated
end
