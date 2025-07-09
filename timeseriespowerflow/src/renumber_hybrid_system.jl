function renumber_hybrid_system(jpc)
    # 提取必要的数据
    busAC = jpc.busAC
    busDC = jpc.busDC
    branchAC = jpc.branchAC
    branchDC = jpc.branchDC
    converter = jpc.converter
    
    # 提取负荷和发电机数据（如果存在）
    loadAC = jpc.loadAC
    loadDC = jpc.loadDC
    genAC = jpc.genAC
    genDC = jpc.genDC
    pvarray = jpc.pv
    
    # 提取交流光伏系统数据（如果存在）
    pv_acsystem = hasfield(typeof(jpc), :pv_acsystem) ? jpc.pv_acsystem : nothing
    
    # 提取储能数据（如果存在）
    storage = hasfield(typeof(jpc), :storage) ? jpc.storage : nothing

    # 为了区分交流和直流节点，我们创建唯一标识符
    # 使用元组 (bus_number, is_ac) 来表示节点
    # is_ac 为 true 表示交流节点，false 表示直流节点
    ac_buses = [(bus, true) for bus in busAC[:, 1]]
    dc_buses = [(bus, false) for bus in busDC[:, 1]]
    all_buses = vcat(ac_buses, dc_buses)
    
    # 创建节点映射（旧编号 -> 新编号）
    # 键是 (bus_number, is_ac)，值是新的编号
    old_to_new = Dict{Tuple{Float64, Bool}, Int}()
    
    # 创建邻接表表示图
    graph = Dict{Tuple{Float64, Bool}, Vector{Tuple{Float64, Bool}}}()
    
    # 初始化图
    for bus in all_buses
        graph[bus] = Tuple{Float64, Bool}[]
    end
    
    # 添加交流支路连接
    for i in 1:size(branchAC, 1)
        from_bus = (branchAC[i, 1], true)  # 交流节点
        to_bus = (branchAC[i, 2], true)    # 交流节点
        push!(graph[from_bus], to_bus)
        push!(graph[to_bus], from_bus)
    end
    
    # 添加直流支路连接
    for i in 1:size(branchDC, 1)
        from_bus = (branchDC[i, 1], false)  # 直流节点
        to_bus = (branchDC[i, 2], false)    # 直流节点
        push!(graph[from_bus], to_bus)
        push!(graph[to_bus], from_bus)
    end
    
    # 添加变流器连接（连接交流和直流节点）
    for i in 1:size(converter, 1)
        ac_bus = (converter[i, 1], true)   # 交流节点
        dc_bus = (converter[i, 2], false)  # 直流节点
        push!(graph[ac_bus], dc_bus)
        push!(graph[dc_bus], ac_bus)
    end
    
    # 找出交流系统的平衡节点（假设类型为3的节点是平衡节点）
    slack_bus_idx = findfirst(x -> x == 3.0, busAC[:, 2])
    if slack_bus_idx === nothing
        slack_bus = (busAC[1, 1], true)  # 如果没有找到平衡节点，使用第一个交流节点
    else
        slack_bus = (busAC[slack_bus_idx, 1], true)
    end
    
    # 使用BFS进行重新编号（将整个系统视为一个连通系统）
    visited = Set{Tuple{Float64, Bool}}()
    queue = [slack_bus]  # 从平衡节点开始BFS
    new_number = 1
    
    while !isempty(queue)
        current = popfirst!(queue)
        
        if current in visited
            continue
        end
        
        push!(visited, current)
        old_to_new[current] = new_number
        new_number += 1
        
        # 将所有未访问的邻居节点加入队列
        for neighbor in graph[current]
            if !(neighbor in visited)
                push!(queue, neighbor)
            end
        end
    end
    
    # 处理可能存在的未连接到主网络的节点
    for bus in all_buses
        if !(bus in keys(old_to_new))
            old_to_new[bus] = new_number
            new_number += 1
        end
    end
    
    # 创建新的数据结构
    new_busAC = copy(busAC)
    new_busDC = copy(busDC)
    new_branchAC = copy(branchAC)
    new_branchDC = copy(branchDC)
    new_converter = copy(converter)
    
    # 复制负荷和发电机数据（如果存在）
    new_loadAC = loadAC !== nothing ? copy(loadAC) : nothing
    new_loadDC = loadDC !== nothing ? copy(loadDC) : nothing
    new_genAC = genAC !== nothing ? copy(genAC) : nothing
    new_genDC = genDC !== nothing ? copy(genDC) : nothing
    new_pvarray = pvarray !== nothing ? copy(pvarray) : nothing
    
    # 复制交流光伏系统数据（如果存在）
    new_pv_acsystem = pv_acsystem !== nothing ? copy(pv_acsystem) : nothing
    
    # 复制储能数据（如果存在）
    new_storage = storage !== nothing ? copy(storage) : nothing
    
    # 更新交流母线编号
    for i in 1:size(busAC, 1)
        old_num = busAC[i, 1]
        new_busAC[i, 1] = old_to_new[(old_num, true)]
    end
    
    # 更新直流母线编号
    for i in 1:size(busDC, 1)
        old_num = busDC[i, 1]
        new_busDC[i, 1] = old_to_new[(old_num, false)]
    end
    
    # 更新交流支路编号并确保from < to
    for i in 1:size(branchAC, 1)
        from_old = branchAC[i, 1]
        to_old = branchAC[i, 2]
        from_new = old_to_new[(from_old, true)]
        to_new = old_to_new[(to_old, true)]
        
        # 确保小的节点号在前
        if from_new < to_new
            new_branchAC[i, 1] = from_new
            new_branchAC[i, 2] = to_new
        else
            new_branchAC[i, 1] = to_new
            new_branchAC[i, 2] = from_new
        end
    end
    
    # 更新直流支路编号并确保from < to
    for i in 1:size(branchDC, 1)
        from_old = branchDC[i, 1]
        to_old = branchDC[i, 2]
        from_new = old_to_new[(from_old, false)]
        to_new = old_to_new[(to_old, false)]
        
        # 确保小的节点号在前
        if from_new < to_new
            new_branchDC[i, 1] = from_new
            new_branchDC[i, 2] = to_new
        else
            new_branchDC[i, 1] = to_new
            new_branchDC[i, 2] = from_new
        end
    end
    
    # 更新变流器编号并确保from < to
    for i in 1:size(converter, 1)
        ac_old = converter[i, 1]
        dc_old = converter[i, 2]
        ac_new = old_to_new[(ac_old, true)]
        dc_new = old_to_new[(dc_old, false)]
        
        # 确保小的节点号在前
        if ac_new < dc_new
            new_converter[i, 1] = ac_new
            new_converter[i, 2] = dc_new
        else
            new_converter[i, 1] = dc_new
            new_converter[i, 2] = ac_new
        end
    end
    
    # 更新交流负荷编号（如果存在）
    if new_loadAC !== nothing
        for i in 1:size(new_loadAC, 1)
            bus_num = new_loadAC[i, 2]  # 第二列是母线编号
            if haskey(old_to_new, (bus_num, true))
                new_loadAC[i, 2] = old_to_new[(bus_num, true)]
            end
        end
    end
    
    # 更新直流负荷编号（如果存在）
    if new_loadDC !== nothing
        for i in 1:size(new_loadDC, 1)
            bus_num = new_loadDC[i, 2]  # 第二列是母线编号
            if haskey(old_to_new, (bus_num, false))
                new_loadDC[i, 2] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # 更新交流发电机编号（如果存在）
    if new_genAC !== nothing
        for i in 1:size(new_genAC, 1)
            bus_num = new_genAC[i, 1]  # 第一列是母线编号
            if haskey(old_to_new, (bus_num, true))
                new_genAC[i, 1] = old_to_new[(bus_num, true)]
            end
        end
    end
    
    # 更新直流发电机编号（如果存在）
    if new_genDC !== nothing
        for i in 1:size(new_genDC, 1)
            bus_num = new_genDC[i, 1]  # 第一列是母线编号
            if haskey(old_to_new, (bus_num, false))
                new_genDC[i, 1] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # 更新光伏阵列编号（如果存在）
    if new_pvarray !== nothing
        for i in 1:size(new_pvarray, 1)
            bus_num = new_pvarray[i, 2]  # 第二列是母线编号
            # 假设光伏阵列连接到直流节点
            if haskey(old_to_new, (bus_num, false))
                new_pvarray[i, 2] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # 更新交流光伏系统编号（如果存在）
    if new_pv_acsystem !== nothing
        for i in 1:size(new_pv_acsystem, 1)
            bus_num = new_pv_acsystem[i, 2]  # 第二列是母线编号
            # 交流光伏系统连接到交流节点
            if haskey(old_to_new, (bus_num, true))
                new_pv_acsystem[i, 2] = old_to_new[(bus_num, true)]
            else
                @warn "找不到ID为 $bus_num 的交流母线，无法更新交流光伏系统编号"
            end
        end
    end
    
    # 更新储能设备编号（如果存在）
    if new_storage !== nothing
        for i in 1:size(new_storage, 1)
            bus_num = new_storage[i, 1]  # 假设第一列是母线编号 (ESS_BUS)
            # 确定储能设备连接到的是交流还是直流母线
            # 这里需要根据实际情况确定，我假设储能设备可能连接到交流或直流母线
            if haskey(old_to_new, (bus_num, false))
                # 连接到直流母线
                new_storage[i, 1] = old_to_new[(bus_num, false)]
            else
                @warn "找不到ID为 $bus_num 的母线，无法更新储能设备编号"
            end
        end
    end
    
    # 对结果进行排序
    sort_idx_ac = sortperm(new_busAC[:, 1])
    sort_idx_dc = sortperm(new_busDC[:, 1])
    new_busAC = new_busAC[sort_idx_ac, :]
    new_busDC = new_busDC[sort_idx_dc, :]
    
    # 对支路进行排序（按照起始节点排序）
    sort_idx_branchAC = sortperm(new_branchAC[:, 1])
    sort_idx_branchDC = sortperm(new_branchDC[:, 1])
    sort_idx_converter = sortperm(new_converter[:, 1])
    new_branchAC = new_branchAC[sort_idx_branchAC, :]
    new_branchDC = new_branchDC[sort_idx_branchDC, :]
    new_converter = new_converter[sort_idx_converter, :]
    
    # 对负荷和发电机数据进行排序（如果存在）
    if new_loadAC !== nothing
        sort_idx_loadAC = sortperm(new_loadAC[:, 1])  # 按照母线编号排序
        new_loadAC = new_loadAC[sort_idx_loadAC, :]
    end
    
    if new_loadDC !== nothing
        sort_idx_loadDC = sortperm(new_loadDC[:, 1])  # 按照母线编号排序
        new_loadDC = new_loadDC[sort_idx_loadDC, :]
    end
    
    if new_genAC !== nothing
        sort_idx_genAC = sortperm(new_genAC[:, 1])  # 按照母线编号排序
        new_genAC = new_genAC[sort_idx_genAC, :]
    end
    
    if new_genDC !== nothing
        sort_idx_genDC = sortperm(new_genDC[:, 1])  # 按照母线编号排序
        new_genDC = new_genDC[sort_idx_genDC, :]
    end
    
    if new_pvarray !== nothing
        sort_idx_pvarray = sortperm(new_pvarray[:, 2])  # 按照母线编号排序
        new_pvarray = new_pvarray[sort_idx_pvarray, :]
    end
    
    # 对交流光伏系统数据进行排序（如果存在）
    if new_pv_acsystem !== nothing
        sort_idx_pv_acsystem = sortperm(new_pv_acsystem[:, 2])  # 按照母线编号排序
        new_pv_acsystem = new_pv_acsystem[sort_idx_pv_acsystem, :]
    end
    
    # 对储能设备数据进行排序（如果存在）
    if new_storage !== nothing
        sort_idx_storage = sortperm(new_storage[:, 1])  # 按照母线编号排序
        new_storage = new_storage[sort_idx_storage, :]
    end
    
    # 创建新的jpc结构
    new_jpc = (
        busAC = new_busAC,
        busDC = new_busDC,
        branchAC = new_branchAC,
        branchDC = new_branchDC,
        converter = new_converter
    )
    
    # 创建完整的映射字典（便于查询）
    complete_mapping = Dict{Tuple{Float64, String}, Int}()
    for ((bus, is_ac), new_num) in old_to_new
        system_type = is_ac ? "AC" : "DC"
        complete_mapping[(bus, system_type)] = new_num
    end
    
    # 添加映射字典到结果中
    new_jpc = merge(new_jpc, (mapping = complete_mapping,))
    
    # 添加负荷和发电机数据（如果存在）
    if loadAC !== nothing
        new_jpc = merge(new_jpc, (loadAC = new_loadAC,))
    end
    
    if loadDC !== nothing
        new_jpc = merge(new_jpc, (loadDC = new_loadDC,))
    end
    
    if genAC !== nothing
        new_jpc = merge(new_jpc, (genAC = new_genAC,))
    end
    
    if genDC !== nothing
        new_jpc = merge(new_jpc, (genDC = new_genDC,))
    end
    
    if pvarray !== nothing
        new_jpc = merge(new_jpc, (pv = new_pvarray,))
    end
    
    # 添加交流光伏系统数据（如果存在）
    if pv_acsystem !== nothing
        new_jpc = merge(new_jpc, (pv_acsystem = new_pv_acsystem,))
    end
    
    # 添加储能数据（如果存在）
    if storage !== nothing
        new_jpc = merge(new_jpc, (storage = new_storage,))
    end
    
    # 创建交流节点和直流节点的映射字典
    ac_node_mapping = Dict{Float64, Int}()
    dc_node_mapping = Dict{Float64, Int}()
    
    # 填充映射字典
    for ((bus, is_ac), new_num) in old_to_new
        if is_ac
            ac_node_mapping[bus] = new_num
        else
            dc_node_mapping[bus] = new_num
        end
    end
    
    return new_jpc, ac_node_mapping, dc_node_mapping
end
