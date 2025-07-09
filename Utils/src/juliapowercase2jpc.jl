"""
    resolve_node_mapping(node_id, node_merge_map)

递归解析节点映射，确保多层合并正确处理。
"""
function resolve_node_mapping(node_id, node_merge_map)
    while haskey(node_merge_map, node_id)
        node_id = node_merge_map[node_id]
    end
    return node_id
end

"""
    merge_virtual_nodes(case::JuliaPowerCase)

合并虚拟节点两侧的节点，去除虚拟节点和虚拟连接，避免在潮流计算中出现奇异导纳矩阵。
返回更新后的case对象，不含虚拟节点和虚拟连接。
"""
function merge_virtual_nodes(case::JuliaPowerCase)
    # 深拷贝case以便修改
    new_case = deepcopy(case)
    
    # 识别虚拟节点（名称中包含"_虚拟节点"的节点）
    virtual_node_ids = Int[]
    virtual_node_map = Dict{Int, String}()  # 虚拟节点ID到名称的映射
    
    for bus in new_case.busesAC
        if occursin("_虚拟节点", bus.name)
            push!(virtual_node_ids, bus.bus_id)
            virtual_node_map[bus.bus_id] = bus.name
        end
    end
    
    # 如果没有虚拟节点，直接返回原始case
    if isempty(virtual_node_ids)
        return new_case
    end
    
    # 识别与虚拟节点相连的断路器和实际节点
    virtual_connections = Dict{Int, Vector{Tuple{Int, Int}}}()  # 虚拟节点ID -> [(连接节点ID, 断路器索引)]
    
    for (i, hvcb) in enumerate(new_case.hvcbs)
        if !hvcb.closed || !hvcb.in_service
            continue
        end
        
        # 检查断路器是否连接到虚拟节点
        if hvcb.bus_from in virtual_node_ids
            if !haskey(virtual_connections, hvcb.bus_from)
                virtual_connections[hvcb.bus_from] = Tuple{Int, Int}[]
            end
            push!(virtual_connections[hvcb.bus_from], (hvcb.bus_to, i))
        elseif hvcb.bus_to in virtual_node_ids
            if !haskey(virtual_connections, hvcb.bus_to)
                virtual_connections[hvcb.bus_to] = Tuple{Int, Int}[]
            end
            push!(virtual_connections[hvcb.bus_to], (hvcb.bus_from, i))
        end
    end
    
    # 为每个虚拟节点，确定合并策略
    node_merge_map = Dict{Int, Int}()  # 原节点ID -> 合并后的节点ID
    
    # 第一步：处理每个虚拟节点，初步确定合并映射
    for virtual_node_id in virtual_node_ids
        if !haskey(virtual_connections, virtual_node_id) || length(virtual_connections[virtual_node_id]) < 2
            # 虚拟节点至少需要连接两个其他节点才能进行合并
            continue
        end
        
        # 获取与虚拟节点相连的所有实际节点
        connected_nodes = [node_id for (node_id, _) in virtual_connections[virtual_node_id]]
        
        # 选择第一个节点作为合并的目标节点
        target_node_id = connected_nodes[1]
        
        # 将其他节点映射到目标节点
        for node_id in connected_nodes[2:end]
            # 检查是否已经有映射，如果有，需要确保一致性
            if haskey(node_merge_map, node_id)
                # 已经有映射，需要将当前的目标节点也映射到同一个节点
                existing_target = resolve_node_mapping(node_id, node_merge_map)
                node_merge_map[target_node_id] = existing_target
                target_node_id = existing_target
            else
                node_merge_map[node_id] = target_node_id
            end
        end
        
        # 将虚拟节点也映射到目标节点
        node_merge_map[virtual_node_id] = target_node_id
    end
    
    # 第二步：解析所有映射，确保映射的一致性
    for node_id in keys(node_merge_map)
        node_merge_map[node_id] = resolve_node_mapping(node_merge_map[node_id], node_merge_map)
    end
    
    # 更新所有元素中的节点引用
    
    # 1. 更新交流线路
    for line in new_case.branchesAC
        if !line.in_service
            continue
        end
        
        line.from_bus = resolve_node_mapping(line.from_bus, node_merge_map)
        line.to_bus = resolve_node_mapping(line.to_bus, node_merge_map)
        
        # 检查是否形成自环（同一节点连接到自身）
        if line.from_bus == line.to_bus
            line.in_service = false  # 禁用自环
        end
    end
    
    # 2. 更新变压器
    for transformer in new_case.transformers_2w_etap
        if !transformer.in_service
            continue
        end
        
        transformer.hv_bus = resolve_node_mapping(transformer.hv_bus, node_merge_map)
        transformer.lv_bus = resolve_node_mapping(transformer.lv_bus, node_merge_map)
        
        # 检查是否形成自环
        if transformer.hv_bus == transformer.lv_bus
            transformer.in_service = false  # 禁用自环
        end
    end
    
    # 3. 更新断路器
    new_hvcbs = []
    for hvcb in new_case.hvcbs
        if !hvcb.closed || !hvcb.in_service
            push!(new_hvcbs, hvcb)
            continue
        end
        
        # 解析节点映射
        new_from = resolve_node_mapping(hvcb.bus_from, node_merge_map)
        new_to = resolve_node_mapping(hvcb.bus_to, node_merge_map)
        
        # 如果断路器不形成自环且不连接到将被移除的虚拟节点，则保留
        if new_from != new_to && !(hvcb.bus_from in virtual_node_ids || hvcb.bus_to in virtual_node_ids)
            hvcb.bus_from = new_from
            hvcb.bus_to = new_to
            push!(new_hvcbs, hvcb)
        end
    end
    new_case.hvcbs = new_hvcbs
    
    # 4. 更新负荷
    for load in new_case.loadsAC
        if !load.in_service
            continue
        end
        
        load.bus = resolve_node_mapping(load.bus, node_merge_map)
    end
    
    # 5. 更新发电机
    for gen in new_case.sgensAC
        if !gen.in_service
            continue
        end
        
        gen.bus = resolve_node_mapping(gen.bus, node_merge_map)
    end
    
    # 6. 更新external grids
    for ext in new_case.ext_grids
        if !ext.in_service
            continue
        end
        
        ext.bus = resolve_node_mapping(ext.bus, node_merge_map)
    end
    
    # 移除虚拟节点
    new_busesAC = []
    for bus in new_case.busesAC
        if !(bus.bus_id in virtual_node_ids)
            push!(new_busesAC, bus)
        end
    end
    new_case.busesAC = new_busesAC
    
    # 更新节点名称到ID的映射
    new_case.bus_name_to_id = Dict{String, Int}()
    for bus in new_case.busesAC
        new_case.bus_name_to_id[bus.name] = bus.bus_id
    end
    
    # 合并同一节点上的负荷
    load_by_bus = Dict{Int, Vector{Int}}()  # 节点ID -> 负荷索引列表
    
    for (i, load) in enumerate(new_case.loadsAC)
        if !load.in_service
            continue
        end
        
        if !haskey(load_by_bus, load.bus)
            load_by_bus[load.bus] = Int[]
        end
        push!(load_by_bus[load.bus], i)
    end
    
    # 对每个有多个负荷的节点，合并负荷
    new_loadsAC = []
    processed_loads = Set{Int}()
    
    for (bus_id, load_indices) in load_by_bus
        if length(load_indices) <= 1
            # 只有一个负荷，直接保留
            for idx in load_indices
                push!(new_loadsAC, new_case.loadsAC[idx])
                push!(processed_loads, idx)
            end
        else
            # 多个负荷，合并
            total_p_mw = 0.0
            total_q_mvar = 0.0
            base_load = nothing
            
            for idx in load_indices
                load = new_case.loadsAC[idx]
                total_p_mw += load.p_mw
                total_q_mvar += load.q_mvar
                
                if base_load === nothing
                    base_load = deepcopy(load)
                end
                
                push!(processed_loads, idx)
            end
            
            # 使用第一个负荷作为基础，更新有功和无功功率
            if base_load !== nothing
                base_load.p_mw = total_p_mw
                base_load.q_mvar = total_q_mvar
                base_load.name = "合并负荷_$(bus_id)"
                push!(new_loadsAC, base_load)
            end
        end
    end
    
    # 添加未处理的负荷（例如已禁用的负荷）
    for (i, load) in enumerate(new_case.loadsAC)
        if i ∉ processed_loads && !load.in_service
            push!(new_loadsAC, load)
        end
    end
    
    new_case.loadsAC = new_loadsAC
    
    # 类似地，可以合并同一节点上的其他元素（如发电机、分路元件等）
    # ...
    
    return new_case
end

function JuliaPowerCase2Jpc(case::Utils.JuliaPowerCase)
    # 1. 合并虚拟节点
    case = merge_virtual_nodes(case)
    
    # 2. 创建JPC对象
    jpc = JPC()
    
    # 3. 设置基本参数
    jpc.baseMVA = case.baseMVA
    
    # 4. 设置节点数据
    JPC_buses_process(case, jpc)

    # 5. 设置直流节点数据
    JPC_dcbuses_process(case, jpc)
    
    # 6. 设置线路数据
    JPC_branches_process(case, jpc)
    
    # 7. 设置直流线路数据
    JPC_dcbranches_process(case, jpc)

    # 8. 设置发电机数据
    JPC_gens_process(case, jpc)

    # 9. 设置直流发电机数据
    JPC_battery_gens_process(case, jpc)

    # 10. 设置soc电池数据
    JPC_battery_soc_process(case, jpc)
    
    # 11. 设置负荷数据
    JPC_loads_process(case, jpc)

    # 12. 设置直流负荷数据
    JPC_dcloads_process(case, jpc)

    # 13. 设置PV阵列数据 
    JPC_pv_process(case, jpc)

    # 14. 设置换流器数据
    JPC_inverters_process(case, jpc)

    # 15. 设置交流光伏系统数据
    JPC_ac_pv_system_process(case, jpc)

    return jpc

end

function JPC_buses_process(case::JuliaPowerCase, jpc::JPC)
    # 获取节点数据并深拷贝防止误操作
    buses = deepcopy(case.busesAC)
    
    # 创建一个空矩阵，行数为节点数，列数为13
    num_buses = length(buses)
    bus_matrix = zeros(num_buses, 13)
    
    for (i, bus) in enumerate(buses)
        # 设置电压初始值（根据序号）
        vm = 1.0
        va = 0.0
        
        # 填充矩阵的每一行
        bus_matrix[i, :] = [
            bus.bus_id,      # 节点ID
            1.0,             # 节点类型(全部设为PQ节点)
            0.0,             # PD (MW) 有功负荷（MW）
            0.0,             # QD (MVAR) 无功负荷（MVAR）
            0.0,             # GS (MW) 有功发电（MW）
            0.0,             # BS (MVAR) 无功发电（MVAR）
            bus.area_id,     # 区域编号
            vm,              # 节点电压幅值（p.u.）
            va,              # 节点电压相角（度）
            bus.vn_kv,       # 节点电压基准（kV）
            bus.zone_id,     # 区域编号
            bus.max_vm_pu,   # 最大电压幅值（p.u.）
            bus.min_vm_pu,   # 最小电压幅值（p.u.）
        ]
    end
    
    # 所有序号的结果都存储到busAC字段
    jpc.busAC = bus_matrix
    
    return jpc
end

function JPC_dcbuses_process(case, jpc)
    # 处理直流节点数据，转换为JPC格式
    dcbuses = deepcopy(case.busesDC)
    
    # 创建一个空矩阵，行数为节点数，列数为13
    num_dcbuses = length(dcbuses)
    dcbus_matrix = zeros(num_dcbuses, 13)
    
    for (i, dcbus) in enumerate(dcbuses)
        # 设置电压初始值（根据序号）
        vm = 1.0
        va = 0.0
        
        # 填充矩阵的每一行
        dcbus_matrix[i, :] = [
            dcbus.bus_id,      # 节点ID
            1.0,               # 节点类型(全部设为PQ节点)
            0.0,               # PD (MW) 有功负荷（MW）
            0.0,               # QD (MVAR) 无功负荷（MVAR）
            0.0,               # GS (MW) 有功发电（MW）
            0.0,               # BS (MVAR) 无功发电（MVAR）
            dcbus.area_id,     # 区域编号
            vm,                # 节点电压幅值（p.u.）
            va,                # 节点电压相角（度）
            dcbus.vn_kv,       # 节点电压基准（kV）
            dcbus.zone_id,     # 区域编号
            dcbus.max_vm_pu,   # 最大电压幅值（p.u.）
            dcbus.min_vm_pu,   # 最小电压幅值（p.u.）
        ]
    end
    
    # 所有序号的结果都存储到busDC字段
    jpc.busDC = dcbus_matrix

    jpc = JPC_battery_bus_process(case, jpc)
    
    return jpc
end

function JPC_battery_bus_process(case::JuliaPowerCase, jpc::JPC)
    # 处理电池节点数据，转换为JPC格式并合并到busDC
    batteries = deepcopy(case.storageetap)
    
    # 获取当前busDC的数据
    busDC = jpc.busDC
    current_size = size(busDC, 1)
    
    # 创建一个矩阵来存储电池节点数据
    num_batteries = length(batteries)
    battery_matrix = zeros(num_batteries, size(busDC, 2))
    
    for (i, battery) in enumerate(batteries)
        # 设置电压初始值
        vm = 1.0
        va = 0.0
        vn_kv = battery.voc
        # 创建虚拟节点数据
        battery_row = zeros(1, size(busDC, 2))
        battery_row[1, :] = [
            current_size + i,    # 为电池分配新的节点ID
            2.0,                 # 节点类型(全部设为slack节点)
            0.0,                 # PD (MW) 有功负荷（MW）
            0.0,                 # QD (MVAR) 无功负荷（MVAR）
            0.0,                 # GS (MW) 有功发电（MW）
            0.0,                 # BS (MVAR) 无功发电（MVAR）
            1.0,                 # 区域编号
            vm,                  # 节点电压幅值（p.u.）
            va,                  # 节点电压相角（度）
            vn_kv,               # 节点电压基准（kV）
            1.0,                 # 区域编号
            1.05,                 # 最大电压幅值（p.u.）
            0.95                  # 最小电压幅值（p.u.）
        ]
        
        # 将电池节点数据存入矩阵
        battery_matrix[i, :] = battery_row
    end
    
    # 将电池虚拟节点合并到busDC中
    jpc.busDC = vcat(busDC, battery_matrix)
    
    # # 同时保存原始电池数据到jpc.battery字段，以便后续处理
    # jpc.battery = battery_matrix

    return jpc
end

function JPC_branches_process(case::JuliaPowerCase, jpc::JPC)
    # if sequence == 1||sequence == 2
        # 处理线路数据，转换为JPC格式
        calculate_line_parameters(case::JuliaPowerCase, jpc)
        # 处理变压器数据，转换为JPC格式
        calculate_transformer2w_parameters(case::JuliaPowerCase, jpc)
        # 处理三相变压器数据，转换为JPC格式
    # else
    #     # 处理支路数据，转换为JPC格式
    #     calculate_branch_JPC_zero(case::JuliaPowerCase, jpc)
    # end
end

function JPC_dcbranches_process(case::JuliaPowerCase, jpc::JPC)
    # 处理直流线路数据，转换为JPC格式
    nbr = length(case.branchesDC)
    branch = zeros(nbr, 14)
    dclines = case.branchesDC

    for (i, dcline) in enumerate(dclines)
        # 获取起始和终止母线编号
        from_bus_idx = dcline.from_bus
        to_bus_idx = dcline.to_bus
        
        # 获取起始母线的基准电压(kV)
        basekv = jpc.busDC[from_bus_idx, BASE_KV]
        
        # 计算基准阻抗
        baseR = (basekv^2) / case.baseMVA
        
        # 计算标幺值阻抗
        r_pu = 2 * dcline.length_km * dcline.r_ohm_per_km / baseR
        x_pu = 0
        
        # 填充branchAC矩阵
        branch[i, F_BUS] = from_bus_idx
        branch[i, T_BUS] = to_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        
        # 设置额定容量
        if hasfield(typeof(dcline), :max_i_ka)
            branch[i, RATE_A] = dcline.max_i_ka * basekv * sqrt(3)  # 额定容量(MVA)
        else
            branch[i, RATE_A] = 100.0  # 默认值
        end
        
        # 设置支路状态
        branch[i, BR_STATUS] = dcline.in_service ? 1.0 : 0.0
        
        # 设置相角限制
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end
    # 将直流线路数据添加到JPC结构体
    if isempty(jpc.branchDC)
        jpc.branchDC = branch
    else
        jpc.branchDC = [jpc.branchDC; branch]
    end
    
    # 处理储能虚拟连接
    jpc = JPC_battery_branch_process(case, jpc)

    return jpc
end

function JPC_battery_branch_process(case::JuliaPowerCase, jpc::JPC)
    # 处理电池虚拟连接，创建从电池虚拟节点到实际连接节点的支路
    batteries = deepcopy(case.storageetap)
    num_batteries = length(batteries)
    
    # 如果没有电池，直接返回
    if num_batteries == 0
        return jpc
    end
    
    # 创建电池虚拟支路矩阵，与branchDC结构相同
    battery_branches = zeros(num_batteries, 14)
    
    # 获取当前busDC的大小，用于确定虚拟节点的编号
    busDC_size = size(jpc.busDC, 1) - num_batteries
    
    for (i, battery) in enumerate(batteries)
        # 获取电池连接的实际节点编号
        actual_bus = battery.bus
        
        # 计算电池虚拟节点编号（基于之前在JPC_battery_process中的编号规则）
        virtual_bus = busDC_size + i
        
        # 获取节点的基准电压(kV)
        basekv = 0.0
        for j in 1:size(jpc.busDC, 1)
            if jpc.busDC[j, 1] == actual_bus
                basekv = jpc.busDC[j, BASE_KV]
                break
            end
        end
        
        # 计算基准阻抗
        baseR = (basekv^2) / case.baseMVA
        
        # 计算标幺值阻抗（使用电池的内阻）
        # r_pu = battery.ra / baseR
        # r_pu = 0.0242/baseR  # 假设电池内阻为0.0242Ω
        r_pu = 0.0252115/baseR  # 假设电池内阻为0.0249Ω
        x_pu = 0  # 直流系统无感抗，设置为一个非常小的值
        
        # 填充虚拟支路矩阵
        battery_branches[i, F_BUS] = virtual_bus       # 虚拟节点
        battery_branches[i, T_BUS] = actual_bus        # 实际连接节点
        battery_branches[i, BR_R] = r_pu               # 标幺值电阻
        battery_branches[i, BR_X] = x_pu               # 标幺值电抗
        
        # 设置额定容量（基于电池参数计算）
        # 假设电池额定容量可以从电池参数计算得到
        rated_capacity = battery.package * battery.voc  # 简化计算，实际可能需要更复杂的公式
        battery_branches[i, RATE_A] = rated_capacity
        
        # 设置支路状态
        battery_branches[i, BR_STATUS] = battery.in_service ? 1.0 : 0.0
        
        # 设置相角限制（直流系统中通常不受限制）
        battery_branches[i, ANGMIN] = -360.0
        battery_branches[i, ANGMAX] = 360.0
    end
    
    # 将电池虚拟支路添加到branchDC中
    if isempty(jpc.branchDC)
        jpc.branchDC = battery_branches
    else
        jpc.branchDC = [jpc.branchDC; battery_branches]
    end
    
    return jpc
end

function JPC_battery_soc_process(case::JuliaPowerCase, jpc::JPC)
    # 处理电池SOC数据，转换为JPC格式
    batteries = deepcopy(case.storages)
    num_batteries = length(batteries)
    
    # 如果没有电池，直接返回
    if num_batteries == 0
        return jpc
    end
    
    # 创建电池SOC矩阵
    battery_soc = zeros(num_batteries, 8)  
    
    for (i, battery) in enumerate(batteries)
        battery_soc[i, 1] = battery.bus  # 电池连接的母线ID
        battery_soc[i, 2] = battery.power_capacity_mw   # 电池的SOC值（标幺值）
        battery_soc[i, 3] = battery.energy_capacity_mwh  # 电池的有功功率（MW）
        battery_soc[i, 4] = battery.soc_init  # 电池的无功功率（MVAR）
        battery_soc[i, 5] = battery.min_soc  # 电池的最大有功功率（MW）
        battery_soc[i, 6] = battery.max_soc  # 电池的最小有功功率（MW）
        battery_soc[i, 7] = battery.efficiency  # 电池的最大无功功率（MVAR）
        battery_soc[i, 8] = battery.in_service ? 1.0 : 0.0  # 电池是否在服务中（1.0表示在服务，0.0表示不在服务）
    end
    for (i, battery) in enumerate(batteries)
        # 获取电池连接的实际节点编号
        bus_id = battery.bus
        # 在JPC的busDC中查找对应的节点
        bus_index = findfirst(x -> x[1] == bus_id, jpc.busDC[:, 1])
        jpc.busDC[bus_index, PD] -= 0.0
        loadDC = zeros(1, 8)  # 创建一个空的负荷矩阵
        nd = size(jpc.busDC, 1)
        loadDC[1, 1] = nd + 1  # 设置负荷对应的母线ID
        loadDC[1, 2] = bus_index
        loadDC[1, 3] = 1 # inservice
        loadDC[1, 4] = 0.0
        loadDC[1, 5] = 0.0
        loadDC[1, 6] = 0.0  
        loadDC[1, 7] = 0.0
        loadDC[1, 8] = 1.0
        # 将负荷数据添加到JPC的负荷矩阵中
        if isempty(jpc.loadDC)
            jpc.loadDC = loadDC
        else
            jpc.loadDC = [jpc.loadDC; loadDC]
        end
    end
    # 将电池SOC数据添加到JPC结构体
    jpc.storage = battery_soc
    
    return jpc
end


function calculate_line_parameters(case::JuliaPowerCase, jpc::JPC)
    # 处理线路数据，转换为JPC格式
    nbr = length(case.branchesAC)
    branch = zeros(nbr, 14)
    lines = case.branchesAC

    for (i, line) in enumerate(lines)
        # 获取起始和终止母线编号
        from_bus_idx = line.from_bus
        to_bus_idx = line.to_bus
        
        # 获取起始母线的基准电压(kV)
        basekv = jpc.busAC[from_bus_idx, BASE_KV]
        
        # 计算基准阻抗
        baseR = (basekv^2) / case.baseMVA
        
        # 考虑并联线路的情况
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # 计算标幺值阻抗
        r_pu = line.length_km * line.r_ohm_per_km / baseR / parallel
        x_pu = line.length_km * line.x_ohm_per_km / baseR / parallel
        
        # 计算并联电纳(p.u.)
        b_pu = 2 * π * case.basef * line.length_km * line.c_nf_per_km * 1e-9 * baseR * parallel
        
        # 计算并联电导(p.u.)
        g_pu = 0.0
        if hasfield(typeof(line), :g_us_per_km)
            g_pu = line.g_us_per_km * 1e-6 * baseR * line.length_km * parallel
        end
        
        # 填充branchAC矩阵
        branch[i, F_BUS] = from_bus_idx
        branch[i, T_BUS] = to_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        branch[i, BR_B] = b_pu
        
        # 设置额定容量
        if hasfield(typeof(line), :max_i_ka)
            branch[i, RATE_A] = line.max_i_ka * basekv * sqrt(3)  # 额定容量(MVA)
        else
            branch[i, RATE_A] = 100.0  # 默认值
        end
        
        # 设置支路状态
        branch[i, BR_STATUS] = line.in_service ? 1.0 : 0.0
        
        # 设置相角限制
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end

    jpc.branchAC = branch
end

function calculate_transformer2w_parameters(case::JuliaPowerCase, jpc::JPC)
    # 处理变压器数据，转换为JPC格式
    transformers = case.transformers_2w_etap
    nbr = length(transformers)
    
    if nbr == 0
        return  # 如果没有变压器，直接返回
    end
    
    # 创建变压器分支矩阵
    branch = zeros(nbr, 14)
    
    for (i, transformer) in enumerate(transformers)
        # 获取高压侧和低压侧母线编号
        hv_bus_idx = transformer.hv_bus
        lv_bus_idx = transformer.lv_bus
        
        # 获取高压侧母线的基准电压(kV)
        hv_basekv = jpc.busAC[hv_bus_idx, BASE_KV]
        
        # 计算阻抗参数
        # 变压器阻抗百分比转换为标幺值
        z_pu = transformer.z_percent
        x_r_ratio = transformer.x_r
        
        # 计算电阻和电抗（考虑基准功率转换）
        s_ratio = transformer.sn_mva / case.baseMVA
        z_pu = z_pu / s_ratio  # 转换到系统基准
        
        r_pu = z_pu / sqrt(1 + x_r_ratio^2)
        x_pu = r_pu * x_r_ratio
        
        # 考虑并联变压器
        parallel = transformer.parallel
        if parallel > 1
            r_pu = r_pu / parallel
            x_pu = x_pu / parallel
        end
        
        # 填充分支矩阵
        branch[i, F_BUS] = hv_bus_idx
        branch[i, T_BUS] = lv_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        branch[i, BR_B] = 0.0  # 变压器通常没有并联电纳
        
        # 设置变比和相移
        branch[i, TAP] = 1.0  # 默认变比为1.0
        branch[i, SHIFT] = 0.0  # 默认相移角度为0.0
        
        # 设置额定容量
        branch[i, RATE_A] = case.baseMVA 
        
        # 设置支路状态
        branch[i, BR_STATUS] = transformer.in_service ? 1.0 : 0.0
        
        # 设置相角限制
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end
    
    # 将变压器分支数据添加到JPC结构体
    if isempty(jpc.branchAC)
        jpc.branchAC = branch
    else
        jpc.branchAC = [jpc.branchAC; branch]
    end
end

function JPC_gens_process(case::JuliaPowerCase, jpc::JPC)
    # 统计各类发电设备数量
    n_gen = length(case.gensAC)
    n_sgen = length(case.sgensAC)
    n_ext = length(case.ext_grids)
    
    # 计算总发电设备数量
    total_gens = n_gen + n_sgen + n_ext
    
    if total_gens == 0
        return  # 如果没有发电设备，直接返回
    end
    
    # 创建发电机矩阵，行数为发电设备数量，列数为26
    gen_data = zeros(total_gens, 26)
    
    # 处理外部电网(通常作为平衡节点/参考节点)
    for (i, ext) in enumerate(case.ext_grids)
        if !ext.in_service
            continue
        end
        
        bus_idx = ext.bus
        
        # 填充发电机数据
        gen_data[i, :] = [
            bus_idx,        # 发电机连接的母线编号
            0.0,            # 有功功率输出(MW)
            0.0,            # 无功功率输出(MVAr)
            9999.0,         # 最大无功功率输出(MVAr)
            -9999.0,        # 最小无功功率输出(MVAr)
            ext.vm_pu,      # 电压幅值设定值(p.u.)
            case.baseMVA,   # 发电机基准容量(MVA)
            1.0,            # 发电机状态(1=运行, 0=停运)
            9999.0,         # 最大有功功率输出(MW)
            -9999.0,        # 最小有功功率输出(MW)
            0.0,            # PQ能力曲线低端有功功率输出(MW)
            0.0,            # PQ能力曲线高端有功功率输出(MW)
            0.0,            # PC1处最小无功功率输出(MVAr)
            0.0,            # PC1处最大无功功率输出(MVAr)
            0.0,            # PC2处最小无功功率输出(MVAr)
            0.0,            # PC2处最大无功功率输出(MVAr)
            0.0,            # AGC调节速率(MW/min)
            0.0,            # 10分钟备用调节速率(MW)
            0.0,            # 30分钟备用调节速率(MW)
            0.0,            # 无功功率调节速率(MVAr/min)
            1.0,            # 区域参与因子
            2.0,            # 发电机模型(2=多项式成本模型)
            0.0,            # 启动成本(美元)
            0.0,            # 关机成本(美元)
            3.0,            # 多项式成本函数系数数量
            0.0             # 成本函数参数(后续需要扩展)
        ]
        
        # 更新母线类型为参考节点(REF/平衡节点)
        jpc.busAC[bus_idx, 2] = 3  # 3表示REF节点
    end
    
    # 处理常规发电机(通常作为PV节点)
    offset = n_ext
    for (i, gen) in enumerate(case.gensAC)
        if !gen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = gen.bus
        
        # 计算无功功率(如果没有直接给出)
        q_mvar = 0.0
        if hasfield(typeof(gen), :q_mvar)
            q_mvar = gen.q_mvar
        else
            # 根据功率因数计算无功功率
            p_mw = gen.p_mw * gen.scaling
            if gen.cos_phi > 0 && p_mw > 0
                q_mvar = p_mw * tan(acos(gen.cos_phi))
            end
        end
        
        # 基准容量
        mbase = gen.sn_mva > 0 ? gen.sn_mva : case.baseMVA
        
        # 爬坡率参数
        ramp_agc = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                   gen.ramp_up_rate_mw_per_min : 
                   (gen.max_p_mw - gen.min_p_mw) / 10
        ramp_10 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 10 : 
                  gen.max_p_mw - gen.min_p_mw
        ramp_30 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 30 : 
                  gen.max_p_mw - gen.min_p_mw
        
        # 填充发电机数据
        gen_data[idx, :] = [
            bus_idx,                               # 发电机连接的母线编号
            gen.p_mw * gen.scaling,                # 有功功率输出(MW)
            q_mvar,                                # 无功功率输出(MVAr)
            gen.max_q_mvar,                        # 最大无功功率输出(MVAr)
            gen.min_q_mvar,                        # 最小无功功率输出(MVAr)
            gen.vm_pu,                             # 电压幅值设定值(p.u.)
            mbase,                                 # 发电机基准容量(MVA)
            1.0,                                   # 发电机状态(1=运行, 0=停运)
            gen.max_p_mw,                          # 最大有功功率输出(MW)
            gen.min_p_mw,                          # 最小有功功率输出(MW)
            gen.min_p_mw,                          # PQ能力曲线低端有功功率输出(MW)
            gen.max_p_mw,                          # PQ能力曲线高端有功功率输出(MW)
            gen.min_q_mvar,                        # PC1处最小无功功率输出(MVAr)
            gen.max_q_mvar,                        # PC1处最大无功功率输出(MVAr)
            gen.min_q_mvar,                        # PC2处最小无功功率输出(MVAr)
            gen.max_q_mvar,                        # PC2处最大无功功率输出(MVAr)
            ramp_agc,                              # AGC调节速率(MW/min)
            ramp_10,                               # 10分钟备用调节速率(MW)
            ramp_30,                               # 30分钟备用调节速率(MW)
            (gen.max_q_mvar - gen.min_q_mvar)/10,  # 无功功率调节速率(MVAr/min)
            1.0,                                   # 区域参与因子
            2.0,                                   # 发电机模型(2=多项式成本模型)
            0.0,                                   # 启动成本(美元)
            0.0,                                   # 关机成本(美元)
            3.0,                                   # 多项式成本函数系数数量
            0.0                                    # 成本函数参数(后续需要扩展)
        ]
        
        # 如果母线尚未设置为参考节点，则设置为PV节点
        if jpc.busAC[bus_idx, 2] != 3  # 3表示REF节点
            jpc.busAC[bus_idx, 2] = 2  # 2表示PV节点
        end
    end
    
    # 处理静态发电机(通常作为PQ节点，但如果有电压控制能力，也可以是PV节点)
    offset = n_ext + n_gen
    for (i, sgen) in enumerate(case.sgensAC)
        if !sgen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = sgen.bus
        
        # 填充发电机数据
        gen_data[idx, :] = [
            bus_idx,                                # 发电机连接的母线编号
            sgen.p_mw * sgen.scaling,               # 有功功率输出(MW)
            sgen.q_mvar * sgen.scaling,             # 无功功率输出(MVAr)
            sgen.max_q_mvar,                        # 最大无功功率输出(MVAr)
            sgen.min_q_mvar,                        # 最小无功功率输出(MVAr)
            1.0,                                    # 电压幅值设定值(p.u.)
            case.baseMVA,                           # 发电机基准容量(MVA)
            1.0,                                    # 发电机状态(1=运行, 0=停运)
            sgen.max_p_mw,                          # 最大有功功率输出(MW)
            sgen.min_p_mw,                          # 最小有功功率输出(MW)
            sgen.min_p_mw,                          # PQ能力曲线低端有功功率输出(MW)
            sgen.max_p_mw,                          # PQ能力曲线高端有功功率输出(MW)
            sgen.min_q_mvar,                        # PC1处最小无功功率输出(MVAr)
            sgen.max_q_mvar,                        # PC1处最大无功功率输出(MVAr)
            sgen.min_q_mvar,                        # PC2处最小无功功率输出(MVAr)
            sgen.max_q_mvar,                        # PC2处最大无功功率输出(MVAr)
            (sgen.max_p_mw - sgen.min_p_mw) / 10,   # AGC调节速率(MW/min)
            sgen.max_p_mw - sgen.min_p_mw,          # 10分钟备用调节速率(MW)
            sgen.max_p_mw - sgen.min_p_mw,          # 30分钟备用调节速率(MW)
            (sgen.max_q_mvar - sgen.min_q_mvar)/10, # 无功功率调节速率(MVAr/min)
            1.0,                                    # 区域参与因子
            2.0,                                    # 发电机模型(2=多项式成本模型)
            0.0,                                    # 启动成本(美元)
            0.0,                                    # 关机成本(美元)
            3.0,                                    # 多项式成本函数系数数量
            0.0                                     # 成本函数参数(后续需要扩展)
        ]
        
        # 如果静态发电机可控且母线尚未设置为REF或PV节点，则可能设置为PV节点
        if sgen.controllable && jpc.busAC[bus_idx, 2] == 1  # 1表示PQ节点
            jpc.busAC[bus_idx, 2] = 2  # 2表示PV节点
        end
    end
    
    # 移除未使用的行(对应未投运的发电设备)
    active_rows = findall(x -> x > 0, gen_data[:, 8])  # 第8列是GEN_STATUS
    gen_data = gen_data[active_rows, :]
    
    # 将发电机数据存储到JPC结构体
    jpc.genAC = gen_data
    
    # 确保至少有一个平衡节点
    if !any(jpc.busAC[:, 2] .== 3) && size(gen_data, 1) > 0  # 3表示REF节点
        # 如果没有平衡节点，选择第一个发电机所在母线作为平衡节点
        first_gen_bus = Int(gen_data[1, 1])
        jpc.busAC[first_gen_bus, 2] = 3  # 3表示REF节点
    end
end

function JPC_battery_gens_process(case::JuliaPowerCase, jpc::JPC)
    # 为电池虚拟节点创建虚拟发电机
    batteries = deepcopy(case.storageetap)
    num_batteries = length(batteries)
    
    # 如果没有电池，直接返回
    if num_batteries == 0
        return jpc
    end
    
    # 获取当前busDC的大小，用于确定虚拟节点的编号
    busDC_size = size(jpc.busDC, 1) - num_batteries
    
    # 创建电池虚拟发电机矩阵
    # genDC矩阵通常包含以下列：
    # [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...]
    # 具体列数应与你的genDC结构一致
    num_gen_cols = size(jpc.genAC, 2)
    if num_gen_cols == 0  # 如果genDC为空，设置一个默认列数
        num_gen_cols = 10
    end
    
    battery_gens = zeros(num_batteries, num_gen_cols)
     # 创建storage矩阵
    num_storage_cols = 5  # 根据idx_ess函数定义的列数
    storage_matrix = zeros(num_batteries, num_storage_cols)
    
    for (i, battery) in enumerate(batteries)
        # 计算电池虚拟节点编号
        virtual_bus = busDC_size + i
        
        # 计算电池的功率容量（基于电池参数）
        # 这里使用简化计算，实际应根据电池特性进行更准确的计算
        power_capacity = battery.package * battery.voc
        
        # 填充虚拟发电机矩阵
        battery_gens[i, 1] = virtual_bus       # GEN_BUS: 发电机连接的节点编号
        battery_gens[i, 2] = 0.0               # PG: 初始有功功率输出(MW)，初始设为0
        battery_gens[i, 3] = 0.0               # QG: 初始无功功率输出(MVAR)，直流系统通常为0
        
        # 设置无功功率限制（直流系统通常不考虑）
        battery_gens[i, 4] = 0.0               # QMAX: 最大无功功率输出
        battery_gens[i, 5] = 0.0               # QMIN: 最小无功功率输出
        
        # 设置电压和基准功率
        battery_gens[i, 6] = 1.0               # VG: 电压设定值(p.u.)
        battery_gens[i, 7] = case.baseMVA      # MBASE: 发电机基准功率(MVA)
        
        # 设置发电机状态
        battery_gens[i, 8] = battery.in_service ? 1.0 : 0.0  # GEN_STATUS: 发电机状态
        
        # 设置有功功率限制（充电为负，放电为正）
        battery_gens[i, 9] = power_capacity    # PMAX: 最大有功功率输出(MW)，放电功率
        battery_gens[i, 10] = -power_capacity  # PMIN: 最小有功功率输出(MW)，充电功率

         # 填充storage矩阵
        storage_matrix[i, ESS_BUS] = virtual_bus               # ESS_BUS: 连接的节点编号
        # storage_matrix[i, ESS_POWER_CAPACITY] = power_capacity # ESS_POWER_CAPACITY: 功率容量(MW)
        # storage_matrix[i, ESS_ENERGY_CAPACITY] = 0
        # storage_matrix[i, ESS_AREA] = 1                        # ESS_AREA: 区域编号，默认为1
        
        # 如果genDC有更多列，根据需要设置其他参数
        if num_gen_cols > 10
            # 例如，设置爬坡率限制、成本系数等
            # 这里需要根据你的系统具体需求进行设置
            for j in 11:num_gen_cols
                battery_gens[i, j] = 0.0  # 默认设为0
            end
        end
    end
    
    # 将电池虚拟发电机添加到genDC中
    if isempty(jpc.genDC)
        jpc.genDC = battery_gens
    else
        jpc.genDC = [jpc.genDC; battery_gens]
    end
    
    # 将存储设备信息添加到storage中
    if !isdefined(jpc, :storageetap) || isempty(jpc.storageetap)
        jpc.storageetap = storage_matrix
    else
        jpc.storageetap = [jpc.storageetap; storage_matrix]
    end

    return jpc
end


function JPC_loads_process(case::JuliaPowerCase, jpc::JPC)
    # 处理负荷数据，转换为JPC格式并更新busAC的PD和QD
    
    # 过滤出投运的负荷
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # 如果没有投运的负荷，直接返回
    if isempty(in_service_loads)
        return
    end
    
    # 创建一个空矩阵，行数为负荷数，列数为8
    num_loads = length(in_service_loads)
    load_matrix = zeros(num_loads, 8)
    
    # 创建一个字典，用于累加连接到同一母线的负荷
    bus_load_sum = Dict{Int, Vector{Float64}}()
    
    for (i, load) in enumerate(in_service_loads)
        # 计算实际的有功和无功负荷（考虑缩放因子）
        # if mode =="1_ph_pf"
            actual_p_mw = load.p_mw * load.scaling
            actual_q_mvar = load.q_mvar * load.scaling
        # else
        #     actual_p_mw = load.p_mw * load.scaling / 3.0
        #     actual_q_mvar = load.q_mvar * load.scaling / 3.0
        # end
        
        # 填充负荷矩阵的每一行
        load_matrix[i, :] = [
            i,              # 负荷连接的母线编号
            load.bus,                     # 负荷编号
            1.0,                   # 负荷状态(1=投运)
            actual_p_mw,           # 有功负荷(MW)
            actual_q_mvar,         # 无功负荷(MVAr)
            load.const_z_percent/100,  # 恒阻抗负荷百分比
            load.const_i_percent/100,  # 恒电流负荷百分比
            load.const_p_percent/100   # 恒功率负荷百分比
        ]
        
        # 累加连接到同一母线的负荷
        bus_idx = load.bus
        if haskey(bus_load_sum, bus_idx)
            bus_load_sum[bus_idx][1] += actual_p_mw
            bus_load_sum[bus_idx][2] += actual_q_mvar
        else
            bus_load_sum[bus_idx] = [actual_p_mw, actual_q_mvar]
        end
    end
    
    # 将负荷数据存储到JPC结构体
    jpc.loadAC = load_matrix
    
    # 更新busAC矩阵中的PD和QD字段
    for (bus_idx, load_values) in bus_load_sum
        # 找到对应的母线行
        bus_row = findfirst(x -> x == bus_idx, jpc.busAC[:, 1])
        
        if !isnothing(bus_row)
            # 更新PD(第3列)和QD(第4列)
            jpc.busAC[bus_row, PD] = load_values[1]  # PD - 有功负荷(MW)
            jpc.busAC[bus_row, QD] = load_values[2]  # QD - 无功负荷(MVAr)
        end
    end
end

function JPC_dcloads_process(case::JuliaPowerCase, jpc::JPC)
    # 处理直流负荷数据，转换为JPC格式并更新busDC的PD和QD
    
    # 过滤出投运的直流负荷
    in_service_dcloads = filter(dcload -> dcload.in_service == true, case.loadsDC)
    
    # 如果没有投运的直流负荷，直接返回
    if isempty(in_service_dcloads)
        return
    end
    
    # 创建一个空矩阵，行数为直流负荷数，列数为8
    num_dcloads = length(in_service_dcloads)
    dcload_matrix = zeros(num_dcloads, 8)
    
    for (i, dcload) in enumerate(in_service_dcloads)
        # 填充直流负荷矩阵的每一行
        dcload_matrix[i, :] = [
            dcload.index,              # 直流负荷编号
            dcload.bus,     # 直流负荷连接的母线编号
            1.0,            # 直流负荷状态(1=投运)
            dcload.p_mw * dcload.scaling,  # 有功负荷(MW)
            0.0,            # 无功负荷(MVAr)
            dcload.const_z_percent/100,           # 恒阻抗百分比(默认0)
            dcload.const_i_percent/100,           # 恒电流百分比(默认0)
            dcload.const_p_percent/100            # 恒功率百分比(默认0)
        ]
        
        # 更新busDC矩阵中的PD和QD字段
        bus_row = findfirst(x -> x == dcload.bus, jpc.busDC[:, 1])
        
        if !isnothing(bus_row)
            jpc.busDC[bus_row, PD] += dcload_matrix[i, 4]  # PD - 有功负荷(MW)
            jpc.busDC[bus_row, QD] += dcload_matrix[i, 5]  # QD - 无功负荷(MVAr)
        end
    end
    
    # 将直流负荷数据存储到JPC结构体
    jpc.loadDC = dcload_matrix

    return jpc
end

function JPC_pv_process(case::JuliaPowerCase, jpc::JPC)
    # 处理光伏发电机数据，转换为JPC格式并更新busAC的PD和QD
    
    # 过滤出投运的光伏发电机
    in_service_pvs = filter(pv -> pv.in_service == true, case.pvarray)
    
    # 如果没有投运的光伏发电机，直接返回
    if isempty(in_service_pvs)
        return
    end
    
    # 创建一个空矩阵，行数为光伏发电机数，列数为8
    num_pvs = length(in_service_pvs)
    pv_matrix = zeros(num_pvs, 9)
    
    for (i, pv) in enumerate(in_service_pvs)
        Voc = (pv.voc + (pv.temperature - 25)*pv.β_voc )* pv.numpanelseries
        Vmpp = pv.vmpp * pv.numpanelseries
        Isc = (pv.isc + (pv.temperature - 25)*pv.α_isc)*(pv.irradiance/1000.0) * pv.numpanelparallel
        Impp = pv.impp * (pv.irradiance/1000.0) * pv.numpanelparallel

        # Voc = (pv.voc + (pv.temperature - 25)*pv.β_voc + pv.γ_voc*log(pv.irradiance/1000.0)) * pv.numpanelseries
        # Vmpp = (pv.vmpp + pv.γ_vmpp*log(pv.irradiance/1000.0)) * pv.numpanelseries
        # Isc = (pv.isc + (pv.temperature - 25)*pv.α_isc)*(pv.irradiance/1000.0) * pv.numpanelparallel
        # Impp = pv.impp*(pv.irradiance/1000.0) * pv.numpanelparallel
        # 填充光伏发电机矩阵的每一行
        pv_matrix[i, :] = [
            i,              # 光伏发电机编号
            pv.bus,         # 光伏发电机连接的母线编号
            Voc,         # 光伏发电机额定电压(V)
            Vmpp,        # 光伏发电机额定电压(V)
            Isc,         # 光伏发电机短路电流(A)
            Impp,        # 光伏发电机额定电流(A)
            pv.irradiance,  # 光伏发电机辐照度(W/m²)
            1.0,            # area
            1.0,            # 光伏发电机状态(1=投运)
        ]
        
        # # 更新busAC矩阵中的PD和QD字段
        # bus_row = findfirst(x -> x == pv.bus, jpc.busAC[:, 1])
        
        # if !isnothing(bus_row)
        #     jpc.busAC[bus_row, PD] += pv_matrix[i, 4]  # PD - 有功负荷(MW)
        #     jpc.busAC[bus_row, QD] += pv_matrix[i, 5]  # QD - 无功负荷(MVAr)
        # end
    end
    
    # 将光伏发电机数据存储到JPC结构体
    jpc.pv = pv_matrix

    return jpc
    
end

function JPC_ac_pv_system_process(case::JuliaPowerCase, jpc::JPC)
    # 处理交流侧光伏系统数据，根据控制模式转换为发电机或负荷
    
    # 过滤出投运的交流侧光伏系统
    in_service_ac_pvs = filter(ac_pv -> ac_pv.in_service == true, case.ACPVSystems)
    
    # 如果没有投运的交流侧光伏系统，直接返回
    if isempty(in_service_ac_pvs)
        return jpc
    end
    
    # 创建一个空矩阵，行数为交流侧光伏系统数，列数为13
    num_ac_pvs = length(in_service_ac_pvs)
    ac_pv_matrix = zeros(num_ac_pvs, 15)
    for (i, ac_pv) in enumerate(in_service_ac_pvs)
        Vmpp = ac_pv.vmpp * ac_pv.numpanelseries
        Voc = (ac_pv.voc + (ac_pv.temperature - 25) * ac_pv.β_voc) * ac_pv.numpanelseries
        Isc = (ac_pv.isc + (ac_pv.temperature - 25) * ac_pv.α_isc) * (ac_pv.irradiance / 1000.0) * ac_pv.numpanelparallel
        Impp = ac_pv.impp * ac_pv.numpanelparallel

        p_max = Vmpp * Impp / 1000000.0 * (1-ac_pv.loss_percent)# 最大有功输出(MW)

        if ac_pv.control_mode == "Voltage Control"
            mode = 1
        else
            mode = 0 
        end
        
        # 填充交流侧光伏系统矩阵的每一行
        ac_pv_matrix[i, :] = [
            i,                # 交流光伏系统编号
            ac_pv.bus,            # 连接母线编号
            Voc,               # 光伏系统额定电压(V)
            Vmpp,              # 光伏系统额定电压(V)
            Isc,               # 光伏系统短路电流(A)
            Impp,              # 光伏系统额定电流(A)
            ac_pv.irradiance,  # 光伏系统辐照度(W/m²)
            ac_pv.loss_percent,
            mode,              # 控制模式(0=无功控制, 1=电压控制)
            ac_pv.p_mw,           # 有功出力(MW)
            ac_pv.q_mvar,         # 无功出力(MVAr)
            ac_pv.max_q_mvar,     # 无功上限(MVAr)
            ac_pv.min_q_mvar,     # 无功下限(MVAr)
            1,            # 区域编号
            ac_pv.in_service ? 1.0 : 0.0  # 光伏系统状态(1=投运, 0=停运)
        ]
    end
    # 将交流侧光伏系统数据存储到JPC结构体
    jpc.pv_acsystem = ac_pv_matrix


    # # 分离不同控制模式的光伏系统
    # voltage_control_pvs = filter(ac_pv -> ac_pv.control_mode == "Voltage Control", in_service_ac_pvs)
    # mvar_control_pvs = filter(ac_pv -> ac_pv.control_mode == "Mvar Control", in_service_ac_pvs)
    
    # # 处理电压控制模式的光伏系统（创建发电机，修改母线类型为PV）
    # if !isempty(voltage_control_pvs)
    #     num_voltage_pvs = length(voltage_control_pvs)
    #     voltage_pv_matrix = zeros(num_voltage_pvs, 26)  # 发电机矩阵列数
        
    #     for (i, ac_pv) in enumerate(voltage_control_pvs)
    #         # 获取当前最大发电机编号
    #         max_gen_id = isempty(jpc.genAC) ? 0 : maximum(jpc.genAC[:, 1])
    #         gen_id = max_gen_id + i

    #         Vmpp = ac_pv.vmpp * ac_pv.numpanelseries
    #         Voc = (ac_pv.voc + (ac_pv.temperature - 25) * ac_pv.β_voc) * ac_pv.numpanelseries
    #         Isc = (ac_pv.isc + (ac_pv.temperature - 25) * ac_pv.α_isc) * (ac_pv.irradiance / 1000.0) * ac_pv.numpanelparallel
    #         Impp = ac_pv.impp * ac_pv.numpanelparallel

    #         p_max = Vmpp * Impp / 1000000.0 * (1-ac_pv.loss_percent)# 最大有功输出(MW)
    #         findfirst_bus = findfirst(x -> x == ac_pv.bus, jpc.busAC[:, 1])
    #         jpc.busAC[findfirst_bus, BUS_TYPE] = 2 # 设置母线类型为PV节点
            
    #         voltage_pv_matrix[i, :] = [
    #             ac_pv.bus,                # 连接母线编号
    #             ac_pv.p_mw,              # 有功出力(MW)
    #             ac_pv.q_mvar,            # 无功出力(MVAr)
    #             ac_pv.max_q_mvar,        # 无功上限(MVAr)
    #             ac_pv.min_q_mvar,        # 无功下限(MVAr)
    #             hasfield(typeof(ac_pv), :vm_ac_pu) ? ac_pv.vm_ac_pu : 1.0, # 电压设定值(p.u.)
    #             case.baseMVA,            # 发电机基准容量(MVA)
    #             1.0,                     # 发电机状态(1=投运)
    #             p_max,                   # 有功上限(MW)
    #             0.0,                     # 有功下限(MW)
    #             0.0,                     # PQ能力曲线低端有功输出(MW)
    #             0.0,                     # PQ能力曲线高端有功输出(MW)
    #             0.0,                     # PC1处最小无功输出(MVAr)
    #             0.0,                     # PC1处最大无功输出(MVAr)
    #             0.0,                     # PC2处最小无功输出(MVAr)
    #             0.0,                     # PC2处最大无功输出(MVAr)
    #             0.0,                     # AGC调节速率(MW/min)
    #             0.0,                     # 10分钟备用调节速率(MW)
    #             0.0,                     # 30分钟备用调节速率(MW)
    #             0.0,                     # 无功功率调节速率(MVAr/min)
    #             0.0,                     # 区域参与因子
    #             0.0,                     # 发电机模型(0=无模型)
    #             0.0,                     # 启动成本(美元)
    #             0.0,                     # 关机成本(美元)
    #             0.0,                     # 多项式成本函数系数数量
    #             0.0                      # 成本函数参数(后续需要扩展)
    #         ]
    #     end
        
    #     # 将电压控制光伏系统添加到发电机矩阵
    #     if isempty(jpc.genAC)
    #         jpc.genAC = voltage_pv_matrix
    #     else
    #         jpc.genAC = vcat(jpc.genAC, voltage_pv_matrix)
    #     end
    # end
    
    # # 处理无功控制模式的光伏系统（创建发电机，但不修改母线类型）
    # if !isempty(mvar_control_pvs)
    #     num_mvar_pvs = length(mvar_control_pvs)
    #     mvar_pv_matrix = zeros(num_mvar_pvs, 26)  # 发电机矩阵列数
        
    #     for (i, ac_pv) in enumerate(mvar_control_pvs)
    #         # 获取当前最大发电机编号
    #         current_max_gen_id = isempty(jpc.genAC) ? 0 : maximum(jpc.genAC[:, 1])
    #         # 如果已经有voltage control的发电机，需要在其基础上继续编号
    #         gen_id = current_max_gen_id + length(voltage_control_pvs) + i

    #         Vmpp = ac_pv.vmpp * ac_pv.numpanelseries
    #         Voc = (ac_pv.voc + (ac_pv.temperature - 25) * ac_pv.β_voc) * ac_pv.numpanelseries
    #         Isc = (ac_pv.isc + (ac_pv.temperature - 25) * ac_pv.α_isc) * (ac_pv.irradiance / 1000.0) * ac_pv.numpanelparallel
    #         Impp = ac_pv.impp * ac_pv.numpanelparallel

    #         p_max = Vmpp * Impp / 1000000.0 * (1-ac_pv.loss_percent)# 最大有功输出(MW)
            
    #         # 注意：MVar Control模式不修改母线类型，保持原有类型（通常为PQ节点）
            
    #         mvar_pv_matrix[i, :] = [
    #             ac_pv.bus,                # 连接母线编号
    #             ac_pv.p_mw,              # 有功出力(MW) - MPPT功率
    #             ac_pv.q_mvar,            # 无功出力(MVAr) - 固定无功
    #             ac_pv.max_q_mvar,        # 无功上限(MVAr)
    #             ac_pv.min_q_mvar,        # 无功下限(MVAr)
    #             1.0,                     # 电压设定值(p.u.) - MVar控制不控制电压
    #             case.baseMVA,            # 发电机基准容量(MVA)
    #             1.0,                     # 发电机状态(1=投运)
    #             p_max,                   # 有功上限(MW)
    #             0.0,                     # 有功下限(MW)
    #             0.0,                     # PQ能力曲线低端有功输出(MW)
    #             0.0,                     # PQ能力曲线高端有功输出(MW)
    #             0.0,                     # PC1处最小无功输出(MVAr)
    #             0.0,                     # PC1处最大无功输出(MVAr)
    #             0.0,                     # PC2处最小无功输出(MVAr)
    #             0.0,                     # PC2处最大无功输出(MVAr)
    #             0.0,                     # AGC调节速率(MW/min)
    #             0.0,                     # 10分钟备用调节速率(MW)
    #             0.0,                     # 30分钟备用调节速率(MW)
    #             0.0,                     # 无功功率调节速率(MVAr/min)
    #             0.0,                     # 区域参与因子
    #             0.0,                     # 发电机模型(0=无模型)
    #             0.0,                     # 启动成本(美元)
    #             0.0,                     # 关机成本(美元)
    #             0.0,                     # 多项式成本函数系数数量
    #             0.0                      # 成本函数参数(后续需要扩展)
    #         ]
    #     end
        
    #     # 将无功控制光伏系统添加到发电机矩阵
    #     if isempty(jpc.genAC)
    #         jpc.genAC = mvar_pv_matrix
    #     else
    #         jpc.genAC = vcat(jpc.genAC, mvar_pv_matrix)
    #     end
    # end
    
    return jpc
end



function JPC_inverters_process(case::JuliaPowerCase, jpc::JPC)
    # 处理逆变器数据，转换为JPC格式并更新busAC和busDC的负荷
    
    # 过滤出投运的逆变器
    in_service_inverters = filter(inverter -> inverter.in_service == true, case.converters)
    
    # 如果没有投运的逆变器，直接返回
    if isempty(in_service_inverters)
        return jpc
    end
    
    # 获取当前负荷数量，用于新增负荷的编号
    nld_ac = size(jpc.loadAC, 1)  # 交流侧负荷数量
    nld_dc = size(jpc.loadDC, 1)  # 直流侧负荷数量
    
    # 创建用于存储需要新增的负荷记录
    # 使用矩阵而不是数组来存储新负荷
    num_cols_ac = size(jpc.loadAC, 2)
    num_cols_dc = size(jpc.loadDC, 2)
    
    # 计算需要添加的最大可能负荷数（每个逆变器最多添加一个负荷）
    max_new_loads = length(in_service_inverters)
    new_loads_ac = zeros(0, num_cols_ac)  # 创建一个空矩阵，行数为0，列数与loadAC相同
    new_loads_dc = zeros(0, num_cols_dc)  # 创建一个空矩阵，行数为0，列数与loadDC相同

    #创建converter空矩阵
    converters = zeros(0, 18)
    
    # 跟踪新增负荷的数量
    new_ac_load_count = 0
    new_dc_load_count = 0

    
    
   for (i, inverter) in enumerate(in_service_inverters)
        #为converter矩阵添加连接关系
        converter = zeros(1, 18)  # 创建一行
        
        # 逆变器的工作模式
        mode = inverter.control_mode
        if mode == "δs_Us"
            converter[1, CONV_MODE] = 1.0  # δs_Us模式的逆变器不需要设置
        elseif mode == "Ps_Qs"
            converter[1, CONV_MODE] = 0.0  # Ps_Qs模式的逆变器
        elseif mode == "Ps_Us"
            converter[1, CONV_MODE] = 3.0  # Ps_Us模式的逆变器
        elseif mode == "Udc_Qs"
            converter[1, CONV_MODE] = 4.0  # Udc_Qs模式的逆变器
        elseif mode == "Udc_Us"
            converter[1, CONV_MODE] = 5.0  # Udc_Us模式的逆变器
        elseif mode == "Droop_Udc_Qs"
            converter[1, CONV_MODE] = 6.0  # Droop_Udc_Qs模式的逆变器
        elseif mode == "Droop_Udc_Us"
            converter[1, CONV_MODE] = 7.0  # Droop_Udc_Us模式的逆变器
        else
            @warn "逆变器 $i 的控制模式 $mode 未知或不支持，默认启用 Ps_Qs"
            converter[1, CONV_MODE] = 0.0  # 设置为默认值
        end

        # 计算交流侧功率
        p_ac = -inverter.p_mw 
        q_ac = -inverter.q_mvar 
        
        # 计算直流侧功率（考虑效率）
        efficiency = 1.0 - inverter.loss_percent   # 转换为小数
        
        if p_ac <= 0  # 交流侧输出功率，直流侧输入功率
            p_dc = -p_ac / efficiency  # 负值，表示直流侧消耗的功率
        else  # 交流侧输入功率，直流侧输出功率
            p_dc = -p_ac * efficiency  # 正值，表示直流侧输出的功率
        end

        converter[1,CONV_ACBUS] = inverter.bus_ac
        converter[1,CONV_DCBUS] = inverter.bus_dc
        converter[1,CONV_INSERVICE] = 1.0
        converter[1,CONV_P_AC] = p_ac
        converter[1,CONV_Q_AC] = q_ac
        converter[1,CONV_P_DC] = p_dc
        converter[1,CONV_EFF] = efficiency
        converter[1,CONV_DROOP_KP] = inverter.droop_kv
        converters = vcat(converters, converter)
        
        # 获取交流和直流母线行索引
        ac_bus_row = findfirst(x -> x == inverter.bus_ac, jpc.busAC[:, 1])
        dc_bus_row = findfirst(x -> x == inverter.bus_dc, jpc.busDC[:, 1])
        
        # 根据控制模式决定是否修改负荷信息
        if mode =="δs_Us"
            # δs_Us模式：不做任何修改
        elseif mode == "Ps_Qs"
            # Ps_Qs模式：既修改交流侧的有功和无功，也修改直流侧的有功
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, PD] += p_ac  # PD - 有功负荷(MW)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - 无功负荷(MVAr)
            end
            
            if !isnothing(dc_bus_row)
                jpc.busDC[dc_bus_row, PD] += p_dc  # PD - 有功负荷(MW)
            end
            
            # 处理交流侧负荷
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # 创建一行
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # 负荷编号
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # 母线编号
                new_load_ac[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_ac[1, LOAD_PD] = p_ac                        # 有功功率(MW)
                new_load_ac[1, LOAD_QD] = q_ac                        # 无功功率(MVAr)
                # 逆变器默认为恒功率负荷
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_ac
                    # 获取原始负荷的功率和ZIP比例
                    orig_p = jpc.loadAC[idx, LOAD_PD]
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    orig_z_percent = jpc.loadAC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadAC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadAC[idx, LOADP_PERCENT]
                    
                    # 计算新的总功率
                    new_p = orig_p + p_ac
                    new_q = orig_q + q_ac
                    
                    # 重新计算ZIP比例（加权平均）
                    # 避免除以零的情况
                    if new_p != 0
                        # 原始负荷的权重
                        w_orig = abs(orig_p) / abs(new_p)
                        # 逆变器负荷的权重（默认为恒功率）
                        w_inv = abs(p_ac) / abs(new_p)
                        
                        # 计算新的ZIP比例
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # 确保比例总和为1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # 如果总和为0，设置为默认值
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # 如果新的总功率为0，保持原始ZIP比例
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # 更新负荷矩阵
                    jpc.loadAC[idx, LOAD_PD] = new_p
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    jpc.loadAC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadAC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadAC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
            
            # 处理直流侧负荷
            existing_load_indices_dc = findall(x -> x == inverter.bus_dc, jpc.loadDC[:, 2])
            
            if isempty(existing_load_indices_dc)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_dc_load_count += 1
                new_load_dc = zeros(1, num_cols_dc)  # 创建一行
                new_load_dc[1, LOAD_I] = nld_dc + new_dc_load_count  # 负荷编号
                new_load_dc[1, LOAD_CND] = inverter.bus_dc             # 母线编号
                new_load_dc[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_dc[1, LOAD_PD] = p_dc                        # 有功功率(MW)
                new_load_dc[1, LOAD_QD] = 0.0                         # 无功功率(MVAr)
                # 直流系统无无功
                # 直流侧默认为恒功率负荷
                new_load_dc[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_dc[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_dc[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_dc = vcat(new_loads_dc, new_load_dc)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_dc
                    # 获取原始负荷的功率和ZIP比例
                    orig_p = jpc.loadDC[idx, 4]
                    orig_z_percent = jpc.loadDC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadDC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadDC[idx, LOADP_PERCENT]
                    
                    # 计算新的总功率
                    new_p = orig_p + p_dc
                    
                    # 重新计算ZIP比例（加权平均）
                    # 避免除以零的情况
                    if new_p != 0
                        # 原始负荷的权重
                        w_orig = abs(orig_p) / abs(new_p)
                        # 逆变器负荷的权重（默认为恒功率）
                        w_inv = abs(p_dc) / abs(new_p)
                        
                        # 计算新的ZIP比例
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # 确保比例总和为1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # 如果总和为0，设置为默认值
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # 如果新的总功率为0，保持原始ZIP比例
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # 更新负荷矩阵
                    jpc.loadDC[idx, LOAD_PD] = new_p
                    jpc.loadDC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadDC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadDC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
        elseif mode == "Ps_Us"
            # Ps_Us模式：不做任何修改
           
        elseif mode == "Udc_Qs"
            # Udc_Qs模式：只修改交流侧的无功
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - 无功负荷(MVAr)
                # 不修改有功
            end
            
            # 处理交流侧负荷 - 只修改无功
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # 创建一行
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # 负荷编号
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # 母线编号
                new_load_ac[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_ac[1, LOAD_PD] = 0.0                         # 有功功率(MW) - 不修改
                new_load_ac[1, LOAD_QD] = q_ac                        # 无功功率(MVAr)
                # 逆变器默认为恒功率负荷
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_ac
                    # 获取原始负荷的功率
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    
                    # 计算新的总功率
                    new_q = orig_q + q_ac
                    
                    # 更新负荷矩阵 - 只修改无功
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    # 不修改ZIP比例，因为它们主要与有功相关
                end
            end
        elseif mode == "Udc_Us"
            # Udc_Us模式：不修改交流侧和直流侧负荷
            # 不做任何修改
        elseif mode == "Droop_Udc_Qs"
             # Udc_Qs模式：只修改交流侧的无功
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - 无功负荷(MVAr)
                # 不修改有功
            end
            
            # 处理交流侧负荷 - 只修改无功
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # 创建一行
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # 负荷编号
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # 母线编号
                new_load_ac[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_ac[1, LOAD_PD] = 0.0                         # 有功功率(MW) - 不修改
                new_load_ac[1, LOAD_QD] = q_ac                        # 无功功率(MVAr)
                # 逆变器默认为恒功率负荷
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_ac
                    # 获取原始负荷的功率
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    
                    # 计算新的总功率
                    new_q = orig_q + q_ac
                    
                    # 更新负荷矩阵 - 只修改无功
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    # 不修改ZIP比例，因为它们主要与有功相关
                end
            end
        elseif mode == "Droop_Udc_Us"
            # Droop_Udc_Us模式：不修改交流侧和直流侧负荷
            # 不做任何修改
        else
            # 未知模式：同时修改交流侧和直流侧负荷
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, PD] += p_ac  # PD - 有功负荷(MW)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - 无功负荷(MVAr)
            end
            
            if !isnothing(dc_bus_row)
                jpc.busDC[dc_bus_row, PD] += p_dc  # PD - 有功负荷(MW)
            end
            
            # 处理交流侧负荷
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # 创建一行
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # 负荷编号
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # 母线编号
                new_load_ac[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_ac[1, LOAD_PD] = p_ac                        # 有功功率(MW)
                new_load_ac[1, LOAD_QD] = q_ac                        # 无功功率(MVAr)
                # 逆变器默认为恒功率负荷
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_ac
                    # 获取原始负荷的功率和ZIP比例
                    orig_p = jpc.loadAC[idx, LOAD_PD]
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    orig_z_percent = jpc.loadAC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadAC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadAC[idx, LOADP_PERCENT]
                    
                    # 计算新的总功率
                    new_p = orig_p + p_ac
                    new_q = orig_q + q_ac
                    
                    # 重新计算ZIP比例（加权平均）
                    # 避免除以零的情况
                    if new_p != 0
                        # 原始负荷的权重
                        w_orig = abs(orig_p) / abs(new_p)
                        # 逆变器负荷的权重（默认为恒功率）
                        w_inv = abs(p_ac) / abs(new_p)
                        
                        # 计算新的ZIP比例
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # 确保比例总和为1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # 如果总和为0，设置为默认值
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # 如果新的总功率为0，保持原始ZIP比例
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # 更新负荷矩阵
                    jpc.loadAC[idx, LOAD_PD] = new_p
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    jpc.loadAC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadAC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadAC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
            
            # 处理直流侧负荷
            existing_load_indices_dc = findall(x -> x == inverter.bus_dc, jpc.loadDC[:, 2])
            
            if isempty(existing_load_indices_dc)
                # 如果不存在连接到相同节点的负荷，创建新的负荷记录
                new_dc_load_count += 1
                new_load_dc = zeros(1, num_cols_dc)  # 创建一行
                new_load_dc[1, LOAD_I] = nld_dc + new_dc_load_count  # 负荷编号
                new_load_dc[1, LOAD_CND] = inverter.bus_dc             # 母线编号
                new_load_dc[1, LOAD_STATUS] = 1.0                         # 状态(1=投运)
                new_load_dc[1, LOAD_PD] = p_dc                        # 有功功率(MW)
                new_load_dc[1, LOAD_QD] = 0.0                         # 无功功率(MVAr)
                # 直流系统无无功
                # 直流侧默认为恒功率负荷
                new_load_dc[1, LOADZ_PERCENT] = 0.0                         # 恒阻抗比例
                new_load_dc[1, LOADI_PERCENT] = 0.0                         # 恒电流比例
                new_load_dc[1, LOADP_PERCENT] = 1.0                         # 恒功率比例
                
                # 添加到新负荷矩阵
                new_loads_dc = vcat(new_loads_dc, new_load_dc)
            else
                # 如果存在连接到相同节点的负荷，更新这些负荷
                for idx in existing_load_indices_dc
                    # 获取原始负荷的功率和ZIP比例
                    orig_p = jpc.loadDC[idx, 4]
                    orig_z_percent = jpc.loadDC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadDC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadDC[idx, LOADP_PERCENT]
                    
                    # 计算新的总功率
                    new_p = orig_p + p_dc
                    
                    # 重新计算ZIP比例（加权平均）
                    # 避免除以零的情况
                    if new_p != 0
                        # 原始负荷的权重
                        w_orig = abs(orig_p) / abs(new_p)
                        # 逆变器负荷的权重（默认为恒功率）
                        w_inv = abs(p_dc) / abs(new_p)
                        
                        # 计算新的ZIP比例
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # 确保比例总和为1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # 如果总和为0，设置为默认值
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # 如果新的总功率为0，保持原始ZIP比例
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # 更新负荷矩阵
                    jpc.loadDC[idx, LOAD_PD] = new_p
                    jpc.loadDC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadDC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadDC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
        end

        # 根据inverter的控制模式，对JPC进行相应的处理
        if mode == "δs_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 3.0  # 设置为平衡节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # 设置为P节点
            inverter_gens_ac = zeros(1, 32)  # 假设发电机信息有26列
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: 发电机连接的节点编号
            inverter_gens_ac[1, 2] = -p_ac  # PG: 初始有功功率输出(MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: 初始无功功率输出(MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens_ac[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: 电压设定值(p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens_ac[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # 添加到genAC
        elseif mode == "Ps_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # 设置为PQ节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # 设置为P节点
        elseif mode == "Ps_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # 设置为PV节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # 设置为P节点

            inverter_gens_ac = zeros(1, 26)  # 假设发电机信息有32列
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: 发电机连接的节点编号
            inverter_gens_ac[1, 2] = -p_ac  # PG: 初始有功功率输出(MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: 初始无功功率输出(MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens_ac[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: 电压设定值(p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens_ac[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # 添加到genAC
        elseif mode == "Udc_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # 设置为PQ节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # 设置为平衡节点
            # 创建inverter 发电机信息
            inverter_gens = zeros(1, 32)  # 假设发电机信息有32列
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: 发电机连接的节点编号
            inverter_gens[1, 2] = p_dc  # PG: 初始有功功率输出(MW)
            inverter_gens[1, 3] = 0.0  # QG: 初始无功功率输出(MVAR)，直流系统通常为0
            inverter_gens[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens[1, 6] = inverter.vm_dc_pu  # VG: 电压设定值(p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # 添加到genDC
        elseif mode == "Udc_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # 设置为PV节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # 设置为平衡节点
            # 创建inverter 发电机信息
            inverter_gens = zeros(1, 32)  # 假设发电机信息有32列
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: 发电机连接的节点编号
            inverter_gens[1, 2] = -p_dc  # PG: 初始有功功率输出(MW)
            inverter_gens[1, 3] = 0.0  # QG: 初始无功功率输出(MVAR)，直流系统通常为0
            inverter_gens[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens[1, 6] = inverter.vm_dc_pu  # VG: 电压设定值(p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # 添加到genDC

            inverter_gens_ac = zeros(1, 26)  # 假设发电机信息有32列
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: 发电机连接的节点编号
            inverter_gens_ac[1, 2] = -p_ac  # PG: 初始有功功率输出(MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: 初始无功功率输出(MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens_ac[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: 电压设定值(p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens_ac[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # 添加到genAC
        elseif mode == "Droop_Udc_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # 设置为PQ节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # 设置为平衡节点
            # 创建inverter 发电机信息
            inverter_gens = zeros(1, 32)  # 假设发电机信息有32列
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: 发电机连接的节点编号
            inverter_gens[1, 2] = p_dc  # PG: 初始有功功率输出(MW)
            inverter_gens[1, 3] = 0.0  # QG: 初始无功功率输出(MVAR)，直流系统通常为0
            inverter_gens[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens[1, 6] = 1.0  # VG: 电压设定值(p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # 添加到genDC
        elseif mode == "Droop_Udc_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # 设置为PV节点
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # 设置为平衡节点
            # 创建inverter 发电机信息
            inverter_gens = zeros(1, 32)  # 假设发电机信息有32列
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: 发电机连接的节点编号
            inverter_gens[1, 2] = -p_dc  # PG: 初始有功功率输出(MW)
            inverter_gens[1, 3] = 0.0  # QG: 初始无功功率输出(MVAR)，直流系统通常为0
            inverter_gens[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens[1, 6] = 1.0  # VG: 电压设定值(p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # 添加到genDC

            inverter_gens_ac = zeros(1, 26)  # 假设发电机信息有32列
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: 发电机连接的节点编号
            inverter_gens_ac[1, 2] = -p_ac  # PG: 初始有功功率输出(MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: 初始无功功率输出(MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: 最大无功功率输出
            inverter_gens_ac[1, 5] = 0.0  # QMIN: 最小无功功率输出
            inverter_gens_ac[1, 6] = 1.0  # VG: 电压设定值(p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: 发电机基准功率(MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: 发电机状态(1=投运)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: 最大有功功率输出(MW)，放电功率
            inverter_gens_ac[1, 10] = 0.0  # PMIN: 最小有功功率输出(MW)，充电功率
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # 添加到genAC
        end
    end


    
    # 将新的负荷添加到现有负荷矩阵中
    if size(new_loads_ac, 1) > 0
        if !isempty(jpc.loadAC)
            jpc.loadAC = vcat(jpc.loadAC, new_loads_ac)
        else
            jpc.loadAC = new_loads_ac
        end
    end
    
    if size(new_loads_dc, 1) > 0
        if !isempty(jpc.loadDC)
            jpc.loadDC = vcat(jpc.loadDC, new_loads_dc)
        else
            jpc.loadDC = new_loads_dc
        end
    end

    jpc.converter = converters

    return jpc
end


# function JPC_add_grid_external_sc_impedance(case::JuliaPowerCase, jpc::JPC)    
#     # 如果没有外部电网，返回空数组
#     if isempty(case.ext_grids)
#         return Float64[], Float64[]
#     end
    
#     # 收集所有外部电网的数据
#     external_buses = Int[]
#     Y_grid_real = Float64[]
#     Y_grid_imag = Float64[]
    
#     # 遍历所有外部电网
#     for ext_grid in case.ext_grids
#         # 获取外部电网连接的母线ID
#         external_bus = ext_grid.bus
        
#         # 计算外部电网的短路阻抗
#         c = 1.1
#         s_sc = ext_grid.s_sc_max_mva / jpc.baseMVA
#         rx = ext_grid.rx_max
#         z_grid = c / (s_sc / 3)
#         x_grid = z_grid / sqrt(1 + rx^2)
#         r_grid = x_grid * rx
        
#         # 计算导纳
#         Y_grid = 1 / (r_grid + 1im * x_grid)
        
#         push!(external_buses, external_bus)
#         push!(Y_grid_real, real(Y_grid))
#         push!(Y_grid_imag, imag(Y_grid))
#     end
    
#     # 对相同母线的导纳值进行汇总
#     buses, gs, bs = PowerFlow.sum_by_group(external_buses, Y_grid_real, Y_grid_imag)
    
#     # 更新jpc中的母线数据
#     for (i, bus_id) in enumerate(buses)
#         # 找到对应母线在jpc.busAC中的索引
#         bus_idx = findfirst(row -> row[1] == bus_id, jpc.busAC[:, 1])
        
#         if !isnothing(bus_idx)
#             # 更新母线的电导和电纳值
#             jpc.busAC[bus_idx, GS] = gs[i] * jpc.baseMVA
#             jpc.busAC[bus_idx, BS] = bs[i] * jpc.baseMVA
#         end
#     end
    
#     # 返回计算出的电导和电纳值（乘以基准功率）
#     return gs .* jpc.baseMVA, bs .* jpc.baseMVA, jpc
# end




