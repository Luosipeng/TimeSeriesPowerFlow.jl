
"""
    resolve_node_mapping(node_id, node_merge_map)

Resolves the final mapping of a node by traversing through a node merge map.

This function follows the chain of node mappings to find the final destination node.
It handles multiple levels of redirection by iteratively looking up each node ID
in the mapping until it finds a node that doesn't have a further mapping.

# Arguments
- `node_id`: The initial node ID to resolve
- `node_merge_map`: A dictionary mapping source nodes to destination nodes

# Returns
The final resolved node ID after following all mappings
"""
function resolve_node_mapping(node_id, node_merge_map)
    while haskey(node_merge_map, node_id)
        node_id = node_merge_map[node_id]
    end
    return node_id
end

"""
    merge_virtual_nodes(case::JuliaPowerCase)

Merges virtual nodes in a power system case and updates all connected elements.

This function identifies virtual nodes (nodes with "_虚拟节点" in their name), 
merges them with their connected real nodes, and updates all references in the system.
The process includes:
1. Identifying virtual nodes and their connections
2. Creating a mapping strategy for node merging
3. Updating all elements (lines, transformers, circuit breakers, loads, generators)
4. Removing virtual nodes and their associated circuit breakers
5. Merging loads connected to the same bus after node consolidation

# Arguments
- `case::JuliaPowerCase`: The original power system case

# Returns
A new JuliaPowerCase with virtual nodes merged and all references updated

# Note
- Virtual nodes are identified by the substring "_虚拟节点" in their names
- At least two connections are required for a virtual node to be merged
- Self-loops created during merging are disabled
- Loads at the same bus after merging are consolidated into a single load
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

"""
    JuliaPowerCase2Jpc(case::Utils.JuliaPowerCase)

Converts a JuliaPowerCase object to a JPC (Julia Power Case) format for power flow analysis.

This function performs a complete transformation of the power system case data,
processing all components and preparing them for power flow calculations.
The conversion process includes:

1. Merging virtual nodes to ensure a proper network topology
2. Creating and populating a new JPC object with the following data:
   - Basic parameters (base MVA)
   - AC buses
   - DC buses
   - AC branches (lines and transformers)
   - DC branches
   - AC generators
   - Battery systems (as DC generators)
   - Battery state of charge data
   - AC loads
   - DC loads
   - PV array data
   - Inverters (AC/DC converters)
   - AC-connected PV systems

# Arguments
- `case::Utils.JuliaPowerCase`: The original power system case data

# Returns
A fully populated JPC object ready for power flow analysis
"""
function JuliaPowerCase2Jpc(case::Utils.JuliaPowerCase)
    # 1. merge virtual nodes
    case = merge_virtual_nodes(case)
    
    # 2. create a new JPC object
    jpc = JPC()
    
    # 3. setting basic parameters
    jpc.baseMVA = case.baseMVA
    
    # 4. setting ac node data
    JPC_buses_process(case, jpc)

    # 5. setting dc node data
    JPC_dcbuses_process(case, jpc)
    
    # 6. setting ac branches data
    JPC_branches_process(case, jpc)
    
    # 7. setting dc branches data
    JPC_dcbranches_process(case, jpc)

    # 8. setting ac generators data
    JPC_gens_process(case, jpc)

    # 9. setting battery systems data
    JPC_battery_gens_process(case, jpc)

    # 10. setting battery state of charge data
    JPC_battery_soc_process(case, jpc)
    
    # 11. setting ac loads data
    JPC_loads_process(case, jpc)

    # 12. setting dc loads data
    JPC_dcloads_process(case, jpc)

    # 13. setting pv array data 
    JPC_pv_process(case, jpc)

    # 14. setting inverters data
    JPC_inverters_process(case, jpc)

    # 15. setting ac-connected pv systems data
    JPC_ac_pv_system_process(case, jpc)

    return jpc

end

"""
    JPC_buses_process(case::JuliaPowerCase, jpc::JPC)

Processes AC bus data from a JuliaPowerCase and converts it to JPC format.

This function extracts AC bus information from the input case, transforms it into 
the standard matrix format required by JPC, and stores it in the JPC object.
Each bus is represented as a row in the matrix with the following columns:

1. Bus ID
2. Bus type (all initialized as PQ nodes with value 1.0)
3. Active power demand (PD) in MW (initialized to 0.0)
4. Reactive power demand (QD) in MVAR (initialized to 0.0)
5. Active power shunt conductance (GS) in MW (initialized to 0.0)
6. Reactive power shunt susceptance (BS) in MVAR (initialized to 0.0)
7. Area ID
8. Voltage magnitude (VM) in p.u. (initialized to 1.0)
9. Voltage angle (VA) in degrees (initialized to 0.0)
10. Base voltage (VN) in kV
11. Zone ID
12. Maximum voltage magnitude in p.u.
13. Minimum voltage magnitude in p.u.

# Arguments
- `case::JuliaPowerCase`: The original power system case containing bus data
- `jpc::JPC`: The JPC object where the processed bus data will be stored

# Returns
The updated JPC object with the busAC field populated
"""
function JPC_buses_process(case::JuliaPowerCase, jpc::JPC)
    # obtain AC node data and convert to JPC format
    buses = deepcopy(case.busesAC)
    
    # crate a matrix to store bus data
    num_buses = length(buses)
    bus_matrix = zeros(num_buses, 13)
    
    for (i, bus) in enumerate(buses)
        # setting initial voltage values (based on index)
        vm = 1.0
        va = 0.0
        
        # fill each row of the matrix
        bus_matrix[i, :] = [
            bus.bus_id,      # node ID
            1.0,             # node type (all set to PQ node)
            0.0,             # PD (MW) active load (MW)
            0.0,             # QD (MVAR) reactive load (MVAR)
            0.0,             # GS (MW) active generation (MW)
            0.0,             # BS (MVAR) reactive generation (MVAR)
            bus.area_id,     # area ID
            vm,              # node voltage magnitude (p.u.)
            va,              # node voltage angle (degrees)
            bus.vn_kv,       # node voltage base (kV)
            bus.zone_id,     # zone ID
            bus.max_vm_pu,   # maximum voltage magnitude (p.u.)
            bus.min_vm_pu,   # minimum voltage magnitude (p.u.
        ]
    end
    
    # store all bus data in the busAC field of jpc
    jpc.busAC = bus_matrix
    
    return jpc
end

"""
    JPC_dcbuses_process(case::JuliaPowerCase, jpc::JPC)

Processes DC bus data from a JuliaPowerCase and converts it to JPC format.

This function extracts DC bus information from the input case, transforms it into 
the standard matrix format required by JPC, and stores it in the JPC object.
Each DC bus is represented as a row in the matrix with the following columns:

1. Bus ID
2. Bus type (all initialized as PQ nodes with value 1.0)
3. Active power demand (PD) in MW (initialized to 0.0)
4. Reactive power demand (QD) in MVAR (initialized to 0.0)
5. Active power shunt conductance (GS) in MW (initialized to 0.0)
6. Reactive power shunt susceptance (BS) in MVAR (initialized to 0.0)
7. Area ID
8. Voltage magnitude (VM) in p.u. (initialized to 1.0)
9. Voltage angle (VA) in degrees (initialized to 0.0)
10. Base voltage (VN) in kV
11. Zone ID
12. Maximum voltage magnitude in p.u.
13. Minimum voltage magnitude in p.u.

After processing the basic DC bus data, this function also calls JPC_battery_bus_process
to handle battery-specific bus information.

# Arguments
- `case::JuliaPowerCase`: The original power system case containing DC bus data
- `jpc::JPC`: The JPC object where the processed DC bus data will be stored

# Returns
The updated JPC object with the busDC field populated and battery bus data processed
"""
function JPC_dcbuses_process(case, jpc)
    # obtain DC node data and convert to JPC format
    dcbuses = deepcopy(case.busesDC)
    
    # create a matrix to store DC bus data
    num_dcbuses = length(dcbuses)
    dcbus_matrix = zeros(num_dcbuses, 13)
    
    for (i, dcbus) in enumerate(dcbuses)
        # setting initial voltage values (based on index)
        vm = 1.0
        va = 0.0
        
        # fill each row of the matrix
        dcbus_matrix[i, :] = [
            dcbus.bus_id,      # node ID
            1.0,               # node type (all set to PQ node)
            0.0,               # PD (MW) active load (MW)
            0.0,               # QD (MVAR) reactive load (MVAR)
            0.0,               # GS (MW) active generation (MW)
            0.0,               # BS (MVAR) reactive generation (MVAR)
            dcbus.area_id,     # area ID
            vm,                # node voltage magnitude (p.u.)
            va,                # node voltage angle (degrees)
            dcbus.vn_kv,       # node voltage base (kV)
            dcbus.zone_id,     # zone ID
            dcbus.max_vm_pu,   # maximum voltage magnitude (p.u.)
            dcbus.min_vm_pu,   # minimum voltage magnitude (p.u.
        ]
    end
    
    # store all DC bus data in the busDC field of jpc
    jpc.busDC = dcbus_matrix

    jpc = JPC_battery_bus_process(case, jpc)
    
    return jpc
end

"""
    JPC_battery_bus_process(case::JuliaPowerCase, jpc::JPC)

Processes battery storage data from a JuliaPowerCase and creates virtual DC buses for them in the JPC format.

This function extracts battery storage information from the input case, creates virtual DC buses
to represent the batteries in the power system model, and appends these buses to the existing
DC bus matrix in the JPC object.

Each battery is represented as a virtual DC bus with the following characteristics:
1. Bus ID assigned sequentially after existing DC buses
2. Bus type set to PV node (value 2.0)
3. No initial power demand or generation
4. Voltage magnitude initialized to 1.0 p.u.
5. Voltage angle initialized to 0.0 degrees
6. Base voltage set to the battery's open-circuit voltage
7. Area and zone IDs set to 1.0
8. Voltage limits set to 0.95 (min) and 1.05 (max) p.u.

# Arguments
- `case::JuliaPowerCase`: The original power system case containing battery storage data
- `jpc::JPC`: The JPC object where the virtual battery buses will be added

# Returns
The updated JPC object with virtual battery buses added to the busDC field
"""
function JPC_battery_bus_process(case::JuliaPowerCase, jpc::JPC)
    # process battery data and create virtual nodes for batteries
    batteries = deepcopy(case.storageetap)
    
    # obtain the current busDC data
    busDC = jpc.busDC
    current_size = size(busDC, 1)
    
    # create a matrix to store battery data
    num_batteries = length(batteries)
    battery_matrix = zeros(num_batteries, size(busDC, 2))
    
    for (i, battery) in enumerate(batteries)
        # setting initial voltage values (based on index)
        vm = 1.0
        va = 0.0
        vn_kv = battery.voc
        # create a new row for the battery node
        battery_row = zeros(1, size(busDC, 2))
        battery_row[1, :] = [
            current_size + i,    # distributed node ID (starting from the end of busDC)
            2.0,                 # node type (set to PV node)
            0.0,                 # PD (MW) active load (MW)
            0.0,                 # QD (MVAR) reactive load (MVAR)
            0.0,                 # GS (MW) active generation (MW)
            0.0,                 # BS (MVAR) reactive generation (MVAR)
            1.0,                 # area ID (set to 1.0)
            vm,                  # node voltage magnitude (p.u.)
            va,                  # node voltage angle (degrees)
            vn_kv,               # node voltage base (kV)
            1.0,                 # zone ID (set to 1.0)
            1.05,                 # maximum voltage magnitude (p.u.)
            0.95                  # minimum voltage magnitude (p.u.)
        ]
        
        # store the battery row in the matrix
        battery_matrix[i, :] = battery_row
    end
    
    # merge the battery data into the busDC field of jpc
    jpc.busDC = vcat(busDC, battery_matrix)
    
    # # save the battery data in the jpc structure
    # jpc.battery = battery_matrix

    return jpc
end

"""
    JPC_branches_process(case::JuliaPowerCase, jpc::JPC)

Processes AC branch data from a JuliaPowerCase and converts it to JPC format.

This function handles the calculation and conversion of various branch elements 
in the power system to the format required by JPC. It processes:

1. AC transmission lines by calling calculate_line_parameters()
2. Two-winding transformers by calling calculate_transformer2w_parameters()

The function is structured to potentially handle different sequence components
(positive, negative, zero), though the conditional logic for sequence selection
is currently commented out. In its current implementation, it processes only the
positive sequence components.

Note: The commented code suggests that for sequence 1 or 2 (positive or negative sequence),
the function processes lines and transformers with their standard parameters, while for
sequence 0 (zero sequence), it would call a different function (calculate_branch_JPC_zero).

# Arguments
- `case::JuliaPowerCase`: The original power system case containing branch data
- `jpc::JPC`: The JPC object where the processed branch data will be stored

# Returns
The updated JPC object with branch data processed (implicitly, as the jpc object is modified in place)
"""
function JPC_branches_process(case::JuliaPowerCase, jpc::JPC)
    # if sequence == 1||sequence == 2
        # process AC branch data, converting to JPC format
        calculate_line_parameters(case::JuliaPowerCase, jpc)
        # process transformer data, converting to JPC format
        calculate_transformer2w_parameters(case::JuliaPowerCase, jpc)
        # process transformer 3-winding data, converting to JPC format
    # else
    #     # process AC branch data, converting to JPC format
    #     calculate_branch_JPC_zero(case::JuliaPowerCase, jpc)
    # end
end

"""
    JPC_dcbranches_process(case::JuliaPowerCase, jpc::JPC)

Processes DC branch data from a JuliaPowerCase and converts it to JPC format.

This function extracts DC branch information from the input case, transforms it into 
the standard matrix format required by JPC, and stores it in the JPC object.
Each DC branch is represented as a row in the matrix with the following columns:

1. From bus index (F_BUS)
2. To bus index (T_BUS)
3. Branch resistance in p.u. (BR_R) - calculated from length and ohms/km
4. Branch reactance in p.u. (BR_X) - set to 0 for DC branches
5. Branch rating in MVA (RATE_A) - calculated from max current or set to default
6. Branch status (BR_STATUS) - 1.0 if in service, 0.0 if not
7-14. Various parameters including angle limits (ANGMIN, ANGMAX)

The function performs the following steps:
1. Initializes a matrix to store DC branch data
2. For each DC branch in the case:
   - Identifies the from and to bus indices
   - Calculates the base impedance using the from bus base voltage
   - Converts resistance from ohms to per unit values
   - Sets the branch parameters including status and ratings
   - Sets angle limits to -360° and 360°
3. Appends the processed branch data to the JPC object
4. Processes battery branches by calling JPC_battery_branch_process

# Arguments
- `case::JuliaPowerCase`: The original power system case containing DC branch data
- `jpc::JPC`: The JPC object where the processed DC branch data will be stored

# Returns
The updated JPC object with the branchDC field populated and battery branch data processed
"""
function JPC_dcbranches_process(case::JuliaPowerCase, jpc::JPC)
    # process DC branch data from a JuliaPowerCase and convert it to JPC format
    nbr = length(case.branchesDC)
    branch = zeros(nbr, 14)
    dclines = case.branchesDC

    for (i, dcline) in enumerate(dclines)
        # obtain the from and to bus indices
        from_bus_idx = dcline.from_bus
        to_bus_idx = dcline.to_bus
        
        # obtain the base voltage for the from bus
        basekv = jpc.busDC[from_bus_idx, BASE_KV]
        
        # calculate the base impedance
        baseR = (basekv^2) / case.baseMVA
        
        # calculate the per unit resistance and reactance
        r_pu = 2 * dcline.length_km * dcline.r_ohm_per_km / baseR
        x_pu = 0
        
        # fill the branch matrix
        branch[i, F_BUS] = from_bus_idx
        branch[i, T_BUS] = to_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        
        # setting the branch parameters
        if hasfield(typeof(dcline), :max_i_ka)
            branch[i, RATE_A] = dcline.max_i_ka * basekv * sqrt(3)  # based capacity in MVA
        else
            branch[i, RATE_A] = 100.0  # default value if not specified
        end
        
        # setting the branch status
        branch[i, BR_STATUS] = dcline.in_service ? 1.0 : 0.0
        
        # setting the branch angle limits
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end
    # add the processed branch data to the jpc structure
    if isempty(jpc.branchDC)
        jpc.branchDC = branch
    else
        jpc.branchDC = [jpc.branchDC; branch]
    end
    
    # process battery branches
    jpc = JPC_battery_branch_process(case, jpc)

    return jpc
end

"""
    JPC_battery_branch_process(case::JuliaPowerCase, jpc::JPC)

Processes battery storage branches from a JuliaPowerCase and converts them to JPC format.

This function creates DC branch connections between virtual battery buses and their 
corresponding actual buses in the power system. Each battery is represented as a DC branch
with the following characteristics:

1. From bus: Virtual battery bus (created in JPC_battery_bus_process)
2. To bus: Actual bus where the battery is physically connected
3. Branch resistance: Calculated from battery internal resistance in per unit
4. Branch reactance: Set to 0 (DC branches have no reactance)
5. Branch rating: Based on battery package size and open-circuit voltage
6. Branch status: 1.0 if in service, 0.0 if not
7. Angle limits: Set to -360° and 360°

The function performs the following steps:
1. Creates a deep copy of the battery storage data
2. If no batteries exist, returns the JPC object unchanged
3. Initializes a matrix to store battery branch data
4. For each battery:
   - Identifies the actual bus where the battery is connected
   - Calculates the virtual bus ID for the battery
   - Retrieves the base voltage from the actual bus
   - Calculates the per unit resistance using the base impedance
   - Sets all branch parameters including status and ratings
5. Appends the processed battery branch data to the JPC object's branchDC field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing battery storage data
- `jpc::JPC`: The JPC object where the battery branch data will be added

# Returns
The updated JPC object with battery branches added to the branchDC field
"""
function JPC_battery_branch_process(case::JuliaPowerCase, jpc::JPC)
    # process battery branches from a JuliaPowerCase and convert them to JPC format
    batteries = deepcopy(case.storageetap)
    num_batteries = length(batteries)
    
    # if there are no batteries, return the jpc object as is
    if num_batteries == 0
        return jpc
    end
    
    # create a matrix to store battery branch data
    battery_branches = zeros(num_batteries, 14)
    
    # obtain the size of the existing busDC matrix
    busDC_size = size(jpc.busDC, 1) - num_batteries
    
    for (i, battery) in enumerate(batteries)
        # obtain the actual bus ID where the battery is connected
        actual_bus = battery.bus
        
        # calculate the virtual bus ID for the battery
        virtual_bus = busDC_size + i
        
        # obtain the base voltage for the actual bus
        basekv = 0.0
        for j in 1:size(jpc.busDC, 1)
            if jpc.busDC[j, 1] == actual_bus
                basekv = jpc.busDC[j, BASE_KV]
                break
            end
        end
        
        # calculate the base impedance
        baseR = (basekv^2) / case.baseMVA
        
        # calculate the per unit resistance and reactance
        # r_pu = battery.ra / baseR
        # r_pu = 0.0242/baseR  # calculated from battery parameters
        r_pu = 0.0252115/baseR  # assumed value based on battery parameters
        x_pu = 0  # dc branches have no reactance
        
        # fill the battery branch matrix
        battery_branches[i, F_BUS] = virtual_bus       # virtual bus ID
        battery_branches[i, T_BUS] = actual_bus        # actual bus ID
        battery_branches[i, BR_R] = r_pu               # per unit resistance
        battery_branches[i, BR_X] = x_pu               # per unit reactance
        
        # setting the branch parameters
        # assumed rated capacity based on battery parameters
        rated_capacity = battery.package * battery.voc  # simplified calculation
        battery_branches[i, RATE_A] = rated_capacity
        
        # setting the branch status
        battery_branches[i, BR_STATUS] = battery.in_service ? 1.0 : 0.0
        
        # setting the branch angle limits
        battery_branches[i, ANGMIN] = -360.0
        battery_branches[i, ANGMAX] = 360.0
    end
    
    # add the processed battery branch data to the jpc structure
    if isempty(jpc.branchDC)
        jpc.branchDC = battery_branches
    else
        jpc.branchDC = [jpc.branchDC; battery_branches]
    end
    
    return jpc
end

"""
    JPC_battery_soc_process(case::JuliaPowerCase, jpc::JPC)

Processes battery state of charge (SOC) data from a JuliaPowerCase and converts it to JPC format.

This function extracts battery storage information from the input case and creates two data structures:
1. A battery SOC matrix containing parameters related to battery energy storage capabilities
2. Load entries in the loadDC matrix to represent battery connections to the power system

For each battery, the following SOC parameters are stored:
1. Bus ID where the battery is connected
2. Power capacity in MW (maximum charge/discharge rate)
3. Energy capacity in MWh
4. Initial state of charge (SOC)
5. Minimum allowable SOC
6. Maximum allowable SOC
7. Round-trip efficiency (0.0-1.0)
8. Status indicator (1.0 if in service, 0.0 if not)

Additionally, for each battery, the function:
1. Creates a corresponding load entry in the loadDC matrix
2. Sets the load ID, bus index, and status parameters
3. Initializes power values to zero

The function performs the following steps:
1. Creates a deep copy of the battery storage data
2. If no batteries exist, returns the JPC object unchanged
3. Initializes a matrix to store battery SOC data
4. Populates the SOC matrix with battery parameters
5. For each battery, creates a corresponding load entry in the loadDC matrix
6. Adds the processed battery SOC data to the JPC object's storage field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing battery storage data
- `jpc::JPC`: The JPC object where the battery SOC data will be stored

# Returns
The updated JPC object with battery SOC data in the storage field and corresponding
load entries in the loadDC field
"""
function JPC_battery_soc_process(case::JuliaPowerCase, jpc::JPC)
    # process battery state of charge (SOC) data from a JuliaPowerCase and convert it to JPC format
    batteries = deepcopy(case.storages)
    num_batteries = length(batteries)
    
    # if there are no batteries, return the jpc object as is
    if num_batteries == 0
        return jpc
    end
    
    # create a matrix to store battery SOC data
    battery_soc = zeros(num_batteries, 8)  
    
    for (i, battery) in enumerate(batteries)
        battery_soc[i, 1] = battery.bus  # battery connected bus ID
        battery_soc[i, 2] = battery.power_capacity_mw   # battery charge/discharge power maximum (MW)
        battery_soc[i, 3] = battery.energy_capacity_mwh  # battery energy capacity (MWh)
        battery_soc[i, 4] = battery.soc_init  # battery initial state of charge (SOC)
        battery_soc[i, 5] = battery.min_soc  # battery minimum state of charge (SOC)
        battery_soc[i, 6] = battery.max_soc  # battery maximum state of charge (SOC)
        battery_soc[i, 7] = battery.efficiency  # battery efficiency (0.0-1.0)
        battery_soc[i, 8] = battery.in_service ? 1.0 : 0.0  # if battery is in service (1.0) or not (0.0)
    end
    for (i, battery) in enumerate(batteries)
        # obtain the bus ID where the battery is connected
        bus_id = battery.bus
        # search for the bus index in the busDC matrix
        bus_index = findfirst(x -> x[1] == bus_id, jpc.busDC[:, 1])
        jpc.busDC[bus_index, PD] -= 0.0
        loadDC = zeros(1, 8)  # create a matrix for load data
        nd = size(jpc.busDC, 1)
        loadDC[1, 1] = nd + 1  # setting the load ID
        loadDC[1, 2] = bus_index
        loadDC[1, 3] = 1 # inservice
        loadDC[1, 4] = 0.0
        loadDC[1, 5] = 0.0
        loadDC[1, 6] = 0.0  
        loadDC[1, 7] = 0.0
        loadDC[1, 8] = 1.0
        # add the load data to the jpc structure
        if isempty(jpc.loadDC)
            jpc.loadDC = loadDC
        else
            jpc.loadDC = [jpc.loadDC; loadDC]
        end
    end
    # add the processed battery SOC data to the jpc structure
    jpc.storage = battery_soc
    
    return jpc
end

"""
JPC_battery_soc_process(case::JuliaPowerCase, jpc::JPC)
Processes battery state of charge (SOC) data from a JuliaPowerCase and converts it to JPC format.
This function extracts battery storage information from the input case and creates two data structures:

A battery SOC matrix containing parameters related to battery energy storage capabilities
Load entries in the loadDC matrix to represent battery connections to the power system

For each battery, the following SOC parameters are stored:

Bus ID where the battery is connected
Power capacity in MW (maximum charge/discharge rate)
Energy capacity in MWh
Initial state of charge (SOC)
Minimum allowable SOC
Maximum allowable SOC
Round-trip efficiency (0.0-1.0)
Status indicator (1.0 if in service, 0.0 if not)

Additionally, for each battery, the function:

Creates a corresponding load entry in the loadDC matrix
Sets the load ID, bus index, and status parameters
Initializes power values to zero

The function performs the following steps:

Creates a deep copy of the battery storage data
If no batteries exist, returns the JPC object unchanged
Initializes a matrix to store battery SOC data
Populates the SOC matrix with battery parameters
For each battery, creates a corresponding load entry in the loadDC matrix
Adds the processed battery SOC data to the JPC object's storage field

Arguments

case::JuliaPowerCase: The original power system case containing battery storage data
jpc::JPC: The JPC object where the battery SOC data will be stored

Returns
The updated JPC object with battery SOC data in the storage field and corresponding
load entries in the loadDC field
"""
function calculate_line_parameters(case::JuliaPowerCase, jpc::JPC)
    # Process line data, convert to JPC format
    nbr = length(case.branchesAC)
    branch = zeros(nbr, 14)
    lines = case.branchesAC

    for (i, line) in enumerate(lines)
        # Get from and to bus numbers
        from_bus_idx = line.from_bus
        to_bus_idx = line.to_bus
        
        # Get base voltage (kV) of the from bus
        basekv = jpc.busAC[from_bus_idx, BASE_KV]
        
        # Calculate base impedance
        baseR = (basekv^2) / case.baseMVA
        
        # Consider parallel lines
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # Calculate per unit impedance
        r_pu = line.length_km * line.r_ohm_per_km / baseR / parallel
        x_pu = line.length_km * line.x_ohm_per_km / baseR / parallel
        
        # Calculate shunt susceptance (p.u.)
        b_pu = 2 * π * case.basef * line.length_km * line.c_nf_per_km * 1e-9 * baseR * parallel
        
        # Calculate shunt conductance (p.u.)
        g_pu = 0.0
        if hasfield(typeof(line), :g_us_per_km)
            g_pu = line.g_us_per_km * 1e-6 * baseR * line.length_km * parallel
        end
        
        # Fill branchAC matrix
        branch[i, F_BUS] = from_bus_idx
        branch[i, T_BUS] = to_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        branch[i, BR_B] = b_pu
        
        # Set rated capacity
        if hasfield(typeof(line), :max_i_ka)
            branch[i, RATE_A] = line.max_i_ka * basekv * sqrt(3)  # Rated capacity (MVA)
        else
            branch[i, RATE_A] = 100.0  # Default value
        end
        
        # Set branch status
        branch[i, BR_STATUS] = line.in_service ? 1.0 : 0.0
        
        # Set angle limits
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end

    jpc.branchAC = branch
end


"""
    calculate_transformer2w_parameters(case::JuliaPowerCase, jpc::JPC)

Processes two-winding transformer data from a JuliaPowerCase and converts it to JPC format.

This function extracts two-winding transformer information from the input case and creates branch entries
in the JPC object to represent transformers in the power system model.

For each transformer, the following parameters are calculated and stored:
1. From bus (high voltage side) and to bus (low voltage side) indices
2. Per-unit resistance and reactance values, adjusted to system base
3. Shunt susceptance (typically zero for transformers)
4. Tap ratio and phase shift angle
5. Power rating and operational status
6. Angle limits

The function performs the following steps:
1. Extracts two-winding transformer data from the input case
2. If no transformers exist, returns without modification
3. Creates a branch matrix to store transformer parameters
4. For each transformer:
   - Identifies the high and low voltage bus connections
   - Calculates impedance parameters based on transformer specifications
   - Adjusts values for the system base power
   - Accounts for parallel transformers if present
   - Sets default values for tap ratio, phase shift, and other parameters
5. Adds the transformer branch data to the JPC object's branchAC field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing transformer data
- `jpc::JPC`: The JPC object where the transformer branch data will be stored

# Returns
The updated JPC object with transformer data added to the branchAC field
"""
function calculate_transformer2w_parameters(case::JuliaPowerCase, jpc::JPC)
    # Process transformer data, convert to JPC format
    transformers = case.transformers_2w_etap
    nbr = length(transformers)
    
    if nbr == 0
        return  # If no transformers, return directly
    end
    
    # Create transformer branch matrix
    branch = zeros(nbr, 14)
    
    for (i, transformer) in enumerate(transformers)
        # Get high voltage and low voltage bus numbers
        hv_bus_idx = transformer.hv_bus
        lv_bus_idx = transformer.lv_bus
        
        # Get high voltage bus base voltage (kV)
        hv_basekv = jpc.busAC[hv_bus_idx, BASE_KV]
        
        # Calculate impedance parameters
        # Convert transformer impedance percentage to per unit value
        z_pu = transformer.z_percent
        x_r_ratio = transformer.x_r
        
        # Calculate resistance and reactance (considering base power conversion)
        s_ratio = transformer.sn_mva / case.baseMVA
        z_pu = z_pu / s_ratio  # Convert to system base
        
        r_pu = z_pu / sqrt(1 + x_r_ratio^2)
        x_pu = r_pu * x_r_ratio
        
        # Consider parallel transformers
        parallel = transformer.parallel
        if parallel > 1
            r_pu = r_pu / parallel
            x_pu = x_pu / parallel
        end
        
        # Fill branch matrix
        branch[i, F_BUS] = hv_bus_idx
        branch[i, T_BUS] = lv_bus_idx
        branch[i, BR_R] = r_pu
        branch[i, BR_X] = x_pu
        branch[i, BR_B] = 0.0  # Transformers typically have no shunt susceptance
        
        # Set tap ratio and phase shift
        branch[i, TAP] = 1.0  # Default tap ratio is 1.0
        branch[i, SHIFT] = 0.0  # Default phase shift angle is 0.0
        
        # Set rated capacity
        branch[i, RATE_A] = case.baseMVA 
        
        # Set branch status
        branch[i, BR_STATUS] = transformer.in_service ? 1.0 : 0.0
        
        # Set angle limits
        branch[i, ANGMIN] = -360.0
        branch[i, ANGMAX] = 360.0
    end
    
    # Add transformer branch data to JPC structure
    if isempty(jpc.branchAC)
        jpc.branchAC = branch
    else
        jpc.branchAC = [jpc.branchAC; branch]
    end
end


"""
    JPC_gens_process(case::JuliaPowerCase, jpc::JPC)

Process generator data from a JuliaPowerCase and convert it to JPC format.

This function processes three types of generation devices:
1. External grids (typically serving as slack/reference nodes)
2. Conventional generators (typically serving as PV nodes)
3. Static generators (typically serving as PQ nodes, but can be PV nodes if controllable)

For each generation device, the function:
- Extracts its parameters from the input case
- Converts them to the format required by JPC
- Assigns appropriate bus types based on generator characteristics
- Ensures proper handling of power outputs, limits, and ramp rates

The function performs the following steps:
1. Counts the number of each type of generation device
2. Creates a matrix to store generator data with 26 columns
3. Processes each type of generator and populates the matrix
4. Updates bus types in the JPC structure based on generator connections
5. Removes entries for out-of-service generators
6. Ensures at least one slack node exists in the system

# Arguments
- `case::JuliaPowerCase`: The original power system case containing generator data
- `jpc::JPC`: The JPC object where the processed generator data will be stored

# Returns
None. The JPC object is modified in-place.
"""
function JPC_gens_process(case::JuliaPowerCase, jpc::JPC)
    # Count different types of generation devices
    n_gen = length(case.gensAC)
    n_sgen = length(case.sgensAC)
    n_ext = length(case.ext_grids)
    
    # Calculate total number of generation devices
    total_gens = n_gen + n_sgen + n_ext
    
    if total_gens == 0
        return  # If no generation devices, return directly
    end
    
    # Create generator matrix with rows equal to number of generation devices and 26 columns
    gen_data = zeros(total_gens, 26)
    
    # Process external grids (usually serving as slack/reference nodes)
    for (i, ext) in enumerate(case.ext_grids)
        if !ext.in_service
            continue
        end
        
        bus_idx = ext.bus
        
        # Fill generator data
        gen_data[i, :] = [
            bus_idx,        # Generator connection bus number
            0.0,            # Active power output (MW)
            0.0,            # Reactive power output (MVAr)
            9999.0,         # Maximum reactive power output (MVAr)
            -9999.0,        # Minimum reactive power output (MVAr)
            ext.vm_pu,      # Voltage magnitude setpoint (p.u.)
            case.baseMVA,   # Generator base capacity (MVA)
            1.0,            # Generator status (1=in service, 0=out of service)
            9999.0,         # Maximum active power output (MW)
            -9999.0,        # Minimum active power output (MW)
            0.0,            # PQ capability curve lower end active power output (MW)
            0.0,            # PQ capability curve upper end active power output (MW)
            0.0,            # Minimum reactive power output at PC1 (MVAr)
            0.0,            # Maximum reactive power output at PC1 (MVAr)
            0.0,            # Minimum reactive power output at PC2 (MVAr)
            0.0,            # Maximum reactive power output at PC2 (MVAr)
            0.0,            # AGC ramp rate (MW/min)
            0.0,            # 10-minute reserve ramp rate (MW)
            0.0,            # 30-minute reserve ramp rate (MW)
            0.0,            # Reactive power ramp rate (MVAr/min)
            1.0,            # Area participation factor
            2.0,            # Generator model (2=polynomial cost model)
            0.0,            # Startup cost (dollars)
            0.0,            # Shutdown cost (dollars)
            3.0,            # Number of polynomial cost function coefficients
            0.0             # Cost function parameters (to be expanded later)
        ]
        
        # Update bus type to reference node (REF/slack node)
        jpc.busAC[bus_idx, 2] = 3  # 3 indicates REF node
    end
    
    # Process conventional generators (usually serving as PV nodes)
    offset = n_ext
    for (i, gen) in enumerate(case.gensAC)
        if !gen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = gen.bus
        
        # Calculate reactive power (if not directly provided)
        q_mvar = 0.0
        if hasfield(typeof(gen), :q_mvar)
            q_mvar = gen.q_mvar
        else
            # Calculate reactive power based on power factor
            p_mw = gen.p_mw * gen.scaling
            if gen.cos_phi > 0 && p_mw > 0
                q_mvar = p_mw * tan(acos(gen.cos_phi))
            end
        end
        
        # Base capacity
        mbase = gen.sn_mva > 0 ? gen.sn_mva : case.baseMVA
        
        # Ramp rate parameters
        ramp_agc = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                   gen.ramp_up_rate_mw_per_min : 
                   (gen.max_p_mw - gen.min_p_mw) / 10
        ramp_10 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 10 : 
                  gen.max_p_mw - gen.min_p_mw
        ramp_30 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 30 : 
                  gen.max_p_mw - gen.min_p_mw
        
        # Fill generator data
        gen_data[idx, :] = [
            bus_idx,                               # Generator connection bus number
            gen.p_mw * gen.scaling,                # Active power output (MW)
            q_mvar,                                # Reactive power output (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power output (MVAr)
            gen.min_q_mvar,                        # Minimum reactive power output (MVAr)
            gen.vm_pu,                             # Voltage magnitude setpoint (p.u.)
            mbase,                                 # Generator base capacity (MVA)
            1.0,                                   # Generator status (1=in service, 0=out of service)
            gen.max_p_mw,                          # Maximum active power output (MW)
            gen.min_p_mw,                          # Minimum active power output (MW)
            gen.min_p_mw,                          # PQ capability curve lower end active power output (MW)
            gen.max_p_mw,                          # PQ capability curve upper end active power output (MW)
            gen.min_q_mvar,                        # Minimum reactive power output at PC1 (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power output at PC1 (MVAr)
            gen.min_q_mvar,                        # Minimum reactive power output at PC2 (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power output at PC2 (MVAr)
            ramp_agc,                              # AGC ramp rate (MW/min)
            ramp_10,                               # 10-minute reserve ramp rate (MW)
            ramp_30,                               # 30-minute reserve ramp rate (MW)
            (gen.max_q_mvar - gen.min_q_mvar)/10,  # Reactive power ramp rate (MVAr/min)
            1.0,                                   # Area participation factor
            2.0,                                   # Generator model (2=polynomial cost model)
            0.0,                                   # Startup cost (dollars)
            0.0,                                   # Shutdown cost (dollars)
            3.0,                                   # Number of polynomial cost function coefficients
            0.0                                    # Cost function parameters (to be expanded later)
        ]
        
        # If the bus is not already set as a reference node, set it as a PV node
        if jpc.busAC[bus_idx, 2] != 3  # 3 indicates REF node
            jpc.busAC[bus_idx, 2] = 2  # 2 indicates PV node
        end
    end
    
    # Process static generators (usually serving as PQ nodes, but can be PV nodes if they have voltage control capability)
    offset = n_ext + n_gen
    for (i, sgen) in enumerate(case.sgensAC)
        if !sgen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = sgen.bus
        
        # Fill generator data
        gen_data[idx, :] = [
            bus_idx,                                # Generator connection bus number
            sgen.p_mw * sgen.scaling,               # Active power output (MW)
            sgen.q_mvar * sgen.scaling,             # Reactive power output (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power output (MVAr)
            sgen.min_q_mvar,                        # Minimum reactive power output (MVAr)
            1.0,                                    # Voltage magnitude setpoint (p.u.)
            case.baseMVA,                           # Generator base capacity (MVA)
            1.0,                                    # Generator status (1=in service, 0=out of service)
            sgen.max_p_mw,                          # Maximum active power output (MW)
            sgen.min_p_mw,                          # Minimum active power output (MW)
            sgen.min_p_mw,                          # PQ capability curve lower end active power output (MW)
            sgen.max_p_mw,                          # PQ capability curve upper end active power output (MW)
            sgen.min_q_mvar,                        # Minimum reactive power output at PC1 (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power output at PC1 (MVAr)
            sgen.min_q_mvar,                        # Minimum reactive power output at PC2 (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power output at PC2 (MVAr)
            (sgen.max_p_mw - sgen.min_p_mw) / 10,   # AGC ramp rate (MW/min)
            sgen.max_p_mw - sgen.min_p_mw,          # 10-minute reserve ramp rate (MW)
            sgen.max_p_mw - sgen.min_p_mw,          # 30-minute reserve ramp rate (MW)
            (sgen.max_q_mvar - sgen.min_q_mvar)/10, # Reactive power ramp rate (MVAr/min)
            1.0,                                    # Area participation factor
            2.0,                                    # Generator model (2=polynomial cost model)
            0.0,                                    # Startup cost (dollars)
            0.0,                                    # Shutdown cost (dollars)
            3.0,                                    # Number of polynomial cost function coefficients
            0.0                                     # Cost function parameters (to be expanded later)
        ]
        
        # If the static generator is controllable and the bus is not already set as REF or PV node, it may be set as PV node
        if sgen.controllable && jpc.busAC[bus_idx, 2] == 1  # 1 indicates PQ node
            jpc.busAC[bus_idx, 2] = 2  # 2 indicates PV node
        end
    end
    
    # Remove unused rows (corresponding to out-of-service generation devices)
    active_rows = findall(x -> x > 0, gen_data[:, 8])  # Column 8 is GEN_STATUS
    gen_data = gen_data[active_rows, :]
    
    # Store generator data in JPC structure
    jpc.genAC = gen_data
    
    # Ensure at least one slack node exists
    if !any(jpc.busAC[:, 2] .== 3) && size(gen_data, 1) > 0  # 3 indicates REF node
        # If no slack node, select the first generator's bus as the slack node
        first_gen_bus = Int(gen_data[1, 1])
        jpc.busAC[first_gen_bus, 2] = 3  # 3 indicates REF node
    end
end


"""
    JPC_battery_gens_process(case::JuliaPowerCase, jpc::JPC)

Process battery storage devices and create virtual generators for them in the JPC structure.

This function creates virtual generators to represent battery storage devices in the power system model.
Each battery is assigned a virtual bus and corresponding generator parameters based on its characteristics.

The function performs the following steps:
1. Makes a deep copy of the storage devices from the input case
2. If no batteries exist, returns without modifications
3. Determines virtual node numbers based on the current busDC matrix size
4. Creates matrices for battery generators and storage information
5. For each battery:
   - Assigns a virtual bus number
   - Calculates power capacity based on battery parameters
   - Sets generator parameters (power limits, voltage setpoint, status, etc.)
   - Configures storage-specific parameters
6. Adds the battery generators to the JPC's genDC field
7. Adds the storage information to the JPC's storageetap field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing battery storage data
- `jpc::JPC`: The JPC object where the virtual generators will be added

# Returns
The updated JPC object with battery virtual generators added to genDC and storage information
added to storageetap
"""
function JPC_battery_gens_process(case::JuliaPowerCase, jpc::JPC)
    # Create virtual generators for battery virtual nodes
    batteries = deepcopy(case.storageetap)
    num_batteries = length(batteries)
    
    # If no batteries, return directly
    if num_batteries == 0
        return jpc
    end
    
    # Get current busDC size to determine virtual node numbers
    busDC_size = size(jpc.busDC, 1) - num_batteries
    
    # Create battery virtual generator matrix
    # genDC matrix typically includes the following columns:
    # [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...]
    # The specific number of columns should match your genDC structure
    num_gen_cols = size(jpc.genAC, 2)
    if num_gen_cols == 0  # If genDC is empty, set a default number of columns
        num_gen_cols = 10
    end
    
    battery_gens = zeros(num_batteries, num_gen_cols)
     # Create storage matrix
    num_storage_cols = 5  # Based on the number of columns defined by idx_ess function
    storage_matrix = zeros(num_batteries, num_storage_cols)
    
    for (i, battery) in enumerate(batteries)
        # Calculate battery virtual node number
        virtual_bus = busDC_size + i
        
        # Calculate battery power capacity (based on battery parameters)
        # This is a simplified calculation; actual calculation should be based on battery characteristics
        power_capacity = battery.package * battery.voc
        
        # Fill virtual generator matrix
        battery_gens[i, 1] = virtual_bus       # GEN_BUS: Generator connection node number
        battery_gens[i, 2] = 0.0               # PG: Initial active power output(MW), initially set to 0
        battery_gens[i, 3] = 0.0               # QG: Initial reactive power output(MVAR), usually 0 for DC systems
        
        # Set reactive power limits (usually not considered for DC systems)
        battery_gens[i, 4] = 0.0               # QMAX: Maximum reactive power output
        battery_gens[i, 5] = 0.0               # QMIN: Minimum reactive power output
        
        # Set voltage and base power
        battery_gens[i, 6] = 1.0               # VG: Voltage setpoint(p.u.)
        battery_gens[i, 7] = case.baseMVA      # MBASE: Generator base power(MVA)
        
        # Set generator status
        battery_gens[i, 8] = battery.in_service ? 1.0 : 0.0  # GEN_STATUS: Generator status
        
        # Set active power limits (charging is negative, discharging is positive)
        battery_gens[i, 9] = power_capacity    # PMAX: Maximum active power output(MW), discharging power
        battery_gens[i, 10] = -power_capacity  # PMIN: Minimum active power output(MW), charging power

         # Fill storage matrix
        storage_matrix[i, ESS_BUS] = virtual_bus               # ESS_BUS: Connected node number
        # storage_matrix[i, ESS_POWER_CAPACITY] = power_capacity # ESS_POWER_CAPACITY: Power capacity(MW)
        # storage_matrix[i, ESS_ENERGY_CAPACITY] = 0
        # storage_matrix[i, ESS_AREA] = 1                        # ESS_AREA: Area number, default is 1
        
        # If genDC has more columns, set other parameters as needed
        if num_gen_cols > 10
            # For example, set ramp rate limits, cost coefficients, etc.
            # This should be set according to your system's specific requirements
            for j in 11:num_gen_cols
                battery_gens[i, j] = 0.0  # Default to 0
            end
        end
    end
    
    # Add battery virtual generators to genDC
    if isempty(jpc.genDC)
        jpc.genDC = battery_gens
    else
        jpc.genDC = [jpc.genDC; battery_gens]
    end
    
    # Add storage device information to storage
    if !isdefined(jpc, :storageetap) || isempty(jpc.storageetap)
        jpc.storageetap = storage_matrix
    else
        jpc.storageetap = [jpc.storageetap; storage_matrix]
    end

    return jpc
end


"""
    JPC_loads_process(case::JuliaPowerCase, jpc::JPC)

Process AC load data from a JuliaPowerCase and convert it to JPC format.

This function extracts information about AC loads from the input case, applies scaling factors,
and stores the processed data in the JPC structure. It also updates the power demand values
in the corresponding AC bus data.

For each AC load, the function:
1. Records its basic properties (index, bus connection)
2. Calculates active and reactive power demand with scaling applied
3. Stores load model composition (ZIP model percentages)
4. Accumulates loads connected to the same bus

The function performs the following steps:
1. Filters the input case to identify in-service AC loads
2. If no in-service loads exist, returns without modifications
3. Creates a matrix to store AC load data with 8 columns
4. For each AC load:
   - Calculates actual power demand with scaling applied
   - Populates the matrix with all required parameters
   - Accumulates loads by bus connection
5. Adds the processed AC load data to the JPC object's loadAC field
6. Updates the busAC matrix with the accumulated load values for each bus

# Arguments
- `case::JuliaPowerCase`: The original power system case containing AC load data
- `jpc::JPC`: The JPC object where the processed AC load data will be stored

# Returns
None. The JPC object is modified in-place.
"""
function JPC_loads_process(case::JuliaPowerCase, jpc::JPC)
    # Process load data, convert to JPC format and update busAC's PD and QD
    
    # Filter out in-service loads
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # If no in-service loads, return directly
    if isempty(in_service_loads)
        return
    end
    
    # Create an empty matrix, with rows equal to number of loads and 8 columns
    num_loads = length(in_service_loads)
    load_matrix = zeros(num_loads, 8)
    
    # Create a dictionary to accumulate loads connected to the same bus
    bus_load_sum = Dict{Int, Vector{Float64}}()
    
    for (i, load) in enumerate(in_service_loads)
        # Calculate actual active and reactive loads (considering scaling factor)
        # if mode =="1_ph_pf"
            actual_p_mw = load.p_mw * load.scaling
            actual_q_mvar = load.q_mvar * load.scaling
        # else
        #     actual_p_mw = load.p_mw * load.scaling / 3.0
        #     actual_q_mvar = load.q_mvar * load.scaling / 3.0
        # end
        
        # Fill each row of the load matrix
        load_matrix[i, :] = [
            i,                      # Load number
            load.bus,               # Bus number where load is connected
            1.0,                    # Load status (1=in service)
            actual_p_mw,            # Active load (MW)
            actual_q_mvar,          # Reactive load (MVAr)
            load.const_z_percent/100,  # Constant impedance load percentage
            load.const_i_percent/100,  # Constant current load percentage
            load.const_p_percent/100   # Constant power load percentage
        ]
        
        # Accumulate loads connected to the same bus
        bus_idx = load.bus
        if haskey(bus_load_sum, bus_idx)
            bus_load_sum[bus_idx][1] += actual_p_mw
            bus_load_sum[bus_idx][2] += actual_q_mvar
        else
            bus_load_sum[bus_idx] = [actual_p_mw, actual_q_mvar]
        end
    end
    
    # Store load data in JPC structure
    jpc.loadAC = load_matrix
    
    # Update PD and QD fields in busAC matrix
    for (bus_idx, load_values) in bus_load_sum
        # Find the corresponding bus row
        bus_row = findfirst(x -> x == bus_idx, jpc.busAC[:, 1])
        
        if !isnothing(bus_row)
            # Update PD (column 3) and QD (column 4)
            jpc.busAC[bus_row, PD] = load_values[1]  # PD - Active load (MW)
            jpc.busAC[bus_row, QD] = load_values[2]  # QD - Reactive load (MVAr)
        end
    end
end


"""
    JPC_dcloads_process(case::JuliaPowerCase, jpc::JPC)

Process DC load data from a JuliaPowerCase and convert it to JPC format.

This function extracts information about DC loads from the input case, applies scaling factors,
and stores the processed data in the JPC structure. It also updates the power demand values
in the corresponding DC bus data.

For each DC load, the function:
1. Records its basic properties (index, bus connection)
2. Calculates active power demand with scaling applied
3. Stores load model composition (ZIP model percentages)
4. Updates the connected DC bus with the additional power demand

The function performs the following steps:
1. Filters the input case to identify in-service DC loads
2. If no in-service loads exist, returns without modifications
3. Creates a matrix to store DC load data with 8 columns
4. For each DC load:
   - Populates the matrix with all required parameters
   - Updates the corresponding bus data with the load's power demand
5. Adds the processed DC load data to the JPC object's loadDC field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing DC load data
- `jpc::JPC`: The JPC object where the processed DC load data will be stored

# Returns
The updated JPC object with DC load data in the loadDC field and updated busDC power demands
"""
function JPC_dcloads_process(case::JuliaPowerCase, jpc::JPC)
    # Process DC load data, convert to JPC format and update busDC's PD and QD
    
    # Filter out in-service DC loads
    in_service_dcloads = filter(dcload -> dcload.in_service == true, case.loadsDC)
    
    # If no in-service DC loads, return directly
    if isempty(in_service_dcloads)
        return
    end
    
    # Create an empty matrix, with rows equal to number of DC loads and 8 columns
    num_dcloads = length(in_service_dcloads)
    dcload_matrix = zeros(num_dcloads, 8)
    
    for (i, dcload) in enumerate(in_service_dcloads)
        # Fill each row of the DC load matrix
        dcload_matrix[i, :] = [
            dcload.index,              # DC load number
            dcload.bus,     # Bus number where DC load is connected
            1.0,            # DC load status (1=in service)
            dcload.p_mw * dcload.scaling,  # Active load (MW)
            0.0,            # Reactive load (MVAr)
            dcload.const_z_percent/100,    # Constant impedance percentage (default 0)
            dcload.const_i_percent/100,    # Constant current percentage (default 0)
            dcload.const_p_percent/100     # Constant power percentage (default 0)
        ]
        
        # Update PD and QD fields in busDC matrix
        bus_row = findfirst(x -> x == dcload.bus, jpc.busDC[:, 1])
        
        if !isnothing(bus_row)
            jpc.busDC[bus_row, PD] += dcload_matrix[i, 4]  # PD - Active load (MW)
            jpc.busDC[bus_row, QD] += dcload_matrix[i, 5]  # QD - Reactive load (MVAr)
        end
    end
    
    # Store DC load data in JPC structure
    jpc.loadDC = dcload_matrix

    return jpc
end


"""
    JPC_pv_process(case::JuliaPowerCase, jpc::JPC)

Process photovoltaic (PV) array data from a JuliaPowerCase and convert it to JPC format.

This function extracts information about PV arrays from the input case, performs necessary
calculations to determine their electrical characteristics based on environmental conditions,
and stores the processed data in the JPC structure.

For each PV array, the function:
1. Calculates key electrical parameters adjusted for environmental conditions:
   - Open-circuit voltage (Voc) adjusted for temperature and series connections
   - Maximum power point voltage (Vmpp) adjusted for series-connected panels
   - Short-circuit current (Isc) adjusted for temperature, irradiance, and parallel connections
   - Maximum power point current (Impp) adjusted for irradiance and parallel-connected panels
2. Stores all parameters in a standardized matrix format

The function performs the following steps:
1. Filters the input case to identify in-service PV arrays
2. If no in-service arrays exist, returns without modifications
3. Creates a matrix to store PV array data with 9 columns
4. For each PV array:
   - Calculates electrical parameters based on temperature and irradiance conditions
   - Populates the matrix with all required parameters
5. Adds the processed PV array data to the JPC object's pv field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing PV array data
- `jpc::JPC`: The JPC object where the processed PV array data will be stored

# Returns
The updated JPC object with PV array data in the pv field
"""
function JPC_pv_process(case::JuliaPowerCase, jpc::JPC)
    # Process PV array data, convert to JPC format and update busAC's PD and QD
    
    # Filter out in-service PV arrays
    in_service_pvs = filter(pv -> pv.in_service == true, case.pvarray)
    
    # If no in-service PV arrays, return directly
    if isempty(in_service_pvs)
        return
    end
    
    # Create an empty matrix, with rows equal to number of PV arrays and 9 columns
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
        # Fill each row of the PV array matrix
        pv_matrix[i, :] = [
            i,              # PV array number
            pv.bus,         # Bus number where PV array is connected
            Voc,            # Open circuit voltage (V)
            Vmpp,           # Maximum power point voltage (V)
            Isc,            # Short circuit current (A)
            Impp,           # Maximum power point current (A)
            pv.irradiance,  # Solar irradiance (W/m²)
            1.0,            # area
            1.0,            # PV array status (1=in service)
        ]
        
        # # Update PD and QD fields in busAC matrix
        # bus_row = findfirst(x -> x == pv.bus, jpc.busAC[:, 1])
        
        # if !isnothing(bus_row)
        #     jpc.busAC[bus_row, PD] += pv_matrix[i, 4]  # PD - Active load (MW)
        #     jpc.busAC[bus_row, QD] += pv_matrix[i, 5]  # QD - Reactive load (MVAr)
        # end
    end
    
    # Store PV array data in JPC structure
    jpc.pv = pv_matrix

    return jpc
    
end


"""
    JPC_ac_pv_system_process(case::JuliaPowerCase, jpc::JPC)

Processes AC-side photovoltaic (PV) systems from a JuliaPowerCase and converts them to JPC format.

This function extracts information about AC-connected PV systems from the input case,
performs necessary calculations to determine their electrical characteristics, and stores
the processed data in the JPC structure.

For each AC-side PV system, the function:
1. Calculates operating voltages and currents based on panel configuration and conditions:
   - Maximum power point voltage (Vmpp) adjusted for series-connected panels
   - Open-circuit voltage (Voc) adjusted for temperature and series connections
   - Short-circuit current (Isc) adjusted for temperature, irradiance, and parallel connections
   - Maximum power point current (Impp) adjusted for parallel-connected panels
2. Calculates maximum power output accounting for system losses
3. Determines the control mode (Voltage Control or MVar Control)
4. Stores all parameters in a standardized matrix format

The function performs the following steps:
1. Filters the input case to identify in-service AC-side PV systems
2. If no in-service systems exist, returns the JPC object unchanged
3. Creates a matrix to store PV system data
4. For each PV system:
   - Calculates electrical parameters based on environmental conditions
   - Sets control mode flags
   - Populates the matrix with all required parameters
5. Adds the processed PV system data to the JPC object's pv_acsystem field

# Arguments
- `case::JuliaPowerCase`: The original power system case containing AC-side PV system data
- `jpc::JPC`: The JPC object where the processed PV system data will be stored

# Returns
The updated JPC object with AC-side PV system data in the pv_acsystem field
"""
function JPC_ac_pv_system_process(case::JuliaPowerCase, jpc::JPC)
    # process AC side PV systems, converting to JPC format and updating busAC
    
    # filter out in-service AC side PV systems
    in_service_ac_pvs = filter(ac_pv -> ac_pv.in_service == true, case.ACPVSystems)
    
    # if no in-service AC side PV systems, return
    if isempty(in_service_ac_pvs)
        return jpc
    end
    
    # create an empty matrix, number of rows = number of AC side PV systems, columns = 15
    num_ac_pvs = length(in_service_ac_pvs)
    ac_pv_matrix = zeros(num_ac_pvs, 15)
    for (i, ac_pv) in enumerate(in_service_ac_pvs)
        Vmpp = ac_pv.vmpp * ac_pv.numpanelseries
        Voc = (ac_pv.voc + (ac_pv.temperature - 25) * ac_pv.β_voc) * ac_pv.numpanelseries
        Isc = (ac_pv.isc + (ac_pv.temperature - 25) * ac_pv.α_isc) * (ac_pv.irradiance / 1000.0) * ac_pv.numpanelparallel
        Impp = ac_pv.impp * ac_pv.numpanelparallel

        p_max = Vmpp * Impp / 1000000.0 * (1-ac_pv.loss_percent)# maximum active power output (MW)

        if ac_pv.control_mode == "Voltage Control"
            mode = 1
        else
            mode = 0 
        end
        
        # fill the AC side PV system matrix for each row
        ac_pv_matrix[i, :] = [
            i,                # ac side PV system index
            ac_pv.bus,            # connected bus index
            Voc,               # ac side PV system rated voltage(V)
            Vmpp,              # base voltage(V)
            Isc,               # short-circuit current(A)
            Impp,              # base current(A)
            ac_pv.irradiance,  # irradiance(W/m²)
            ac_pv.loss_percent,# loss percentage
            mode,              # control mode (0=MVar Control, 1=Voltage Control)
            ac_pv.p_mw,           # active power output(MW)
            ac_pv.q_mvar,         # reactive power output(MVAr)
            ac_pv.max_q_mvar,     # reactive power upper limit(MVAr)
            ac_pv.min_q_mvar,     # reactive power lower limit(MVAr)
            1,            # area (default 1)
            ac_pv.in_service ? 1.0 : 0.0  # pv system status (1=operational, 0=out of service)
        ]
    end
    # store the AC side PV system data in the JPC structure
    jpc.pv_acsystem = ac_pv_matrix
    
    return jpc
end


"""
    JPC_inverters_process(case::JuliaPowerCase, jpc::JPC)

Process converters/inverters from a JuliaPowerCase and integrate them into the JPC format.

This function handles the conversion of power electronic converters (inverters) from the 
JuliaPowerCase format to the JPC format, updating the relevant bus types and power flow parameters
according to each inverter's control mode.

For each in-service inverter, the function:
1. Determines the control mode and sets appropriate flags
2. Calculates AC and DC side power flows considering efficiency losses
3. Updates bus types for both AC and DC sides based on control mode
4. Creates generator entries when needed for voltage control
5. Updates load records for power injection/consumption
6. Handles ZIP load model parameters through weighted averaging

The function supports the following control modes:
- δs_Us: Slack bus on AC side, P node on DC side
- Ps_Qs: PQ node on AC side, P node on DC side
- Ps_Us: PV node on AC side, P node on DC side
- Udc_Qs: PQ node on AC side, slack node on DC side
- Udc_Us: PV node on AC side, slack node on DC side
- Droop_Udc_Qs: PQ node on AC side with droop control, slack node on DC side
- Droop_Udc_Us: PV node on AC side with droop control, slack node on DC side

# Arguments
- `case::JuliaPowerCase`: The original power system case containing converter data
- `jpc::JPC`: The JPC object where the processed converter data will be stored

# Returns
The updated JPC object with converter data integrated into the appropriate fields
"""
function JPC_inverters_process(case::JuliaPowerCase, jpc::JPC)
    # process inverters, converting to JPC format and updating busAC
    
    # filter out in-service inverters
    in_service_inverters = filter(inverter -> inverter.in_service == true, case.converters)
    
    # if no in-service inverters, return
    if isempty(in_service_inverters)
        return jpc
    end
    
    # obtain the number of external buses and existing generators
    nld_ac = size(jpc.loadAC, 1)  # ac side load count
    nld_dc = size(jpc.loadDC, 1)  # dc side load count
    
    # create an empty matrix for inverters, number of rows = number of in-service inverters, columns = 18
    # use the same number of columns as jpc.loadAC and jpc.loadDC
    num_cols_ac = size(jpc.loadAC, 2)
    num_cols_dc = size(jpc.loadDC, 2)
    
    # calculate the maximum number of new loads
    max_new_loads = length(in_service_inverters)
    new_loads_ac = zeros(0, num_cols_ac)  # create an empty matrix, rows = 0, columns = same as loadAC
    new_loads_dc = zeros(0, num_cols_dc)  # create an empty matrix, rows = 0, columns = same as loadDC

    # create an empty matrix for converters, number of rows = 0, columns = 18
    converters = zeros(0, 18)
    
    # follow the same column order as jpc.loadAC and jpc.loadDC
    new_ac_load_count = 0
    new_dc_load_count = 0

    
    
   for (i, inverter) in enumerate(in_service_inverters)
        # add a new row for each inverter
        converter = zeros(1, 18)  # create a new row with 18 columns
        
        # control mode
        mode = inverter.control_mode
        if mode == "δs_Us"
            converter[1, CONV_MODE] = 1.0  # δs_Us mode inverter
        elseif mode == "Ps_Qs"
            converter[1, CONV_MODE] = 0.0  # Ps_Qs mode inverter
        elseif mode == "Ps_Us"
            converter[1, CONV_MODE] = 3.0  # Ps_Us mode inverter
        elseif mode == "Udc_Qs"
            converter[1, CONV_MODE] = 4.0  # Udc_Qs mode inverter
        elseif mode == "Udc_Us"
            converter[1, CONV_MODE] = 5.0  # Udc_Us mode inverter
        elseif mode == "Droop_Udc_Qs"
            converter[1, CONV_MODE] = 6.0  # Droop_Udc_Qs mode inverter
        elseif mode == "Droop_Udc_Us"
            converter[1, CONV_MODE] = 7.0  # Droop_Udc_Us mode inverter
        else
            @warn "Control mode $mode of inverter $i is unknown or not supported, Ps_Qs is enabled by default"
            converter[1, CONV_MODE] = 0.0  # setting default to Ps_Qs mode
        end

        # calculate AC and DC side power
        p_ac = -inverter.p_mw 
        q_ac = -inverter.q_mvar 
        
        # calculate DC side power based on AC side power and efficiency
        efficiency = 1.0 - inverter.loss_percent   # transform loss percentage to efficiency
        
        if p_ac <= 0  # output power on AC side, input power on DC side
            p_dc = -p_ac / efficiency  # negative value, indicating input power on DC side
        else  # output power on DC side, input power on AC side
            p_dc = -p_ac * efficiency  # positive value, indicating output power on DC side
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
        
        # obtain the bus indices for AC and DC sides
        ac_bus_row = findfirst(x -> x == inverter.bus_ac, jpc.busAC[:, 1])
        dc_bus_row = findfirst(x -> x == inverter.bus_dc, jpc.busDC[:, 1])
        
        # decide how to handle the inverter based on its control mode
        if mode =="δs_Us"
            # δs_Us mode : without any operations
        elseif mode == "Ps_Qs"
            # Ps_Qs mode : update the AC side bus PD and QD
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, PD] += p_ac  # PD - ac side active power (MW)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - reactive power (MVAr)
            end
            
            if !isnothing(dc_bus_row)
                jpc.busDC[dc_bus_row, PD] += p_dc  # PD - active power (MW)
            end
            
            # process AC side loads
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # if no loads connected to the same AC bus, create a new load record
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # create a new row
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # laod index
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # node index
                new_load_ac[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_ac[1, LOAD_PD] = p_ac                        # active power (MW)
                new_load_ac[1, LOAD_QD] = q_ac                        # reactive power (MVAr)
                # default to ZIP load with no reactive power
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # if loads connected to the same AC bus, update these loads
                for idx in existing_load_indices_ac
                    # obtain the original load's power and ZIP percentages
                    orig_p = jpc.loadAC[idx, LOAD_PD]
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    orig_z_percent = jpc.loadAC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadAC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadAC[idx, LOADP_PERCENT]
                    
                    # calculate the new total power
                    new_p = orig_p + p_ac
                    new_q = orig_q + q_ac
                    
                    # re calculate the ZIP percentages (weighted average)
                    # avoid division by zero
                    if new_p != 0
                        # weight of the original load
                        w_orig = abs(orig_p) / abs(new_p)
                        # weight of the inverter load (default to constant power)
                        w_inv = abs(p_ac) / abs(new_p)
                        
                        # re calculate the new ZIP percentages
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # assure the percentages sum to 1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # if the sum is 0, set to default values
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # if the new total power is 0, keep the original ZIP percentages
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # update the load matrix
                    jpc.loadAC[idx, LOAD_PD] = new_p
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    jpc.loadAC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadAC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadAC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
            
            # process DC side loads
            existing_load_indices_dc = findall(x -> x == inverter.bus_dc, jpc.loadDC[:, 2])
            
            if isempty(existing_load_indices_dc)
                # if no loads connected to the same DC bus, create a new load record
                new_dc_load_count += 1
                new_load_dc = zeros(1, num_cols_dc)  # create a new row
                new_load_dc[1, LOAD_I] = nld_dc + new_dc_load_count  # load index
                new_load_dc[1, LOAD_CND] = inverter.bus_dc             # node index
                new_load_dc[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_dc[1, LOAD_PD] = p_dc                        # active power (MW)
                new_load_dc[1, LOAD_QD] = 0.0                         # reactive power (MVAr)
                # reactive power is not considered in DC side
                # constant impedance, current, and power percentages
                new_load_dc[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_dc[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_dc[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_dc = vcat(new_loads_dc, new_load_dc)
            else
                # if loads connected to the same DC bus, update these loads
                for idx in existing_load_indices_dc
                    # obtain the original load's power and ZIP percentages
                    orig_p = jpc.loadDC[idx, 4]
                    orig_z_percent = jpc.loadDC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadDC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadDC[idx, LOADP_PERCENT]
                    
                    # calculate the new total power
                    new_p = orig_p + p_dc
                    
                    # update the ZIP percentages (weighted average)
                    # avoid division by zero
                    if new_p != 0
                        # weight of the original load
                        w_orig = abs(orig_p) / abs(new_p)
                        # weight of the inverter load (default to constant power)
                        w_inv = abs(p_dc) / abs(new_p)
                        
                        # calculate the new ZIP percentages
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # assure the percentages sum to 1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # if the sum is 0, set to default values
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # if the new total power is 0, keep the original ZIP percentages
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # update the load matrix
                    jpc.loadDC[idx, LOAD_PD] = new_p
                    jpc.loadDC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadDC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadDC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
        elseif mode == "Ps_Us"
            # Ps_Us mode : without any operations
           
        elseif mode == "Udc_Qs"
            # Udc_Qs mode : only update the AC side bus QD
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - reactive power (MVAr)
                # do not modify PD
            end
            
            # process AC side loads - only modify reactive power
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # if no loads connected to the same AC bus, create a new load record
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # create a new row
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # load index
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # bus index
                new_load_ac[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_ac[1, LOAD_PD] = 0.0                         # active power (MW) - do not modify
                new_load_ac[1, LOAD_QD] = q_ac                        # reactive power (MVAr)
                # default to ZIP load with no active power
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # if loads connected to the same AC bus, update these loads
                for idx in existing_load_indices_ac
                    # obtain the original load's power
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    
                    # calculate the new total reactive power
                    new_q = orig_q + q_ac
                    
                    # update the load matrix - only modify reactive power
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    # do not modify active power or ZIP percentages
                end
            end
        elseif mode == "Udc_Us"
            # Udc_Us mode : without any operations
        elseif mode == "Droop_Udc_Qs"
             # Udc_Qs mode : update the AC side bus QD and DC side bus PD
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - reactive power (MVAr)
                # do not modify PD
            end
            
            # process AC side loads - only modify reactive power
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # if no loads connected to the same AC bus, create a new load record
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # load index
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # bus index
                new_load_ac[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_ac[1, LOAD_PD] = 0.0                         # active power (MW) - do not modify
                new_load_ac[1, LOAD_QD] = q_ac                        # reactive power (MVAr)
                # default to ZIP load with no active power
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # if loads connected to the same AC bus, update these loads
                for idx in existing_load_indices_ac
                    # obtain the original load's power
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    
                    # calculate the new total reactive power
                    new_q = orig_q + q_ac
                    
                    # update the load matrix - only modify reactive power
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    # do not modify active power or ZIP percentages
                end
            end
        elseif mode == "Droop_Udc_Us"
            # Droop_Udc_Us mode : update the AC side bus PD and QD, and DC side bus PD
        else
            # unknown control mode, default to Ps_Qs
            if !isnothing(ac_bus_row)
                jpc.busAC[ac_bus_row, PD] += p_ac  # PD - active power (MW)
                jpc.busAC[ac_bus_row, QD] += q_ac  # QD - reactive power (MVAr)
            end
            
            if !isnothing(dc_bus_row)
                jpc.busDC[dc_bus_row, PD] += p_dc  # PD - active power (MW)
            end
            
            # process AC side loads
            existing_load_indices_ac = findall(x -> x == inverter.bus_ac, jpc.loadAC[:, 2])
            
            if isempty(existing_load_indices_ac)
                # if no loads connected to the same AC bus, create a new load record
                new_ac_load_count += 1
                new_load_ac = zeros(1, num_cols_ac)  # create a new row
                new_load_ac[1, LOAD_I] = nld_ac + new_ac_load_count  # load index
                new_load_ac[1, LOAD_CND] = inverter.bus_ac             # bus index
                new_load_ac[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_ac[1, LOAD_PD] = p_ac                        # active power (MW)
                new_load_ac[1, LOAD_QD] = q_ac                        # reactive power (MVAr)
                # default to ZIP load with no reactive power
                new_load_ac[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_ac[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_ac[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_ac = vcat(new_loads_ac, new_load_ac)
            else
                # if loads connected to the same AC bus, update these loads
                for idx in existing_load_indices_ac
                    # obtain the original load's power and ZIP percentages
                    orig_p = jpc.loadAC[idx, LOAD_PD]
                    orig_q = jpc.loadAC[idx, LOAD_QD]
                    orig_z_percent = jpc.loadAC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadAC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadAC[idx, LOADP_PERCENT]
                    
                    # calculate the new total power
                    new_p = orig_p + p_ac
                    new_q = orig_q + q_ac
                    
                    # re calculate the ZIP percentages (weighted average)
                    # avoid division by zero
                    if new_p != 0
                        # weight of the original load
                        w_orig = abs(orig_p) / abs(new_p)
                        # weight of the inverter load (default to constant power)
                        w_inv = abs(p_ac) / abs(new_p)
                        
                        # calculate the new ZIP percentages
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # assure the percentages sum to 1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # if the sum is 0, set to default values
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # if the new total power is 0, keep the original ZIP percentages
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # update the load matrix
                    jpc.loadAC[idx, LOAD_PD] = new_p
                    jpc.loadAC[idx, LOAD_QD] = new_q
                    jpc.loadAC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadAC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadAC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
            
            # process DC side loads
            existing_load_indices_dc = findall(x -> x == inverter.bus_dc, jpc.loadDC[:, 2])
            
            if isempty(existing_load_indices_dc)
                # if no loads connected to the same DC bus, create a new load record
                new_dc_load_count += 1
                new_load_dc = zeros(1, num_cols_dc)  # create a new row
                new_load_dc[1, LOAD_I] = nld_dc + new_dc_load_count  # laod index
                new_load_dc[1, LOAD_CND] = inverter.bus_dc             # bus index
                new_load_dc[1, LOAD_STATUS] = 1.0                         # status(1=on)
                new_load_dc[1, LOAD_PD] = p_dc                        # active power (MW)
                new_load_dc[1, LOAD_QD] = 0.0                         # reactive power (MVAr)
                # dc side does not consider reactive power
                # dc side constant impedance, current, and power percentages
                new_load_dc[1, LOADZ_PERCENT] = 0.0                         # constant impedance percentage
                new_load_dc[1, LOADI_PERCENT] = 0.0                         # constant current percentage
                new_load_dc[1, LOADP_PERCENT] = 1.0                         # constant power percentage
                
                # add to new loads matrix
                new_loads_dc = vcat(new_loads_dc, new_load_dc)
            else
                # if loads connected to the same DC bus, update these loads
                for idx in existing_load_indices_dc
                    # obtain the original load's power and ZIP percentages
                    orig_p = jpc.loadDC[idx, 4]
                    orig_z_percent = jpc.loadDC[idx, LOADZ_PERCENT]
                    orig_i_percent = jpc.loadDC[idx, LOADI_PERCENT]
                    orig_p_percent = jpc.loadDC[idx, LOADP_PERCENT]
                    
                    # calculate the new total power
                    new_p = orig_p + p_dc
                    
                    # re calculate the ZIP percentages (weighted average)
                    # avoid division by zero
                    if new_p != 0
                        # 原始负荷的权重
                        w_orig = abs(orig_p) / abs(new_p)
                        # weight of the inverter load (default to constant power)
                        w_inv = abs(p_dc) / abs(new_p)
                        
                        # calculate the new ZIP percentages
                        new_z_percent = orig_z_percent * w_orig + 0.0 * w_inv
                        new_i_percent = orig_i_percent * w_orig + 0.0 * w_inv
                        new_p_percent = orig_p_percent * w_orig + 1.0 * w_inv
                        
                        # assure the percentages sum to 1
                        sum_percent = new_z_percent + new_i_percent + new_p_percent
                        if sum_percent != 0
                            new_z_percent /= sum_percent
                            new_i_percent /= sum_percent
                            new_p_percent /= sum_percent
                        else
                            # if the sum is 0, set to default values
                            new_z_percent = 0.0
                            new_i_percent = 0.0
                            new_p_percent = 1.0
                        end
                    else
                        # if the new total power is 0, keep the original ZIP percentages
                        new_z_percent = orig_z_percent
                        new_i_percent = orig_i_percent
                        new_p_percent = orig_p_percent
                    end
                    
                    # update the load matrix
                    jpc.loadDC[idx, LOAD_PD] = new_p
                    jpc.loadDC[idx, LOADZ_PERCENT] = new_z_percent
                    jpc.loadDC[idx, LOADI_PERCENT] = new_i_percent
                    jpc.loadDC[idx, LOADP_PERCENT] = new_p_percent
                end
            end
        end

        # Process JPC according to inverter control mode
        if mode == "δs_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 3.0  # Set as slack bus
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # Set as P node
            inverter_gens_ac = zeros(1, 32)  # Assume generator info has 26 columns
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: Bus number generator is connected to
            inverter_gens_ac[1, 2] = -p_ac  # PG: Initial active power output (MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: Initial reactive power output (MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens_ac[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: Voltage setpoint (p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens_ac[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # Add to genAC
        elseif mode == "Ps_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # Set as PQ node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # Set as P node
        elseif mode == "Ps_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # Set as PV node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 1.0  # Set as P node

            inverter_gens_ac = zeros(1, 26)  # Assume generator info has 32 columns
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: Bus number generator is connected to
            inverter_gens_ac[1, 2] = -p_ac  # PG: Initial active power output (MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: Initial reactive power output (MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens_ac[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: Voltage setpoint (p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens_ac[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # Add to genAC
        elseif mode == "Udc_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # Set as PQ node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # Set as slack node
            # Create inverter generator information
            inverter_gens = zeros(1, 32)  # Assume generator info has 32 columns
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: Bus number generator is connected to
            inverter_gens[1, 2] = p_dc  # PG: Initial active power output (MW)
            inverter_gens[1, 3] = 0.0  # QG: Initial reactive power output (MVAR), usually 0 for DC systems
            inverter_gens[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens[1, 6] = inverter.vm_dc_pu  # VG: Voltage setpoint (p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # Add to genDC
        elseif mode == "Udc_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # Set as PV node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # Set as slack node
            # Create inverter generator information
            inverter_gens = zeros(1, 32)  # Assume generator info has 32 columns
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: Bus number generator is connected to
            inverter_gens[1, 2] = -p_dc  # PG: Initial active power output (MW)
            inverter_gens[1, 3] = 0.0  # QG: Initial reactive power output (MVAR), usually 0 for DC systems
            inverter_gens[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens[1, 6] = inverter.vm_dc_pu  # VG: Voltage setpoint (p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # Add to genDC

            inverter_gens_ac = zeros(1, 26)  # Assume generator info has 32 columns
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: Bus number generator is connected to
            inverter_gens_ac[1, 2] = -p_ac  # PG: Initial active power output (MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: Initial reactive power output (MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens_ac[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens_ac[1, 6] = inverter.vm_ac_pu  # VG: Voltage setpoint (p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens_ac[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # Add to genAC
        elseif mode == "Droop_Udc_Qs"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 1.0  # Set as PQ node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # Set as slack node
            # Create inverter generator information
            inverter_gens = zeros(1, 32)  # Assume generator info has 32 columns
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: Bus number generator is connected to
            inverter_gens[1, 2] = p_dc  # PG: Initial active power output (MW)
            inverter_gens[1, 3] = 0.0  # QG: Initial reactive power output (MVAR), usually 0 for DC systems
            inverter_gens[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens[1, 6] = 1.0  # VG: Voltage setpoint (p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # Add to genDC
        elseif mode == "Droop_Udc_Us"
            jpc.busAC[ac_bus_row, BUS_TYPE] = 2.0  # Set as PV node
            jpc.busDC[dc_bus_row, BUS_TYPE] = 2.0  # Set as slack node
            # Create inverter generator information
            inverter_gens = zeros(1, 32)  # Assume generator info has 32 columns
            inverter_gens[1, 1] = inverter.bus_dc  # GEN_BUS: Bus number generator is connected to
            inverter_gens[1, 2] = -p_dc  # PG: Initial active power output (MW)
            inverter_gens[1, 3] = 0.0  # QG: Initial reactive power output (MVAR), usually 0 for DC systems
            inverter_gens[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens[1, 6] = 1.0  # VG: Voltage setpoint (p.u.)
            inverter_gens[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genDC = vcat(jpc.genDC, inverter_gens)  # Add to genDC

            inverter_gens_ac = zeros(1, 26)  # Assume generator info has 32 columns
            inverter_gens_ac[1, 1] = inverter.bus_ac  # GEN_BUS: Bus number generator is connected to
            inverter_gens_ac[1, 2] = -p_ac  # PG: Initial active power output (MW)
            inverter_gens_ac[1, 3] = -q_ac  # QG: Initial reactive power output (MVAR)
            inverter_gens_ac[1, 4] = 0.0  # QMAX: Maximum reactive power output
            inverter_gens_ac[1, 5] = 0.0  # QMIN: Minimum reactive power output
            inverter_gens_ac[1, 6] = 1.0  # VG: Voltage setpoint (p.u.)
            inverter_gens_ac[1, 7] = jpc.baseMVA  # MBASE: Generator base power (MVA)
            inverter_gens_ac[1, 8] = 1.0  # GEN_STATUS: Generator status (1=in service)
            inverter_gens_ac[1, 9] = 0.0  # PMAX: Maximum active power output (MW), discharge power
            inverter_gens_ac[1, 10] = 0.0  # PMIN: Minimum active power output (MW), charge power
            jpc.genAC = vcat(jpc.genAC, inverter_gens_ac)  # Add to genAC
        end

    end


    
    # add new loads to the JPC structure
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
