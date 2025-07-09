"""
    renumber_hybrid_system(jpc)

Renumber buses and branches in a hybrid AC-DC power system to create a more organized numbering scheme.

# Arguments
- `jpc`: A structure containing the hybrid power system data with fields:
  - `busAC`: Matrix of AC bus data
  - `busDC`: Matrix of DC bus data
  - `branchAC`: Matrix of AC branch data
  - `branchDC`: Matrix of DC branch data
  - `converter`: Matrix of AC-DC converter data
  - Optional fields: `loadAC`, `loadDC`, `genAC`, `genDC`, `pv`, `pv_acsystem`, `storage`

# Returns
- `new_jpc`: A structure with renumbered system data
- `ac_node_mapping`: Dictionary mapping old AC bus numbers to new numbers
- `dc_node_mapping`: Dictionary mapping old DC bus numbers to new numbers

# Description
This function performs a breadth-first search (BFS) based renumbering of buses in a hybrid AC-DC power system.
Starting from the AC slack bus, it assigns new sequential numbers to buses while preserving the system topology.
The function handles both AC and DC subsystems and their interconnections through converters.

The renumbering process:
1. Creates a unified graph representation of the hybrid system
2. Identifies the AC slack bus as the starting point
3. Performs BFS traversal to assign new sequential numbers
4. Updates all component references to use the new numbering scheme
5. Sorts components based on the new numbering for better organization

The function also ensures that for branches and converters, the bus with the smaller number
is always listed first, which helps standardize the system representation.

All system components (buses, branches, converters, loads, generators, PV arrays, and storage devices)
are updated to use the new numbering scheme, and the results are sorted for better organization.
"""

function renumber_hybrid_system(jpc)
    # Extract necessary data
    busAC = jpc.busAC
    busDC = jpc.busDC
    branchAC = jpc.branchAC
    branchDC = jpc.branchDC
    converter = jpc.converter
    
    # Extract load and generator data (if exists)
    loadAC = jpc.loadAC
    loadDC = jpc.loadDC
    genAC = jpc.genAC
    genDC = jpc.genDC
    pvarray = jpc.pv
    
    # Extract AC PV system data (if exists)
    pv_acsystem = hasfield(typeof(jpc), :pv_acsystem) ? jpc.pv_acsystem : nothing
    
    # Extract storage data (if exists)
    storage = hasfield(typeof(jpc), :storage) ? jpc.storage : nothing

    # To distinguish between AC and DC nodes, we create unique identifiers
    # Using tuple (bus_number, is_ac) to represent nodes
    # is_ac is true for AC nodes, false for DC nodes
    ac_buses = [(bus, true) for bus in busAC[:, 1]]
    dc_buses = [(bus, false) for bus in busDC[:, 1]]
    all_buses = vcat(ac_buses, dc_buses)
    
    # Create node mapping (old number -> new number)
    # Key is (bus_number, is_ac), value is the new number
    old_to_new = Dict{Tuple{Float64, Bool}, Int}()
    
    # Create adjacency list to represent the graph
    graph = Dict{Tuple{Float64, Bool}, Vector{Tuple{Float64, Bool}}}()
    
    # Initialize graph
    for bus in all_buses
        graph[bus] = Tuple{Float64, Bool}[]
    end
    
    # Add AC branch connections
    for i in 1:size(branchAC, 1)
        from_bus = (branchAC[i, 1], true)  # AC node
        to_bus = (branchAC[i, 2], true)    # AC node
        push!(graph[from_bus], to_bus)
        push!(graph[to_bus], from_bus)
    end
    
    # Add DC branch connections
    for i in 1:size(branchDC, 1)
        from_bus = (branchDC[i, 1], false)  # DC node
        to_bus = (branchDC[i, 2], false)    # DC node
        push!(graph[from_bus], to_bus)
        push!(graph[to_bus], from_bus)
    end
    
    # Add converter connections (connecting AC and DC nodes)
    for i in 1:size(converter, 1)
        ac_bus = (converter[i, 1], true)   # AC node
        dc_bus = (converter[i, 2], false)  # DC node
        push!(graph[ac_bus], dc_bus)
        push!(graph[dc_bus], ac_bus)
    end
    
    # Find the slack bus of the AC system (assuming type 3 nodes are slack buses)
    slack_bus_idx = findfirst(x -> x == 3.0, busAC[:, 2])
    if slack_bus_idx === nothing
        slack_bus = (busAC[1, 1], true)  # If no slack bus is found, use the first AC node
    else
        slack_bus = (busAC[slack_bus_idx, 1], true)
    end
    
    # Use BFS for renumbering (treating the entire system as a connected system)
    visited = Set{Tuple{Float64, Bool}}()
    queue = [slack_bus]  # Start BFS from the slack bus
    new_number = 1
    
    while !isempty(queue)
        current = popfirst!(queue)
        
        if current in visited
            continue
        end
        
        push!(visited, current)
        old_to_new[current] = new_number
        new_number += 1
        
        # Add all unvisited neighbor nodes to the queue
        for neighbor in graph[current]
            if !(neighbor in visited)
                push!(queue, neighbor)
            end
        end
    end
    
    # Handle nodes that may not be connected to the main network
    for bus in all_buses
        if !(bus in keys(old_to_new))
            old_to_new[bus] = new_number
            new_number += 1
        end
    end
    
    # Create new data structures
    new_busAC = copy(busAC)
    new_busDC = copy(busDC)
    new_branchAC = copy(branchAC)
    new_branchDC = copy(branchDC)
    new_converter = copy(converter)
    
    # Copy load and generator data (if exists)
    new_loadAC = loadAC !== nothing ? copy(loadAC) : nothing
    new_loadDC = loadDC !== nothing ? copy(loadDC) : nothing
    new_genAC = genAC !== nothing ? copy(genAC) : nothing
    new_genDC = genDC !== nothing ? copy(genDC) : nothing
    new_pvarray = pvarray !== nothing ? copy(pvarray) : nothing
    
    # Copy AC PV system data (if exists)
    new_pv_acsystem = pv_acsystem !== nothing ? copy(pv_acsystem) : nothing
    
    # Copy storage data (if exists)
    new_storage = storage !== nothing ? copy(storage) : nothing
    
    # Update AC bus numbers
    for i in 1:size(busAC, 1)
        old_num = busAC[i, 1]
        new_busAC[i, 1] = old_to_new[(old_num, true)]
    end
    
    # Update DC bus numbers
    for i in 1:size(busDC, 1)
        old_num = busDC[i, 1]
        new_busDC[i, 1] = old_to_new[(old_num, false)]
    end
    
    # Update AC branch numbers and ensure from < to
    for i in 1:size(branchAC, 1)
        from_old = branchAC[i, 1]
        to_old = branchAC[i, 2]
        from_new = old_to_new[(from_old, true)]
        to_new = old_to_new[(to_old, true)]
        
        # Ensure smaller node number comes first
        if from_new < to_new
            new_branchAC[i, 1] = from_new
            new_branchAC[i, 2] = to_new
        else
            new_branchAC[i, 1] = to_new
            new_branchAC[i, 2] = from_new
        end
    end
    
    # Update DC branch numbers and ensure from < to
    for i in 1:size(branchDC, 1)
        from_old = branchDC[i, 1]
        to_old = branchDC[i, 2]
        from_new = old_to_new[(from_old, false)]
        to_new = old_to_new[(to_old, false)]
        
        # Ensure smaller node number comes first
        if from_new < to_new
            new_branchDC[i, 1] = from_new
            new_branchDC[i, 2] = to_new
        else
            new_branchDC[i, 1] = to_new
            new_branchDC[i, 2] = from_new
        end
    end
    
    # Update converter numbers and ensure from < to
    for i in 1:size(converter, 1)
        ac_old = converter[i, 1]
        dc_old = converter[i, 2]
        ac_new = old_to_new[(ac_old, true)]
        dc_new = old_to_new[(dc_old, false)]
        
        # Ensure smaller node number comes first
        if ac_new < dc_new
            new_converter[i, 1] = ac_new
            new_converter[i, 2] = dc_new
        else
            new_converter[i, 1] = dc_new
            new_converter[i, 2] = ac_new
        end
    end
    
    # Update AC load numbers (if exists)
    if new_loadAC !== nothing
        for i in 1:size(new_loadAC, 1)
            bus_num = new_loadAC[i, 2]  # Second column is the bus number
            if haskey(old_to_new, (bus_num, true))
                new_loadAC[i, 2] = old_to_new[(bus_num, true)]
            end
        end
    end
    
    # Update DC load numbers (if exists)
    if new_loadDC !== nothing
        for i in 1:size(new_loadDC, 1)
            bus_num = new_loadDC[i, 2]  # Second column is the bus number
            if haskey(old_to_new, (bus_num, false))
                new_loadDC[i, 2] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # Update AC generator numbers (if exists)
    if new_genAC !== nothing
        for i in 1:size(new_genAC, 1)
            bus_num = new_genAC[i, 1]  # First column is the bus number
            if haskey(old_to_new, (bus_num, true))
                new_genAC[i, 1] = old_to_new[(bus_num, true)]
            end
        end
    end
    
    # Update DC generator numbers (if exists)
    if new_genDC !== nothing
        for i in 1:size(new_genDC, 1)
            bus_num = new_genDC[i, 1]  # First column is the bus number
            if haskey(old_to_new, (bus_num, false))
                new_genDC[i, 1] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # Update PV array numbers (if exists)
    if new_pvarray !== nothing
        for i in 1:size(new_pvarray, 1)
            bus_num = new_pvarray[i, 2]  # Second column is the bus number
            # Assume PV arrays are connected to DC nodes
            if haskey(old_to_new, (bus_num, false))
                new_pvarray[i, 2] = old_to_new[(bus_num, false)]
            end
        end
    end
    
    # Update AC PV system numbers (if exists)
    if new_pv_acsystem !== nothing
        for i in 1:size(new_pv_acsystem, 1)
            bus_num = new_pv_acsystem[i, 2]  # Second column is the bus number
            # AC PV systems are connected to AC nodes
            if haskey(old_to_new, (bus_num, true))
                new_pv_acsystem[i, 2] = old_to_new[(bus_num, true)]
            else
                @warn "Cannot find AC bus with ID $bus_num, unable to update AC PV system number"
            end
        end
    end
    
    # Update storage device numbers (if exists)
    if new_storage !== nothing
        for i in 1:size(new_storage, 1)
            bus_num = new_storage[i, 1]  # Assume first column is the bus number (ESS_BUS)
            # Determine whether the storage device is connected to AC or DC bus
            # This needs to be determined based on actual situation, I assume storage devices may connect to AC or DC buses
            if haskey(old_to_new, (bus_num, false))
                # Connected to DC bus
                new_storage[i, 1] = old_to_new[(bus_num, false)]
            else
                @warn "Cannot find bus with ID $bus_num, unable to update storage device number"
            end
        end
    end
    
    # Sort the results
    sort_idx_ac = sortperm(new_busAC[:, 1])
    sort_idx_dc = sortperm(new_busDC[:, 1])
    new_busAC = new_busAC[sort_idx_ac, :]
    new_busDC = new_busDC[sort_idx_dc, :]
    
    # Sort branches (by starting node)
    sort_idx_branchAC = sortperm(new_branchAC[:, 1])
    sort_idx_branchDC = sortperm(new_branchDC[:, 1])
    sort_idx_converter = sortperm(new_converter[:, 1])
    new_branchAC = new_branchAC[sort_idx_branchAC, :]
    new_branchDC = new_branchDC[sort_idx_branchDC, :]
    new_converter = new_converter[sort_idx_converter, :]
    
    # Sort load and generator data (if exists)
    if new_loadAC !== nothing
        sort_idx_loadAC = sortperm(new_loadAC[:, 1])  # Sort by bus number
        new_loadAC = new_loadAC[sort_idx_loadAC, :]
    end
    
    if new_loadDC !== nothing
        sort_idx_loadDC = sortperm(new_loadDC[:, 1])  # Sort by bus number
        new_loadDC = new_loadDC[sort_idx_loadDC, :]
    end
    
    if new_genAC !== nothing
        sort_idx_genAC = sortperm(new_genAC[:, 1])  # Sort by bus number
        new_genAC = new_genAC[sort_idx_genAC, :]
    end
    
    if new_genDC !== nothing
        sort_idx_genDC = sortperm(new_genDC[:, 1])  # Sort by bus number
        new_genDC = new_genDC[sort_idx_genDC, :]
    end
    
    if new_pvarray !== nothing
        sort_idx_pvarray = sortperm(new_pvarray[:, 2])  # Sort by bus number
        new_pvarray = new_pvarray[sort_idx_pvarray, :]
    end
    
    # Sort AC PV system data (if exists)
    if new_pv_acsystem !== nothing
        sort_idx_pv_acsystem = sortperm(new_pv_acsystem[:, 2])  # Sort by bus number
        new_pv_acsystem = new_pv_acsystem[sort_idx_pv_acsystem, :]
    end
    
    # Sort storage device data (if exists)
    if new_storage !== nothing
        sort_idx_storage = sortperm(new_storage[:, 1])  # Sort by bus number
        new_storage = new_storage[sort_idx_storage, :]
    end
    
    # Create new jpc structure
    new_jpc = (
        busAC = new_busAC,
        busDC = new_busDC,
        branchAC = new_branchAC,
        branchDC = new_branchDC,
        converter = new_converter
    )
    
    # Create complete mapping dictionary (for easy lookup)
    complete_mapping = Dict{Tuple{Float64, String}, Int}()
    for ((bus, is_ac), new_num) in old_to_new
        system_type = is_ac ? "AC" : "DC"
        complete_mapping[(bus, system_type)] = new_num
    end
    
    # Add mapping dictionary to the result
    new_jpc = merge(new_jpc, (mapping = complete_mapping,))
    
    # Add load and generator data (if exists)
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
    
    # Add AC PV system data (if exists)
    if pv_acsystem !== nothing
        new_jpc = merge(new_jpc, (pv_acsystem = new_pv_acsystem,))
    end
    
    # Add storage data (if exists)
    if storage !== nothing
        new_jpc = merge(new_jpc, (storage = new_storage,))
    end
    
    # Create mapping dictionaries for AC nodes and DC nodes
    ac_node_mapping = Dict{Float64, Int}()
    dc_node_mapping = Dict{Float64, Int}()
    
    # Fill mapping dictionaries
    for ((bus, is_ac), new_num) in old_to_new
        if is_ac
            ac_node_mapping[bus] = new_num
        else
            dc_node_mapping[bus] = new_num
        end
    end
    
    return new_jpc, ac_node_mapping, dc_node_mapping
end
