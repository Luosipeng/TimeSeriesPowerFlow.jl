"""
    find_islands(jpc::JPC)

Find electrically isolated islands in an AC power system.

This function identifies connected components in the power network by analyzing
the branch connections between buses. It also identifies isolated PQ buses and
islands that contain only PQ buses (no generators or reference buses).

# Arguments
- `jpc::JPC`: Power system case data structure

# Returns
- `Vector{Vector{Int}}`: List of islands (each containing bus IDs)
- `Vector{Int}`: List of isolated buses (PQ buses with no connections)
"""
function find_islands(jpc::JPC)
    # Get all bus numbers
    bus_ids = jpc.busAC[:, 1]
    n_bus = length(bus_ids)
    
    # Create mapping from bus numbers to indices
    bus_id_to_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_ids)
        bus_id_to_idx[Int(bus_id)] = idx
    end
    
    # Get bus type indices
    PQ, PV, REF = idx_bus()[1:3]
    
    # Create adjacency matrix
    adj = zeros(Int, n_bus, n_bus)
    
    # Only consider branches with status = 1
    for i in 1:size(jpc.branchAC, 1)
        branch = jpc.branchAC[i, :]
        if branch[11] == 1  # BR_STATUS = 11
            f_bus_id = Int(branch[1])
            t_bus_id = Int(branch[2])
            
            # Use mapping to get correct indices
            f_idx = bus_id_to_idx[f_bus_id]
            t_idx = bus_id_to_idx[t_bus_id]
            
            adj[f_idx, t_idx] = 1
            adj[t_idx, f_idx] = 1
        end
    end
    
    # Use DFS to find connected components
    visited = falses(n_bus)
    groups = Vector{Vector{Int}}()
    
    # First find all connected components
    for i in 1:n_bus
        if !visited[i]
            component = Int[]
            stack = [i]
            visited[i] = true
            
            while !isempty(stack)
                node_idx = pop!(stack)
                node_id = Int(bus_ids[node_idx])  # Convert back to original bus number
                push!(component, node_id)
                
                # Check adjacent nodes
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
    
    # Find isolated nodes and all-PQ node groups
    isolated = Int[]
    groups_to_remove = Int[]
    
    # Check each group
    for (idx, group) in enumerate(groups)
        if length(group) == 1
            # Handle single-node groups
            node_id = group[1]
            node_idx = bus_id_to_idx[node_id]
            connections = sum(adj[node_idx, :])
            
            # Find the corresponding bus row in jpc.busAC
            bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
            bus_type = Int(jpc.busAC[bus_row, 2])
            
            if connections == 0 && bus_type == PQ
                # Isolated PQ nodes move to isolated list
                push!(isolated, node_id)
                push!(groups_to_remove, idx)
            end
        else
            # Check if multi-node group consists of only PQ nodes
            all_pq = true
            for node_id in group
                # Find the corresponding bus row in jpc.busAC
                bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
                bus_type = Int(jpc.busAC[bus_row, 2])
                
                if bus_type == PV || bus_type == REF
                    all_pq = false
                    break
                end
            end
            
            if all_pq
                # If all are PQ nodes, move the entire group to isolated
                append!(isolated, group)
                push!(groups_to_remove, idx)
            end
        end
    end
    
    # Delete groups from back to front to avoid index change issues
    sort!(groups_to_remove, rev=true)
    for idx in groups_to_remove
        deleteat!(groups, idx)
    end
      
    # Return results
    return groups, isolated
end


"""
    find_islands_acdc(jpc::JPC)

Find electrically isolated islands in a hybrid AC-DC power system.

This function extends the island detection to handle both AC and DC subsystems,
including converter connections between them. It identifies connected components
across the entire network and classifies islands based on their generation capabilities.

# Arguments
- `jpc::JPC`: Power system case data structure containing both AC and DC components

# Returns
- `Vector{Vector{Int}}`: List of AC islands (each containing bus IDs)
- `Vector{Int}`: List of isolated AC buses (PQ buses with no viable connections)

# Notes
- Islands with only PQ buses and no generation capability (either directly or through
  DC connections) are considered isolated
- If no DC buses exist, the function falls back to the standard AC island detection
"""
function find_islands_acdc(jpc::JPC)
    # Get all AC and DC bus numbers
    bus_ac_ids = jpc.busAC[:, 1]
    bus_dc_ids = size(jpc.busDC, 1) > 0 ? jpc.busDC[:, 1] : Int[]
    
    n_bus_ac = length(bus_ac_ids)
    n_bus_dc = length(bus_dc_ids)
    n_bus_total = n_bus_ac + n_bus_dc
    
    # If there are no DC buses, use the original find_islands function
    if n_bus_dc == 0
        return find_islands(jpc)
    end
    
    # Create unified numbering system
    unified_ids = vcat(bus_ac_ids, bus_dc_ids)
    
    # Create mapping from bus numbers to indices
    bus_id_to_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(unified_ids)
        bus_id_to_idx[Int(bus_id)] = idx
    end
    
    # Create mapping from indices to original system (AC/DC)
    idx_to_system = Dict{Int, Symbol}()
    for i in 1:n_bus_ac
        idx_to_system[i] = :AC
    end
    for i in (n_bus_ac+1):n_bus_total
        idx_to_system[i] = :DC
    end
    
    # Create mapping from original numbers to unified indices
    ac_id_to_unified_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_ac_ids)
        ac_id_to_unified_idx[Int(bus_id)] = idx
    end
    
    dc_id_to_unified_idx = Dict{Int, Int}()
    for (idx, bus_id) in enumerate(bus_dc_ids)
        dc_id_to_unified_idx[Int(bus_id)] = n_bus_ac + idx
    end
    
    # Get bus type indices
    PQ, PV, REF = idx_bus()[1:3]
    
    # Create adjacency matrix
    adj = zeros(Int, n_bus_total, n_bus_total)
    
    # Process AC branches
    for i in 1:size(jpc.branchAC, 1)
        branch = jpc.branchAC[i, :]
        if branch[11] == 1  # BR_STATUS = 11
            f_bus_id = Int(branch[1])
            t_bus_id = Int(branch[2])
            
            # Use mapping to get correct indices
            f_idx = ac_id_to_unified_idx[f_bus_id]
            t_idx = ac_id_to_unified_idx[t_bus_id]
            
            adj[f_idx, t_idx] = 1
            adj[t_idx, f_idx] = 1
        end
    end
    
    # Process DC branches
    if size(jpc.branchDC, 1) > 0
        for i in 1:size(jpc.branchDC, 1)
            branch = jpc.branchDC[i, :]
            if branch[11] == 1  # Assume DC branches also have status field at same position as AC
                f_bus_id = Int(branch[1])
                t_bus_id = Int(branch[2])
                
                # Use mapping to get correct indices
                if haskey(dc_id_to_unified_idx, f_bus_id) && haskey(dc_id_to_unified_idx, t_bus_id)
                    f_idx = dc_id_to_unified_idx[f_bus_id]
                    t_idx = dc_id_to_unified_idx[t_bus_id]
                    
                    adj[f_idx, t_idx] = 1
                    adj[t_idx, f_idx] = 1
                end
            end
        end
    end
    
    # Process converters (connecting AC and DC systems)
    if size(jpc.converter, 1) > 0
        for i in 1:size(jpc.converter, 1)
            converter = jpc.converter[i, :]
            # Assume converters are in service by default, unless explicitly set to 0
            is_active = true
            if size(converter, 1) >= 3
                is_active = converter[3] != 0
            end
            
            if is_active
                ac_bus_id = Int(converter[1])
                dc_bus_id = Int(converter[2])
                
                # Check if these buses exist
                if haskey(ac_id_to_unified_idx, ac_bus_id) && haskey(dc_id_to_unified_idx, dc_bus_id)
                    ac_idx = ac_id_to_unified_idx[ac_bus_id]
                    dc_idx = dc_id_to_unified_idx[dc_bus_id]
                    
                    adj[ac_idx, dc_idx] = 1
                    adj[dc_idx, ac_idx] = 1
                end
            end
        end
    end
    
    # Use DFS to find connected components
    visited = falses(n_bus_total)
    unified_groups = Vector{Vector{Int}}()
    
    # First find all connected components (including AC and DC nodes)
    for i in 1:n_bus_total
        if !visited[i]
            component = Int[]
            stack = [i]
            visited[i] = true
            
            while !isempty(stack)
                node_idx = pop!(stack)
                node_id = Int(unified_ids[node_idx])  # Get original bus number
                push!(component, node_idx)  # Store index rather than ID
                
                # Check adjacent nodes
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
    
    # Convert unified indices back to original bus numbers, and separate AC and DC parts
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
    
    # Find isolated nodes and all-PQ node groups
    isolated = Int[]
    groups_to_remove = Int[]
    
    # Check each group
    for (idx, group) in enumerate(groups)
        if length(group) == 1
            # Handle single-node groups
            node_id = group[1]
            node_idx = ac_id_to_unified_idx[node_id]
            
            # Check if this node has connections (including to DC system)
            connections = sum(adj[node_idx, :])
            
            # Find the corresponding bus row in jpc.busAC
            bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, 2])
                
                if connections == 0 && bus_type == PQ
                    # Isolated PQ nodes move to isolated list
                    push!(isolated, node_id)
                    push!(groups_to_remove, idx)
                end
            end
        else
            # Check if multi-node group consists of only PQ nodes
            all_pq = true
            has_generator = false
            
            for node_id in group
                # Find the corresponding bus row in jpc.busAC
                bus_row = findfirst(x -> Int(x) == node_id, jpc.busAC[:, 1])
                if bus_row !== nothing
                    bus_type = Int(jpc.busAC[bus_row, 2])
                    
                    if bus_type == PV || bus_type == REF
                        all_pq = false
                        break
                    end
                end
                
                # Check if there are generators connected to this bus
                for j in 1:size(jpc.genAC, 1)
                    gen_bus = Int(jpc.genAC[j, 1])  # Assume column 2 is bus number
                    if gen_bus == node_id && jpc.genAC[j, 8] == 1  # Assume column 8 is status
                        has_generator = true
                        all_pq = false
                        break
                    end
                end
                
                if !all_pq
                    break
                end
            end
            
            # Find DC buses connected to this group
            connected_dc_buses = Set{Int}()
            for ac_bus_id in group
                if haskey(ac_id_to_unified_idx, ac_bus_id)
                    ac_idx = ac_id_to_unified_idx[ac_bus_id]
                    
                    # Look for DC nodes directly connected to this AC node
                    for dc_idx in (n_bus_ac+1):n_bus_total
                        if adj[ac_idx, dc_idx] == 1
                            dc_bus_id = unified_ids[dc_idx]
                            push!(connected_dc_buses, dc_bus_id)
                        end
                    end
                end
            end
            
            # Check if these DC buses have generators
            if !isempty(connected_dc_buses) && size(jpc.sgenDC, 1) > 0
                for j in 1:size(jpc.sgenDC, 1)
                    sgen_bus = Int(jpc.sgenDC[j, 1])  # Assume column 1 is bus number
                    # Assume column 3 of sgenDC is status field, if not present assume in service
                    is_active = size(jpc.sgenDC, 2) >= 3 ? jpc.sgenDC[j, 3] != 0 : true
                    
                    if sgen_bus in connected_dc_buses && is_active
                        has_generator = true
                        all_pq = false
                        break
                    end
                end
            end
            
            if all_pq && !has_generator
                # If all are PQ nodes and no generators, move the entire group to isolated
                append!(isolated, group)
                push!(groups_to_remove, idx)
            end
        end
    end
    
    # Delete groups from back to front to avoid index change issues
    sort!(groups_to_remove, rev=true)
    for idx in groups_to_remove
        deleteat!(groups, idx)
    end
      
    # Return results
    return groups, isolated
end
