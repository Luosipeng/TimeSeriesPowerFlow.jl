"""
    create_node_mapping(case::JuliaPowerCase)

Create numbering mapping for nodes, detect duplicate nodes.
Returns a dictionary mapping node names to IDs.
"""
function create_node_mapping(case::JuliaPowerCase)
    # Extract all node names from the case
    nodes = [bus.name for bus in case.busesAC]
    
    count_dict = Dict{String, Int}()
    duplicates = String[]
    
    # Count node occurrences, detect duplicate nodes
    for node in nodes
        node_str = string(node)
        count_dict[node_str] = get(count_dict, node_str, 0) + 1
        if count_dict[node_str] > 1 && !(node_str in duplicates)
            push!(duplicates, node_str)
        end
    end
    
    # Output warnings for duplicate nodes
    if !isempty(duplicates)
        println("\nWarning: Duplicate node names found:")
        for name in duplicates
            println(" - ", name, " (appears ", count_dict[name], " times)")
        end
        println()
    end
    
    # Create node to ID mapping
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

Filter out connections with status value of true.
Returns a list of active connections.
"""
function filter_active_edges(case::JuliaPowerCase)
    # Extract all edges from the case
    active_edges = []
    
    # Process AC lines
    for line in case.branchesAC
        if line.in_service
            push!(active_edges, (line.from_bus, line.to_bus, line))
        end
    end
    
    # Process transformers
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            push!(active_edges, (transformer.hv_bus, transformer.lv_bus, transformer))
        end
    end
    
    # Process circuit breakers
    for hvcb in case.hvcbs
        if hvcb.closed
            push!(active_edges, (hvcb.bus_from, hvcb.bus_to, hvcb))
        end
    end
    
    return active_edges
end

"""
    update_edge_ids(edges, node_dict)

Update the start and end IDs of the connection table.
Returns a list of ID pairs for edges.
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

Determine if an element is a switch type.
"""
function is_switch_element(element)
    # Check if the element is a circuit breaker type
    if typeof(element) == DistributionSystem.HighVoltageCircuitBreaker
        return true
    end
    
    # Check for other possible switch types
    if hasfield(typeof(element), :type)
        type_str = lowercase(strip(string(element.type)))
        switch_types = ["switch", "circuit breaker", "isolating switch", "load switch", "knife switch", "disconnector"]
        is_switch = any(type -> type == type_str || occursin(type, type_str), switch_types)
        is_not_room = !occursin("room", type_str) && !occursin("chamber", type_str) && !occursin("station", type_str)
        return is_switch && is_not_room
    end
    
    return false
end

"""
    is_substation(element)

Determine if an element is a substation type.
"""
function is_substation(element)
    if hasfield(typeof(element), :type)
        type_str = lowercase(strip(string(element.type)))
        return occursin("substation", type_str)
    end
    return false
end

"""
    count_switch_ports(case::JuliaPowerCase)

Count the number of switch ports.
Returns a dictionary mapping node IDs to port counts.
"""
function count_switch_ports(case::JuliaPowerCase)
    port_counts = Dict{Int, Int}()
    
    # Initialize all nodes' connection counts to 0
    for bus in case.busesAC
        port_counts[bus.bus_id] = 0
    end
    
    # Count the connections for each node
    # Process AC lines
    for line in case.branchesAC
        if line.in_service
            port_counts[line.from_bus] = get(port_counts, line.from_bus, 0) + 1
            port_counts[line.to_bus] = get(port_counts, line.to_bus, 0) + 1
        end
    end
    
    # Process transformers
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            port_counts[transformer.hv_bus] = get(port_counts, transformer.hv_bus, 0) + 1
            port_counts[transformer.lv_bus] = get(port_counts, transformer.lv_bus, 0) + 1
        end
    end
    
    # Process circuit breakers
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

Create virtual node data.
"""
function create_virtual_node(case::JuliaPowerCase, bus_id::Int, virtual_name::String)
    # Find the original node
    original_bus = nothing
    for bus in case.busesAC
        if bus.bus_id == bus_id
            original_bus = bus
            break
        end
    end
    
    if original_bus === nothing
        error("Cannot find node with ID $bus_id")
    end
    
    # Create virtual node (copy attributes from original node)
    virtual_bus = deepcopy(original_bus)
    virtual_bus.index = length(case.busesAC) + 1  # New ID
    virtual_bus.bus_id = length(case.busesAC) + 1  # New ID
    virtual_bus.name = virtual_name
    
    return virtual_bus
end

"""
    create_virtual_connection(case::JuliaPowerCase, from_bus::Int, to_bus::Int)

Create virtual connection data.
"""
function create_virtual_connection(case::JuliaPowerCase, from_bus::Int, to_bus::Int)
    # Create a new circuit breaker as a virtual connection
    virtual_cb = DistributionSystem.HighVoltageCircuitBreaker(
        length(case.hvcbs) + 1,  # New ID
        "VIRTUAL_CB_$(from_bus)_$(to_bus)",  # Name
        from_bus,  # Start node
        to_bus,    # End node
        "l",       # Type
        0,         # Sequence number
        true,     # Closed status
        "CB",      # Type
        0.0,       # Rated current
        true       # Operating status
    )
    
    return virtual_cb
end

"""
    process_switch_nodes(case::JuliaPowerCase)

Process switch nodes, add virtual nodes.
Returns the updated case object and virtual connection markers.
"""
function process_switch_nodes(case::JuliaPowerCase)
    # Deep copy the case for modification
    new_case = deepcopy(case)
    
    # Create a set to mark virtual connections
    virtual_connections = Set{Tuple{Int, Int}}()
    
    # Get port counts
    port_counts = count_switch_ports(case)
    
    # Create a set to track processed circuit breakers
    processed_hvcbs = Set{Int}()
    
    # Find all switch elements
    switch_elements = []
    
    # Check circuit breakers
    for hvcb in case.hvcbs
        if hvcb.closed && !(hvcb.index in processed_hvcbs)
            bus_from = hvcb.bus_from
            bus_to = hvcb.bus_to
            
            # Check if both end nodes meet the conditions
            from_eligible = get(port_counts, bus_from, 0) > 1
            to_eligible = get(port_counts, bus_to, 0) > 1
            
            if from_eligible && to_eligible
                # If both ends meet the conditions, choose the end with more ports
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
            
            # Mark this circuit breaker as processed
            push!(processed_hvcbs, hvcb.index)
        end
    end
    
    # Process each switch node
    for (bus_id, element) in switch_elements
        port_count = get(port_counts, bus_id, 0)
        
        # Process cases with any number of ports (at least 2)
        if port_count >= 2
            # Get original node name
            original_name = ""
            for bus in case.busesAC
                if bus.bus_id == bus_id
                    original_name = bus.name
                    break
                end
            end
            
            # Create a virtual node
            virtual_name = string(original_name) * "_virtual_node"
            virtual_bus = create_virtual_node(new_case, bus_id, virtual_name)
            push!(new_case.busesAC, virtual_bus)
            
            # Connect original node and virtual node
            virtual_cb = create_virtual_connection(new_case, bus_id, virtual_bus.bus_id)
            push!(new_case.hvcbs, virtual_cb)
            
            # Add first virtual connection (original node-virtual node)
            if bus_id < virtual_bus.bus_id
                push!(virtual_connections, (bus_id, virtual_bus.bus_id))
            else
                push!(virtual_connections, (virtual_bus.bus_id, bus_id))
            end
            
            # Get the other end node of the switch element
            other_bus_id = element.bus_from == bus_id ? element.bus_to : element.bus_from
            
            # Create new circuit breaker connecting virtual node and the other end node
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
            
            # Add second virtual connection (virtual node-other end node)
            if virtual_bus.bus_id < other_bus_id
                push!(virtual_connections, (virtual_bus.bus_id, other_bus_id))
            else
                push!(virtual_connections, (other_bus_id, virtual_bus.bus_id))
            end
            
            # Set the original circuit breaker to not closed
            for j in 1:length(new_case.hvcbs)
                if new_case.hvcbs[j].index == element.index
                    new_case.hvcbs[j].closed = false
                    break
                end
            end
        end
    end
    
    # Update node name to ID mapping
    new_case.bus_name_to_id = Dict{String, Int}()
    for bus in new_case.busesAC
        new_case.bus_name_to_id[bus.name] = bus.bus_id
    end
    
    return new_case, virtual_connections
end



"""
    identify_partitions(case::JuliaPowerCase)

Identify partitions in the network.
Returns node partitions and edge partitions.
"""
function identify_partitions(case::JuliaPowerCase)
    # Create graph structure
    G = SimpleGraph(length(case.busesAC))
    
    # Add edges
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
    
    # Find connected components
    components = connected_components(G)
    
    # Create node partitions and edge partitions
    node_partitions = Dict{Int, Int}()
    edge_partitions = Dict{Tuple{Int, Int}, Int}()
    
    # Assign partitions to each node
    for (partition_id, component) in enumerate(components)
        for node_id in component
            node_partitions[node_id] = partition_id
        end
    end
    
    # Assign partitions to each edge
    for line in case.branchesAC
        if line.in_service
            from_partition = get(node_partitions, line.from_bus, 0)
            to_partition = get(node_partitions, line.to_bus, 0)
            
            # Edge should belong to the same partition
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

Get the source and target nodes of an edge.
"""
function get_edge_endpoints(e)
    # Try different methods to get edge endpoints
    try
        return (e.src, e.dst)
    catch
        try
            # If edge is a tuple, return it directly
            return e
        catch
            # If all above methods fail, try to convert edge to string and parse
            edge_str = string(e)
            # Assume format is "Edge 1 => 2" or similar
            m = match(r"Edge\s+(\d+)\s*=>\s*(\d+)", edge_str)
            if m !== nothing
                return (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
            end
            # If all methods fail, throw error
            error("Cannot get endpoints for edge $e")
        end
    end
end

"""
    dfs_tree(graph, root)

Build a spanning tree using depth-first search.
"""
function dfs_tree(graph, root)
    tree_edges = Vector{Tuple{Int, Int}}()
    visited = falses(nv(graph))
    parent = zeros(Int, nv(graph))
    
    # Helper function, perform DFS
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
    
    # Perform DFS for each connected component
    for v in 1:nv(graph)
        if !visited[v]
            dfs_helper(v)
        end
    end
    
    return tree_edges, parent
end

"""
    find_path_to_lca(u, v, parent)

Find the path from nodes u and v to their lowest common ancestor.
"""
function find_path_to_lca(u, v, parent)
    # Path from u to root
    path_u_to_root = Vector{Int}()
    current = u
    while current != 0
        push!(path_u_to_root, current)
        current = parent[current]
    end
    
    # From v up to the lowest common ancestor
    path_v_to_lca = Vector{Int}()
    current = v
    lca_found = false
    lca = 0
    
    while current != 0
        if current in path_u_to_root
            # Found the lowest common ancestor
            lca = current
            lca_found = true
            break
        end
        
        push!(path_v_to_lca, current)
        current = parent[current]
    end
    
    if !lca_found
        # If no common ancestor found, return empty path
        return Vector{Int}(), false
    end
    
    # Truncate path_u_to_root to lca
    path_u_to_lca = Vector{Int}()
    for node in path_u_to_root
        if node == lca
            break
        end
        push!(path_u_to_lca, node)
    end
    
    # Build path from u to v
    path_u_to_v = copy(path_u_to_lca)
    
    # Reverse path_v_to_lca and add to path_u_to_v
    for i in length(path_v_to_lca):-1:1
        push!(path_u_to_v, path_v_to_lca[i])
    end
    
    return path_u_to_v, true
end

"""
    find_fundamental_cycles(G::SimpleGraph)

Find fundamental cycles in an undirected graph using a depth-first search based method.
"""
function find_fundamental_cycles(G::SimpleGraph)
    n = nv(G)  # Get number of nodes in the graph
    
    # Create a mapping from actual node IDs to consecutive indices
    node_to_index = Dict{Int, Int}()
    index_to_node = Dict{Int, Int}()
    
    # Fill the mapping
    index = 1
    for v in vertices(G)
        node_to_index[v] = index
        index_to_node[index] = v
        index += 1
    end
    
    # Create a new graph using consecutive indices
    G_continuous = SimpleGraph(n)
    
    # Add edges from original graph, but using consecutive indices
    for e in edges(G)
        src_node, dst_node = get_edge_endpoints(e)
        add_edge!(G_continuous, node_to_index[src_node], node_to_index[dst_node])
    end
    
    # Now find cycles on the graph with continuous indices
    cycles = Vector{Vector{Int}}()
    
    # For each connected component, find fundamental cycles
    components = connected_components(G_continuous)
    
    for component in components
        if length(component) > 0
            root = component[1]
            tree_edges, parent = dfs_tree(G_continuous, root)
            
            # Collect non-tree edges
            non_tree_edges = Vector{Tuple{Int, Int}}()
            for e in edges(G_continuous)
                src_idx, dst_idx = get_edge_endpoints(e)
                if !((src_idx, dst_idx) in tree_edges || (dst_idx, src_idx) in tree_edges)
                    push!(non_tree_edges, (src_idx, dst_idx))
                end
            end
            
            # For each non-tree edge, find a cycle
            for (u, v) in non_tree_edges
                path_u_to_v, success = find_path_to_lca(u, v, parent)
                
                if !success
                    # If no common ancestor found, graph might be disconnected
                    continue
                end
                
                # Add edge (v,u) to complete the cycle
                push!(path_u_to_v, u)
                
                # Map consecutive indices back to original node IDs
                original_cycle = [index_to_node[idx] for idx in path_u_to_v]
                push!(cycles, original_cycle)
            end
        end
    end
    
    return cycles
end

"""
    create_and_plot_graph_by_partition(case::JuliaPowerCase, node_partitions, edge_partitions, virtual_connections)

Create graphs for each partition and find cycles.
"""
function create_and_plot_graph_by_partition(case::JuliaPowerCase, node_partitions, edge_partitions, virtual_connections)
    # Get unique partition IDs
    partition_ids = unique(values(node_partitions))
    
    # Create and plot graphs for each partition
    cycles_by_partition = Dict{Int, Vector{Vector{Int}}}()
    
    for partition_id in partition_ids
        # Create subgraph for this partition
        partition_nodes = [node_id for (node_id, part_id) in node_partitions if part_id == partition_id]
        partition_edges = [(from, to) for ((from, to), part_id) in edge_partitions if part_id == partition_id]
        
        if isempty(partition_nodes)
            println("Partition $partition_id has no nodes, skipping")
            continue
        end
        
        # Create graph using actual node IDs
        G = SimpleGraph()
        
        # Add nodes
        for node_id in partition_nodes
            add_vertex!(G)
        end
        
        # Create mapping from node ID to graph index
        node_to_vertex = Dict{Int, Int}()
        vertex_to_node = Dict{Int, Int}()
        
        for (i, node_id) in enumerate(partition_nodes)
            node_to_vertex[node_id] = i
            vertex_to_node[i] = node_id
        end
        
        # Add edges using graph indices
        for (from, to) in partition_edges
            if haskey(node_to_vertex, from) && haskey(node_to_vertex, to)
                from_vertex = node_to_vertex[from]
                to_vertex = node_to_vertex[to]
                add_edge!(G, from_vertex, to_vertex)
            else
                println("Warning: One or both endpoints of edge ($from, $to) are not in partition $partition_id")
            end
        end
        
        # Find cycles
        cycles = Vector{Vector{Int}}()
        
        try
            # Use custom cycle detection function
            graph_cycles = find_fundamental_cycles(G)
            
            # Map graph indices back to original node IDs
            for cycle in graph_cycles
                original_cycle = [vertex_to_node[v] for v in cycle]
                push!(cycles, original_cycle)
            end
            
            println("Partition $partition_id: Found $(length(cycles)) cycles")
        catch e
            println("Warning: Cycle detection error for partition $partition_id: $(typeof(e)): $(e)")
            println(stacktrace())
        end
        
        cycles_by_partition[partition_id] = cycles
    end
    
    return cycles_by_partition
end



"""
    write_results_with_partitions(output_file, case::JuliaPowerCase, cycles_by_partition, virtual_connections, node_partitions, edge_partitions)

Write results including partition information. Automatically overwrites existing file.
"""
function write_results_with_partitions(output_file, case::JuliaPowerCase, cycles_by_partition, virtual_connections, node_partitions, edge_partitions)
    # Create result dataframes
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
    
    # Fill node data
    for bus in case.busesAC
        push!(nodes_df, [
            bus.index,
            bus.name,
            "Bus",
            get(node_partitions, bus.index, 0)
        ])
    end
    
    # Fill edge data
    # Process AC lines
    for line in case.branchesAC
        if line.in_service
            from_name = ""
            to_name = ""
            
            # Find node names
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
    
    # Process transformers
    for transformer in case.transformers_2w_etap
        if transformer.in_service
            from_name = ""
            to_name = ""
            
            # Find node names
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
    
    # Process circuit breakers
    for hvcb in case.hvcbs
        if hvcb.closed
            from_name = ""
            to_name = ""
            
            # Find node names
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
    
    # Fill cycle data
    for (partition_id, cycles) in cycles_by_partition
        for (cycle_id, cycle) in enumerate(cycles)
            # Convert node IDs in cycle to names
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
                join(node_names, " -> ") * " -> " * node_names[1]  # Close the cycle
            ])
        end
    end
    
    # If file exists, delete it first
    if isfile(output_file)
        rm(output_file)
    end
    
    # Create a new XLSX file
    XLSX.writetable(output_file, 
        Nodes = (collect(eachcol(nodes_df)), names(nodes_df)),
        Edges = (collect(eachcol(edges_df)), names(edges_df)),
        Cycles = (collect(eachcol(cycles_df)), names(cycles_df))
    )
    
    println("Results saved to $output_file")
    
    return Dict(
        "nodes" => nodes_df,
        "edges" => edges_df,
        "cycles" => cycles_df
    )
end

"""
    generate_partition_report(output_file, case::JuliaPowerCase, node_partitions, edge_partitions)

Generate detailed partition report.
"""
function generate_partition_report(output_file, case::JuliaPowerCase, node_partitions, edge_partitions)
    # Get unique partition IDs
    partition_ids = unique(values(node_partitions))
    
    # Create partition statistics dataframe
    partition_stats = DataFrame(
                Partition_ID = Int[],
        Node_Count = Int[],
        Edge_Count = Int[],
        Load_Count = Int[],
        Total_Load_MW = Float64[],
        Has_Cycles = Bool[]
    )
    
    # Calculate statistics for each partition
    for partition_id in partition_ids
        # Calculate number of nodes in this partition
        nodes_in_partition = [node_id for (node_id, part_id) in node_partitions if part_id == partition_id]
        node_count = length(nodes_in_partition)
        
        # Calculate number of edges in this partition
        edges_in_partition = [(from, to) for ((from, to), part_id) in edge_partitions if part_id == partition_id]
        edge_count = length(edges_in_partition)
        
        # Calculate number of loads and total load in this partition
        load_count = 0
        total_load = 0.0
        
        for load in case.loadsAC
            if load.in_service && load.bus in nodes_in_partition
                load_count += 1
                total_load += load.p_mw
            end
        end
        
        # Check if this partition has cycles
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
    
    # Create partition node details dataframe
    partition_nodes = DataFrame(
        Partition_ID = Int[],
        Node_ID = Int[],
        Node_Name = String[],
        Node_Type = String[],
        Has_Load = Bool[],
        Load_MW = Float64[]
    )
    
    # Fill partition node details
    for (node_id, partition_id) in node_partitions
        node_name = ""
        node_type = ""
        has_load = false
        load_mw = 0.0
        
        # Find node name and type
        for bus in case.busesAC
            if bus.index == node_id
                node_name = bus.name
                node_type = "Bus"
                break
            end
        end
        
        # Check if this node has load
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
    
    # Create partition edge details dataframe
    partition_edges = DataFrame(
        Partition_ID = Int[],
        From_ID = Int[],
        To_ID = Int[],
        From_Name = String[],
        To_Name = String[],
        Edge_Type = String[],
        Is_Virtual = Bool[]
    )
    
    # Fill partition edge details
    for ((from_id, to_id), partition_id) in edge_partitions
        from_name = ""
        to_name = ""
        edge_type = ""
        is_virtual = false
        
        # Find node names
        for bus in case.busesAC
            if bus.index == from_id
                from_name = bus.name
            end
            if bus.index == to_id
                to_name = bus.name
            end
        end
        
        # Determine edge type
        # Check if it's a line
        for line in case.branchesAC
            if line.in_service && ((line.from_bus == from_id && line.to_bus == to_id) || 
                                   (line.from_bus == to_id && line.to_bus == from_id))
                edge_type = "Line"
                break
            end
        end
        
        # Check if it's a transformer
        if edge_type == ""
            for transformer in case.transformers_2w_etap
                if transformer.in_service && ((transformer.hv_bus == from_id && transformer.lv_bus == to_id) || 
                                             (transformer.hv_bus == to_id && transformer.lv_bus == from_id))
                    edge_type = "Transformer"
                    break
                end
            end
        end
        
        # Check if it's a circuit breaker
        if edge_type == ""
            for hvcb in case.hvcbs
                if hvcb.closed && ((hvcb.bus_from == from_id && hvcb.bus_to == to_id) || 
                                  (hvcb.bus_from == to_id && hvcb.bus_to == from_id))
                    edge_type = "CircuitBreaker"
                    
                    # Check if it's a virtual connection
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
    
    # Generate report filename
    report_file = replace(output_file, ".xlsx" => "_partition_report.xlsx")
    
    # If file exists, delete it first
    if isfile(report_file)
        rm(report_file)
    end
    
    # Write report to Excel file
    XLSX.writetable(report_file, 
        Partition_Stats = (collect(eachcol(partition_stats)), names(partition_stats)),
        Partition_Nodes = (collect(eachcol(partition_nodes)), names(partition_nodes)),
        Partition_Edges = (collect(eachcol(partition_edges)), names(partition_edges))
    )
    
    println("Partition report saved to $report_file")
    
    return report_file
end


"""
    extract_edges_from_case(case::JuliaPowerCase)

Extract all edge information from the case.
Returns a list of edges, each edge is a tuple (from_bus, to_bus, edge_object).
"""
function extract_edges_from_case(case::JuliaPowerCase)
    edges = []
    
    # Process AC lines
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
    
    # Process transformers
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
    
    # Process circuit breakers
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

Main function, executes the complete processing workflow.
"""
function topology_analysis(case::JuliaPowerCase; output_file = "./output_result.xlsx", debug=false)
    # Create node mapping
    node_dict = create_node_mapping(case)
    
    # Extract all edges
    edges = extract_edges_from_case(case)
    
    # Count switch ports
    port_counts = count_switch_ports(case)
    
    # Process switch nodes, add virtual nodes, and get virtual connection markers
    new_case, virtual_connections = process_switch_nodes(case)
    
    # Identify partitions
    node_partitions, edge_partitions = identify_partitions(new_case)
    
    # Create and save graphs, while getting cycle information, draw by partition
    cycles = create_and_plot_graph_by_partition(new_case, node_partitions, edge_partitions, virtual_connections)
    
    # Write results
    results = write_results_with_partitions(output_file, new_case, cycles, virtual_connections, node_partitions, edge_partitions)
    
    # Generate detailed partition report
    generate_partition_report(output_file, new_case, node_partitions, edge_partitions)
    
    return results, new_case
end

# Helper function: Merge multiple dataframes into one, add "table" column
function gather_tables_as_one(tables::Vector{Pair{String, DataFrame}})::DataFrame
    # Collect all possible columns
    allcols = Symbol[]
    for (_, df) in tables
        for c in names(df)
            if c ∉ allcols
                push!(allcols, c)
            end
        end
    end
    
    # For each table, ensure it has all columns, fill missing values, add "table" column
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
