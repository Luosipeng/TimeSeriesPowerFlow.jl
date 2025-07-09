"""
    calculate_bus(net, jpc, sequence, slack_bus, opt)

Calculate bus parameters for the specified sequence network.

# Arguments
- `net`: Network data structure containing bus information
- `jpc`: MATPOWER-style power flow case structure
- `sequence`: Sequence type (1 for positive, 2 for negative, 0 for zero sequence)
- `slack_bus`: Index of the slack bus
- `opt`: Options dictionary containing power flow parameters

# Returns
- `jpc`: Updated jpc structure with bus parameters
"""
function calculate_bus(net, jpc, sequence, slack_bus, opt)

    # Initialize
    nb = size(net["bus"], 1)
    bus = zeros(nb, 13)

    if sequence == 1
        # Positive sequence network
        # Process bus data
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # If all elements are the same, assign 1 uniformly
            bus[:,ZONE] .= 1
        else
            # If elements are not all the same, group identical values and assign 1, 2, 3... sequentially
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # Map each zone to the corresponding number (1, 2, 3...)
            bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
        end
        bus[:,BASE_KV] = bus1.vn_kv 
        bus[:,VM] .= 1.0
        bus[:,VA] .= 0.0
        bus[:,VMAX] = bus1.max_vm_pu
        bus[:,VMIN] = bus1.min_vm_pu
        
    elseif sequence == 2
        # Negative sequence network
        # Process bus data
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # If all elements are the same, assign 1 uniformly
            bus[:,ZONE] .= 1
        else
            # If elements are not all the same, group identical values and assign 1, 2, 3... sequentially
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # Map each zone to the corresponding number (1, 2, 3...)
            bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
        end
        bus[:,BASE_KV] = bus1.vn_kv 
        bus[:,VM].= 0.0
        bus[:,VA] .= 0.0
        bus[:,VMAX] = bus1.max_vm_pu
        bus[:,VMIN] = bus1.min_vm_pu
        
    else 
        # Zero sequence network
        # Process bus data
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # If all elements are the same, assign 1 uniformly
            bus[:,ZONE] .= 1
        else
            # If elements are not all the same, group identical values and assign 1, 2, 3... sequentially
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # Map each zone to the corresponding number (1, 2, 3...)
            bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
        end
        bus[:,BASE_KV] = bus1.vn_kv 
        bus[:,VM].= 0.0
        bus[:,VA] .= 0.0
        bus[:,VMAX] = bus1.max_vm_pu
        bus[:,VMIN] = bus1.min_vm_pu

    end
    jpc["bus"] = bus

    return jpc
end

"""
    calculate_bus(net, jpc, slack_bus, opt)

Calculate bus parameters for the default network model (positive sequence).

# Arguments
- `net`: Network data structure containing bus information
- `jpc`: MATPOWER-style power flow case structure
- `slack_bus`: Index of the slack bus
- `opt`: Options dictionary containing power flow parameters

# Returns
- `jpc`: Updated jpc structure with bus parameters
"""
function calculate_bus(net, jpc, slack_bus, opt)

    # Initialize
    nb = size(net["bus"], 1)
    bus = zeros(nb, 13)

    # Positive sequence network
    # Process bus data
    bus1 = deepcopy(net["bus"])  
    bus[:, BUS_I] = bus1.index
    bus[:, BUS_TYPE] .= 1
    bus[slack_bus, BUS_TYPE] .= 3
    if all(x -> x == bus1.zone[1], bus1.zone)
        # If all elements are the same, assign 1 uniformly
        bus[:,ZONE] .= 1
    else
        # If elements are not all the same, group identical values and assign 1, 2, 3... sequentially
        unique_zones = unique(bus1.zone)
        zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
        
        # Map each zone to the corresponding number (1, 2, 3...)
        bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
    end
    bus[:,BASE_KV] = bus1.vn_kv 
    bus[:,VM] .= 1.0
    bus[:,VA] .= 0.0
    bus[:,VMAX] = bus1.max_vm_pu
    bus[:,VMIN] = bus1.min_vm_pu
    
    jpc["bus"] = bus

    return jpc
end

"""
    add_grid_external_sc_impedance(jpc_new, external_grid)

Add external grid short-circuit impedance to the network model.

# Arguments
- `jpc_new`: MATPOWER-style power flow case structure
- `external_grid`: External grid data containing short circuit parameters

# Returns
- Tuple of (gs, bs): Conductance and susceptance values added to the external grid bus
"""
function add_grid_external_sc_impedance(jpc_new, external_grid)
    
    external_bus = external_grid.bus
    c = 1.1
    s_sc = external_grid.s_sc_max_mva/jpc_new["baseMVA"]
    rx = external_grid.rx_max
    z_grid = c / (s_sc/3)
    x_grid = z_grid./sqrt.(1 .+rx.^2)
    r_grid = x_grid * rx

    Y_grid = 1 ./ (r_grid .+ 1im*x_grid)
    buses, gs, bs = sum_by_group(external_bus, real(Y_grid), imag(Y_grid))
    jpc_new["bus"][external_bus, GS] .= gs * jpc_new["baseMVA"]
    jpc_new["bus"][external_bus, BS] .= bs * jpc_new["baseMVA"]

    return gs * jpc_new["baseMVA"], bs * jpc_new["baseMVA"]
end
