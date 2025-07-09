function calculate_bus(net, jpc ,sequence, slack_bus, opt)

    # 初始化
    nb = size(net["bus"], 1)
    bus = zeros(nb, 13)


    if sequence == 1
        # 正序网络
        # 处理母线数据
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # 如果所有元素都相同，统一赋值为1
            bus[:,ZONE] .= 1
        else
            # 如果元素不都相同，将相同的值分组并依次赋值为1、2、3...
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # 将每个zone映射到对应的数字（1、2、3...）
            bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
        end
        bus[:,BASE_KV] = bus1.vn_kv 
        bus[:,VM] .= 1.0
        bus[:,VA] .= 0.0
        bus[:,VMAX] = bus1.max_vm_pu
        bus[:,VMIN] = bus1.min_vm_pu
        
    elseif sequence == 2
        # 负序网络
        # 处理母线数据
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # 如果所有元素都相同，统一赋值为1
            bus[:,ZONE] .= 1
        else
            # 如果元素不都相同，将相同的值分组并依次赋值为1、2、3...
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # 将每个zone映射到对应的数字（1、2、3...）
            bus[:,ZONE] = [zone_mapping[zone] for zone in bus1.zone]
        end
        bus[:,BASE_KV] = bus1.vn_kv 
        bus[:,VM].= 0.0
        bus[:,VA] .= 0.0
        bus[:,VMAX] = bus1.max_vm_pu
        bus[:,VMIN] = bus1.min_vm_pu
        
    else 
        # 零序网络
        # 处理母线数据
        bus1 = deepcopy(net["bus"])  
        bus[:, BUS_I] = bus1.index
        bus[:, BUS_TYPE] .= 1
        bus[slack_bus, BUS_TYPE] .= 3
        if all(x -> x == bus1.zone[1], bus1.zone)
            # 如果所有元素都相同，统一赋值为1
            bus[:,ZONE] .= 1
        else
            # 如果元素不都相同，将相同的值分组并依次赋值为1、2、3...
            unique_zones = unique(bus1.zone)
            zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
            
            # 将每个zone映射到对应的数字（1、2、3...）
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

function calculate_bus(net, jpc , slack_bus, opt)

    # 初始化
    nb = size(net["bus"], 1)
    bus = zeros(nb, 13)

    # 正序网络
    # 处理母线数据
    bus1 = deepcopy(net["bus"])  
    bus[:, BUS_I] = bus1.index
    bus[:, BUS_TYPE] .= 1
    bus[slack_bus, BUS_TYPE] .= 3
    if all(x -> x == bus1.zone[1], bus1.zone)
        # 如果所有元素都相同，统一赋值为1
        bus[:,ZONE] .= 1
    else
        # 如果元素不都相同，将相同的值分组并依次赋值为1、2、3...
        unique_zones = unique(bus1.zone)
        zone_mapping = Dict(zone => i for (i, zone) in enumerate(unique_zones))
        
        # 将每个zone映射到对应的数字（1、2、3...）
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

function add_grid_external_sc_impedance(jpc_new,external_grid)
    
    external_bus = external_grid.bus
    c=1.1
    s_sc = external_grid.s_sc_max_mva/jpc_new["baseMVA"]
    rx = external_grid.rx_max
    z_grid = c / (s_sc/3)
    x_grid = z_grid./sqrt.(1 .+rx.^2)
    r_grid = x_grid * rx

    Y_grid = 1 ./ (r_grid .+ 1im*x_grid)
    buses,gs,bs = sum_by_group(external_bus,real(Y_grid),imag(Y_grid))
    jpc_new["bus"][external_bus,GS] .= gs * jpc_new["baseMVA"]
    jpc_new["bus"][external_bus,BS] .= bs * jpc_new["baseMVA"]

    return gs * jpc_new["baseMVA"], bs * jpc_new["baseMVA"]
end

