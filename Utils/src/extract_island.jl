function extract_islands(jpc::JPC)
    
    # 设置连接矩阵
    nb = size(jpc.busAC, 1)     # 母线数量
    nl = size(jpc.branchAC, 1)  # 支路数量
    ng = size(jpc.genAC, 1)     # 发电机数量
    nld = size(jpc.loadAC, 1)   # 负载数量

    # 创建外部母线编号到内部索引的映射
    i2e = jpc.busAC[:, BUS_I]
    e2i = Dict{Int, Int}()
    for i in 1:nb
        e2i[Int(i2e[i])] = i
    end
    
    # 找出所有岛屿
    groups, isolated = find_islands(jpc)
    
    # 提取每个岛屿
    jpc_list = JPC[]
    
    # 创建一个集合，记录所有在有效岛屿中的母线编号
    all_valid_island_buses = Set{Int}()
    
    # 处理每个岛屿
    for i in eachindex(groups)
        # 获取岛屿i中的外部母线编号
        b_external = groups[i]
        
        # 检查岛屿是否有发电机或参考节点（即是否带电）
        has_power_source = false
        for bus_id in b_external
            # 找到对应的母线在jpc.busAC中的行
            bus_row = findfirst(x -> Int(x) == bus_id, jpc.busAC[:, BUS_I])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, BUS_TYPE])
                if bus_type == PV || bus_type == REF
                    has_power_source = true
                    break
                end
            end
            
            # 检查是否有连接到该母线的发电机
            for j in 1:ng
                gen_bus = Int(jpc.genAC[j, GEN_BUS])
                if gen_bus == bus_id && jpc.genAC[j, GEN_STATUS] == 1
                    has_power_source = true
                    break
                end
            end
            
            if has_power_source
                break
            end
        end
        
        # 如果岛屿没有电源，将其母线添加到isolated中，并跳过后续处理
        if !has_power_source
            append!(isolated, b_external)
            continue
        end
        
        # 记录有效岛屿中的母线
        union!(all_valid_island_buses, b_external)
        
        # 将外部母线编号转换为内部索引
        b_internal = Int[]
        for bus_id in b_external
            if haskey(e2i, bus_id)
                push!(b_internal, e2i[bus_id])
            else
                @warn "母线编号 $bus_id 不在系统中"
            end
        end
        
        # 找出两端都在岛屿i中的支路
        ibr = Int[]
        for j in 1:nl
            f_bus = Int(jpc.branchAC[j, F_BUS])
            t_bus = Int(jpc.branchAC[j, T_BUS])
            if (f_bus in b_external) && (t_bus in b_external)
                push!(ibr, j)
            end
        end
        
        # 找出连接到岛屿i中母线的发电机
        ig = Int[]
        for j in 1:ng
            gen_bus = Int(jpc.genAC[j, GEN_BUS])
            if gen_bus in b_external
                push!(ig, j)
            end
        end
        
        # 找出连接到岛屿i中母线的负载
        ild = Int[]
        for j in 1:nld
            load_bus = Int(jpc.loadAC[j, LOAD_CND])  # 使用LOAD_CND常量表示负载连接的母线
            if load_bus in b_external
                push!(ild, j)
            end
        end
        
        # 找出连接到岛屿i中母线的灵活负载
        ild_flex = Int[]
        if size(jpc.loadAC_flex, 1) > 0
            for j in 1:size(jpc.loadAC_flex, 1)
                load_bus = Int(jpc.loadAC_flex[j, 1])  # 假设第一列是母线编号
                if load_bus in b_external
                    push!(ild_flex, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的非对称负载
        ild_asymm = Int[]
        if size(jpc.loadAC_asymm, 1) > 0
            for j in 1:size(jpc.loadAC_asymm, 1)
                load_bus = Int(jpc.loadAC_asymm[j, 1])  # 假设第一列是母线编号
                if load_bus in b_external
                    push!(ild_asymm, j)
                end
            end
        end
        
        # 找出两端都在岛屿i中的三相支路
        ibr3ph = Int[]
        if size(jpc.branch3ph, 1) > 0
            for j in 1:size(jpc.branch3ph, 1)
                f_bus = Int(jpc.branch3ph[j, 1])  # 假设第一列是起始母线
                t_bus = Int(jpc.branch3ph[j, 2])  # 假设第二列是终止母线
                if (f_bus in b_external) && (t_bus in b_external)
                    push!(ibr3ph, j)
                end
            end
        end
        
        # 找出岛屿i中的DC母线
        bdc_external = Int[]
        bdc_internal = Int[]
        if size(jpc.busDC, 1) > 0
            for j in 1:size(jpc.busDC, 1)
                bus_id = Int(jpc.busDC[j, 1])  # 假设第一列是DC母线编号
                # 这里需要确定DC母线与AC母线的关联关系
                if bus_id in b_external
                    push!(bdc_external, bus_id)
                    push!(bdc_internal, j)
                end
            end
        end
        
        # 找出两端都在岛屿i中的DC支路
        ibrdc = Int[]
        if size(jpc.branchDC, 1) > 0
            for j in 1:size(jpc.branchDC, 1)
                f_bus = Int(jpc.branchDC[j, 1])  # 假设第一列是起始母线
                t_bus = Int(jpc.branchDC[j, 2])  # 假设第二列是终止母线
                if (f_bus in bdc_external) && (t_bus in bdc_external)
                    push!(ibrdc, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的AC分布式发电
        isgen = Int[]
        if size(jpc.sgenAC, 1) > 0
            for j in 1:size(jpc.sgenAC, 1)
                sgen_bus = Int(jpc.sgenAC[j, 1])  # 假设第一列是母线编号
                if sgen_bus in b_external
                    push!(isgen, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的储能系统
        istorage = Int[]
        if size(jpc.storage, 1) > 0
            for j in 1:size(jpc.storage, 1)
                storage_bus = Int(jpc.storage[j, 1])  # 假设第一列是母线编号
                if storage_bus in b_external
                    push!(istorage, j)
                end
            end
        end
        
        # 找出连接到岛屿i中DC母线的DC分布式发电
        isgendc = Int[]
        if size(jpc.sgenDC, 1) > 0
            for j in 1:size(jpc.sgenDC, 1)
                sgen_bus = Int(jpc.sgenDC[j, 1])  # 假设第一列是DC母线编号
                if sgen_bus in bdc_external
                    push!(isgendc, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的转换器
        iconv = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])  # 假设第一列是AC母线编号
                dc_bus = Int(jpc.converter[j, 2])  # 假设第二列是DC母线编号
                if (ac_bus in b_external) && (dc_bus in bdc_external)
                    push!(iconv, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的外部电网
        iext = Int[]
        if size(jpc.ext_grid, 1) > 0
            for j in 1:size(jpc.ext_grid, 1)
                ext_bus = Int(jpc.ext_grid[j, 1])  # 假设第一列是母线编号
                if ext_bus in b_external
                    push!(iext, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的高压断路器
        ihvcb = Int[]
        if size(jpc.hvcb, 1) > 0
            for j in 1:size(jpc.hvcb, 1)
                f_bus = Int(jpc.hvcb[j, 1])  # 假设第一列是起始母线
                t_bus = Int(jpc.hvcb[j, 2])  # 假设第二列是终止母线
                if (f_bus in b_external) && (t_bus in b_external)
                    push!(ihvcb, j)
                end
            end
        end
        
        # 找出岛屿i中的微电网
        img = Int[]
        if size(jpc.microgrid, 1) > 0
            for j in 1:size(jpc.microgrid, 1)
                mg_bus = Int(jpc.microgrid[j, 1])  # 假设第一列是母线编号
                if mg_bus in b_external
                    push!(img, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的光伏系统
        ipv = Int[]
        if isdefined(jpc, :pv) && size(jpc.pv, 1) > 0
            for j in 1:size(jpc.pv, 1)
                pv_bus = Int(jpc.pv[j, 2])  # 假设第二列是母线编号
                if pv_bus in b_external
                    # 如果有状态字段，检查设备是否在运行
                    is_in_service = size(jpc.pv, 2) >= 9 ? jpc.pv[j, 9] == 1 : true
                    
                    if is_in_service
                        push!(ipv, j)
                    end
                end
            end
        end
        
        # 找出连接到岛屿i中母线的交流光伏系统
        ipv_ac = Int[]
        if isdefined(jpc, :pv_acsystem) && size(jpc.pv_acsystem, 1) > 0
            for j in 1:size(jpc.pv_acsystem, 1)
                pv_bus = Int(jpc.pv_acsystem[j, 2])  # 假设第二列是母线编号
                if pv_bus in b_external
                    # 如果有状态字段(假设是第9列)，检查设备是否在运行
                    is_in_service = size(jpc.pv_acsystem, 2) >= 9 ? jpc.pv_acsystem[j, 9] == 1 : true
                    
                    if is_in_service
                        push!(ipv_ac, j)
                    end
                end
            end
        end
        
        # 创建这个岛屿的JPC副本
        jpck = JPC(jpc.version, jpc.baseMVA)
        
        # 复制相关数据
        if !isempty(b_internal)
            jpck.busAC = jpc.busAC[b_internal, :]
        end
        
        if !isempty(ibr)
            jpck.branchAC = jpc.branchAC[ibr, :]
        end
        
        if !isempty(ig)
            jpck.genAC = jpc.genAC[ig, :]
        end
        
        if !isempty(ild)
            jpck.loadAC = jpc.loadAC[ild, :]
        end
        
        if !isempty(ild_flex)
            jpck.loadAC_flex = jpc.loadAC_flex[ild_flex, :]
        end
        
        if !isempty(ild_asymm)
            jpck.loadAC_asymm = jpc.loadAC_asymm[ild_asymm, :]
        end
        
        if !isempty(ibr3ph)
            jpck.branch3ph = jpc.branch3ph[ibr3ph, :]
        end
        
        if !isempty(bdc_internal)
            jpck.busDC = jpc.busDC[bdc_internal, :]
        end
        
        if !isempty(ibrdc)
            jpck.branchDC = jpc.branchDC[ibrdc, :]
        end
        
        if !isempty(isgen)
            jpck.sgenAC = jpc.sgenAC[isgen, :]
        end
        
        if !isempty(istorage)
            jpck.storage = jpc.storage[istorage, :]
        end
        
        if !isempty(isgendc)
            jpck.sgenDC = jpc.sgenDC[isgendc, :]
        end
        
        if !isempty(iconv)
            jpck.converter = jpc.converter[iconv, :]
        end
        
        if !isempty(iext)
            jpck.ext_grid = jpc.ext_grid[iext, :]
        end
        
        if !isempty(ihvcb)
            jpck.hvcb = jpc.hvcb[ihvcb, :]
        end
        
        if !isempty(img)
            jpck.microgrid = jpc.microgrid[img, :]
        end
        
        # 复制光伏数据
        if !isempty(ipv) && isa(jpc.pv, Array)
            jpck.pv = jpc.pv[ipv, :]
        else
            # 如果没有光伏设备连接到这个岛屿，创建一个空数组
            if isa(jpc.pv, Array)
                jpck.pv = similar(jpc.pv, 0, size(jpc.pv, 2))
            else
                jpck.pv = deepcopy(jpc.pv)  # 如果pv不是数组，直接复制
            end
        end
        
        # 复制交流光伏系统数据
        if !isempty(ipv_ac) && isdefined(jpc, :pv_acsystem)
            jpck.pv_acsystem = jpc.pv_acsystem[ipv_ac, :]
        else
            # 如果没有交流光伏系统连接到这个岛屿，创建一个空数组
            if isdefined(jpc, :pv_acsystem) && isa(jpc.pv_acsystem, Array)
                jpck.pv_acsystem = similar(jpc.pv_acsystem, 0, size(jpc.pv_acsystem, 2))
            elseif isdefined(jpc, :pv_acsystem)
                jpck.pv_acsystem = deepcopy(jpc.pv_acsystem)  # 如果pv_acsystem不是数组，直接复制
            end
        end
        
        push!(jpc_list, jpck)
    end
    
    # 检查是否有不在任何有效岛屿中的母线（即不在all_valid_island_buses中且不在isolated中的母线）
    for i in 1:nb
        bus_id = Int(jpc.busAC[i, BUS_I])
        if !(bus_id in all_valid_island_buses) && !(bus_id in isolated)
            push!(isolated, bus_id)
        end
    end
    
    return jpc_list, isolated
end


function extract_islands_acdc(jpc::JPC)
    # 设置连接矩阵
    nb_ac = size(jpc.busAC, 1)     # AC母线数量
    nb_dc = size(jpc.busDC, 1)     # DC母线数量
    nl_ac = size(jpc.branchAC, 1)  # AC支路数量
    nl_dc = size(jpc.branchDC, 1)  # DC支路数量
    ng = size(jpc.genAC, 1)        # 发电机数量
    nld = size(jpc.loadAC, 1)      # 负载数量
    
    # 创建外部母线编号到内部索引的映射
    i2e_ac = jpc.busAC[:, BUS_I]
    e2i_ac = Dict{Int, Int}()
    for i in 1:nb_ac
        e2i_ac[Int(i2e_ac[i])] = i
    end
    
    i2e_dc = nb_dc > 0 ? jpc.busDC[:, 1] : Int[]  # 假设第一列是DC母线编号
    e2i_dc = Dict{Int, Int}()
    for i in 1:nb_dc
        e2i_dc[Int(i2e_dc[i])] = i
    end
    
    # 找出所有岛屿
    groups, isolated = find_islands_acdc(jpc)
    
    # 提取每个岛屿
    jpc_list = JPC[]
    
    # 创建一个集合，记录所有在有效岛屿中的母线编号
    all_valid_island_buses = Set{Int}()
    
    # 处理每个岛屿
    for i in eachindex(groups)
        # 获取岛屿i中的外部AC母线编号
        b_external_ac = groups[i]
        
        # 检查岛屿是否有发电机或参考节点（即是否带电）
        has_power_source = false
        for bus_id in b_external_ac
            # 找到对应的母线在jpc.busAC中的行
            bus_row = findfirst(x -> Int(x) == bus_id, jpc.busAC[:, BUS_I])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, BUS_TYPE])
                if bus_type == PV || bus_type == REF
                    has_power_source = true
                    break
                end
            end
            
            # 检查是否有连接到该母线的发电机
            for j in 1:ng
                gen_bus = Int(jpc.genAC[j, GEN_BUS])
                if gen_bus == bus_id && jpc.genAC[j, GEN_STATUS] == 1
                    has_power_source = true
                    break
                end
            end
            
            if has_power_source
                break
            end
        end
        
        # 找出与这些AC母线相连的DC母线
        b_external_dc = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])
                dc_bus = Int(jpc.converter[j, 2])
                # 假设第3列是状态字段，如果没有状态字段，则假设为投运
                is_active = size(jpc.converter, 2) >= 3 ? jpc.converter[j, 3] == 1 : true
                
                if ac_bus in b_external_ac && is_active
                    push!(b_external_dc, dc_bus)
                end
            end
        end
        
        # 通过DC支路扩展DC母线集合
        if !isempty(b_external_dc) && nl_dc > 0
            # 使用广度优先搜索找到所有连通的DC母线
            visited_dc = Set(b_external_dc)
            queue = copy(b_external_dc)
            
            while !isempty(queue)
                current_bus = popfirst!(queue)
                
                for j in 1:nl_dc
                    branch = jpc.branchDC[j, :]
                    # 假设DC支路结构与AC支路相似
                    if branch[11] != 1  # 检查状态
                        continue
                    end
                    
                    f_bus = Int(branch[1])
                    t_bus = Int(branch[2])
                    
                    if f_bus == current_bus && !(t_bus in visited_dc)
                        push!(visited_dc, t_bus)
                        push!(queue, t_bus)
                    elseif t_bus == current_bus && !(f_bus in visited_dc)
                        push!(visited_dc, f_bus)
                        push!(queue, f_bus)
                    end
                end
            end
            
            b_external_dc = collect(visited_dc)
        end
        
        # 检查DC系统是否有电源
        if !has_power_source && !isempty(b_external_dc)
            # 检查DC母线中是否有参考节点（类型为2）
            for bus_id in b_external_dc
                # 找到对应的母线在jpc.busDC中的行
                bus_row = findfirst(x -> Int(x) == bus_id, jpc.busDC[:, 1])
                if bus_row !== nothing && size(jpc.busDC, 2) >= 2
                    bus_type = Int(jpc.busDC[bus_row, 2])  # 假设第2列是母线类型
                    if bus_type == 2  # DC参考节点类型为2
                        has_power_source = true
                        break
                    end
                end
            end
            
            # 检查DC分布式发电
            if !has_power_source && size(jpc.sgenDC, 1) > 0
                for j in 1:size(jpc.sgenDC, 1)
                    sgen_bus = Int(jpc.sgenDC[j, 1])
                    # 假设第3列是状态字段
                    is_active = size(jpc.sgenDC, 2) >= 3 ? jpc.sgenDC[j, 3] == 1 : true
                    
                    if sgen_bus in b_external_dc && is_active
                        has_power_source = true
                        break
                    end
                end
            end
            
            # 检查DC发电机
            if !has_power_source && isdefined(jpc, :genDC) && size(jpc.genDC, 1) > 0
                for j in 1:size(jpc.genDC, 1)
                    gen_bus = Int(jpc.genDC[j, 1])
                    is_active = size(jpc.genDC, 2) >= 3 ? jpc.genDC[j, 3] == 1 : true
                    
                    if gen_bus in b_external_dc && is_active
                        has_power_source = true
                        break
                    end
                end
            end
        end
        
        # 如果岛屿没有电源，将其母线添加到isolated中，并跳过后续处理
        if !has_power_source
            append!(isolated, b_external_ac)
            continue
        end
        
        # 记录有效岛屿中的母线
        union!(all_valid_island_buses, b_external_ac)
        
        # 将外部母线编号转换为内部索引
        b_internal_ac = Int[]
        for bus_id in b_external_ac
            if haskey(e2i_ac, bus_id)
                push!(b_internal_ac, e2i_ac[bus_id])
            else
                @warn "AC母线编号 $bus_id 不在系统中"
            end
        end
        
        b_internal_dc = Int[]
        for bus_id in b_external_dc
            if haskey(e2i_dc, bus_id)
                push!(b_internal_dc, e2i_dc[bus_id])
            else
                @warn "DC母线编号 $bus_id 不在系统中"
            end
        end
        
        # 找出两端都在岛屿i中的AC支路
        ibr_ac = Int[]
        for j in 1:nl_ac
            f_bus = Int(jpc.branchAC[j, F_BUS])
            t_bus = Int(jpc.branchAC[j, T_BUS])
            if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                push!(ibr_ac, j)
            end
        end
        
        # 找出两端都在岛屿i中的DC支路
        ibr_dc = Int[]
        for j in 1:nl_dc
            f_bus = Int(jpc.branchDC[j, 1])  # 假设第一列是起始母线
            t_bus = Int(jpc.branchDC[j, 2])  # 假设第二列是终止母线
            if (f_bus in b_external_dc) && (t_bus in b_external_dc)
                push!(ibr_dc, j)
            end
        end
        
        # 找出连接岛屿i中AC和DC母线的转换器
        iconv = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])
                dc_bus = Int(jpc.converter[j, 2])
                if (ac_bus in b_external_ac) && (dc_bus in b_external_dc)
                    push!(iconv, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的发电机
        ig = Int[]
        for j in 1:ng
            gen_bus = Int(jpc.genAC[j, GEN_BUS])
            if gen_bus in b_external_ac
                push!(ig, j)
            end
        end
        
        # 找出连接到岛屿i中DC母线的DC发电机
        igen_dc = Int[]
        if isdefined(jpc, :genDC) && size(jpc.genDC, 1) > 0
            for j in 1:size(jpc.genDC, 1)
                gen_bus = Int(jpc.genDC[j, 1])  # 假设第一列是母线编号
                if gen_bus in b_external_dc
                    push!(igen_dc, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的负载
        ild = Int[]
        for j in 1:nld
            load_bus = Int(jpc.loadAC[j, LOAD_CND])  # 使用LOAD_CND常量表示负载连接的母线
            if load_bus in b_external_ac
                push!(ild, j)
            end
        end
        
        # 找出连接到岛屿i中母线的光伏设备
        ipv = Int[]
        if isa(jpc.pv, Array) && size(jpc.pv, 1) > 0
            for j in 1:size(jpc.pv, 1)
                pv_bus = Int(jpc.pv[j, 2])  # 使用第二列作为母线编号
                
                # 检查光伏设备是否连接到直流母线
                if pv_bus in b_external_dc
                    # 如果有状态字段(假设是第8列或其他位置)，检查设备是否在运行
                    is_in_service = size(jpc.pv, 2) >= 8 ? jpc.pv[j, 8] == 1 : true
                    
                    if is_in_service
                        push!(ipv, j)
                    end
                end
            end
        end
        
        # 找出连接到岛屿i中母线的交流光伏系统
        ipv_ac = Int[]
        if isdefined(jpc, :pv_acsystem) && size(jpc.pv_acsystem, 1) > 0
            for j in 1:size(jpc.pv_acsystem, 1)
                pv_bus = Int(jpc.pv_acsystem[j, PV_AC_BUS])  # 使用PV_AC_BUS常量获取母线编号
                if pv_bus in b_external_ac
                    # 如果有状态字段，检查设备是否在运行
                    is_in_service = jpc.pv_acsystem[j, PV_AC_IN_SERVICE] == 1
                    
                    if is_in_service
                        push!(ipv_ac, j)
                    end
                end
            end
        end
        
        # 找出连接到岛屿i中母线的灵活负载
        ild_flex = Int[]
        if size(jpc.loadAC_flex, 1) > 0
            for j in 1:size(jpc.loadAC_flex, 1)
                load_bus = Int(jpc.loadAC_flex[j, 1])  # 假设第一列是母线编号
                if load_bus in b_external_ac
                    push!(ild_flex, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的非对称负载
        ild_asymm = Int[]
        if size(jpc.loadAC_asymm, 1) > 0
            for j in 1:size(jpc.loadAC_asymm, 1)
                load_bus = Int(jpc.loadAC_asymm[j, 1])  # 假设第一列是母线编号
                if load_bus in b_external_ac
                    push!(ild_asymm, j)
                end
            end
        end
        
        # 找出两端都在岛屿i中的三相支路
        ibr3ph = Int[]
        if size(jpc.branch3ph, 1) > 0
            for j in 1:size(jpc.branch3ph, 1)
                f_bus = Int(jpc.branch3ph[j, 1])  # 假设第一列是起始母线
                t_bus = Int(jpc.branch3ph[j, 2])  # 假设第二列是终止母线
                if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                    push!(ibr3ph, j)
                end
            end
        end
        
        # 找出连接到岛屿i中DC母线的DC负载
        ild_dc = Int[]
        if size(jpc.loadDC, 1) > 0
            for j in 1:size(jpc.loadDC, 1)
                load_bus = Int(jpc.loadDC[j, 2])  # 假设第二列是母线编号
                if load_bus in b_external_dc
                    push!(ild_dc, j)
                end
            end
        end
        
        # 找出连接到岛屿i中DC母线的DC分布式发电
        isgen_dc = Int[]
        if size(jpc.sgenDC, 1) > 0
            for j in 1:size(jpc.sgenDC, 1)
                sgen_bus = Int(jpc.sgenDC[j, 1])  # 假设第一列是母线编号
                if sgen_bus in b_external_dc
                    push!(isgen_dc, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的AC分布式发电
        isgen_ac = Int[]
        if size(jpc.sgenAC, 1) > 0
            for j in 1:size(jpc.sgenAC, 1)
                sgen_bus = Int(jpc.sgenAC[j, 1])  # 假设第一列是母线编号
                if sgen_bus in b_external_ac
                    push!(isgen_ac, j)
                end
            end
        end
        
        # 找出连接到岛屿i中DC母线的储能系统
        istorage = Int[]
        if size(jpc.storage, 1) > 0
            # 确定储能系统连接的DC母线编号所在的列索引
            ess_bus_col = -1
            # 尝试查找名为ESS_BUS的列（如果存在列名）
            if isdefined(jpc.storage, :colnames) && :ESS_BUS in jpc.storage.colnames
                ess_bus_col = findfirst(x -> x == :ESS_BUS, jpc.storage.colnames)
            else
                # 如果没有列名，假设第1列是储能系统连接的DC母线编号
                ess_bus_col = 1
            end
            
            if ess_bus_col > 0 && ess_bus_col <= size(jpc.storage, 2)
                for j in 1:size(jpc.storage, 1)
                    storage_bus = Int(jpc.storage[j, ess_bus_col])
                    if storage_bus in b_external_dc
                        push!(istorage, j)
                    end
                end
            end
        end
        
        # 找出连接到岛屿i中母线的外部电网
        iext = Int[]
        if size(jpc.ext_grid, 1) > 0
            for j in 1:size(jpc.ext_grid, 1)
                ext_bus = Int(jpc.ext_grid[j, 1])  # 假设第一列是母线编号
                if ext_bus in b_external_ac
                    push!(iext, j)
                end
            end
        end
        
        # 找出连接到岛屿i中母线的高压断路器
        ihvcb = Int[]
        if size(jpc.hvcb, 1) > 0
            for j in 1:size(jpc.hvcb, 1)
                f_bus = Int(jpc.hvcb[j, 1])  # 假设第一列是起始母线
                t_bus = Int(jpc.hvcb[j, 2])  # 假设第二列是终止母线
                if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                    push!(ihvcb, j)
                end
            end
        end
        
        # 找出岛屿i中的微电网
        img = Int[]
        if size(jpc.microgrid, 1) > 0
            for j in 1:size(jpc.microgrid, 1)
                mg_bus = Int(jpc.microgrid[j, 1])  # 假设第一列是母线编号
                if mg_bus in b_external_ac
                    push!(img, j)
                end
            end
        end
        
        # 创建这个岛屿的JPC副本
        jpck = JPC(jpc.version, jpc.baseMVA)
        
        # 复制相关数据
        if !isempty(b_internal_ac)
            jpck.busAC = deepcopy(jpc.busAC[b_internal_ac, :])
        end
        
        if !isempty(b_internal_dc)
            jpck.busDC = deepcopy(jpc.busDC[b_internal_dc, :])
            
            # 确保DC系统中的参考节点类型保持为2而不是3
            if size(jpck.busDC, 1) > 0 && size(jpck.busDC, 2) >= 2
                for j in 1:size(jpck.busDC, 1)
                    if jpck.busDC[j, 2] == 3  # 如果被错误地设置为3
                        jpck.busDC[j, 2] = 2  # 将其改回2
                    end
                end
            end
        end
        
        if !isempty(ibr_ac)
            jpck.branchAC = deepcopy(jpc.branchAC[ibr_ac, :])
        end
        
        if !isempty(ibr_dc)
            jpck.branchDC = deepcopy(jpc.branchDC[ibr_dc, :])
        end
        
        if !isempty(iconv)
            jpck.converter = deepcopy(jpc.converter[iconv, :])
        end
        
        if !isempty(ig)
            jpck.genAC = deepcopy(jpc.genAC[ig, :])
        end
        
        if !isempty(igen_dc) && isdefined(jpc, :genDC)
            jpck.genDC = deepcopy(jpc.genDC[igen_dc, :])
        end
        
        if !isempty(ild)
            jpck.loadAC = deepcopy(jpc.loadAC[ild, :])
        end
        
        if !isempty(ild_dc)
            jpck.loadDC = deepcopy(jpc.loadDC[ild_dc, :])
        end
        
        # 复制光伏数据
        if !isempty(ipv) && isa(jpc.pv, Array)
            jpck.pv = deepcopy(jpc.pv[ipv, :])
        else
            # 如果没有光伏设备连接到这个岛屿，创建一个空数组
            if isa(jpc.pv, Array)
                jpck.pv = similar(jpc.pv, 0, size(jpc.pv, 2))
            else
                jpck.pv = deepcopy(jpc.pv)  # 如果pv不是数组，直接复制
            end
        end
        
        # 复制交流光伏系统数据
        if !isempty(ipv_ac) && isdefined(jpc, :pv_acsystem)
            jpck.pv_acsystem = deepcopy(jpc.pv_acsystem[ipv_ac, :])
        else
            # 如果没有交流光伏系统连接到这个岛屿，创建一个空数组
            if isdefined(jpc, :pv_acsystem) && isa(jpc.pv_acsystem, Array)
                jpck.pv_acsystem = similar(jpc.pv_acsystem, 0, size(jpc.pv_acsystem, 2))
            elseif isdefined(jpc, :pv_acsystem)
                jpck.pv_acsystem = deepcopy(jpc.pv_acsystem)  # 如果pv_acsystem不是数组，直接复制
            end
        end
        
        if !isempty(ild_flex)
            jpck.loadAC_flex = deepcopy(jpc.loadAC_flex[ild_flex, :])
        end
        
        if !isempty(ild_asymm)
            jpck.loadAC_asymm = deepcopy(jpc.loadAC_asymm[ild_asymm, :])
        end
        
        if !isempty(ibr3ph)
            jpck.branch3ph = deepcopy(jpc.branch3ph[ibr3ph, :])
        end
        
        if !isempty(isgen_ac)
            jpck.sgenAC = deepcopy(jpc.sgenAC[isgen_ac, :])
        end
        
        if !isempty(isgen_dc)
            jpck.sgenDC = deepcopy(jpc.sgenDC[isgen_dc, :])
        end
        
        if !isempty(istorage)
            jpck.storage = deepcopy(jpc.storage[istorage, :])
        end
        
        if !isempty(iext)
            jpck.ext_grid = deepcopy(jpc.ext_grid[iext, :])
        end
        
        if !isempty(ihvcb)
            jpck.hvcb = deepcopy(jpc.hvcb[ihvcb, :])
        end
        
        if !isempty(img)
            jpck.microgrid = deepcopy(jpc.microgrid[img, :])
        end
        
        push!(jpc_list, jpck)
    end
    
    # 检查是否有不在任何有效岛屿中的母线（即不在all_valid_island_buses中且不在isolated中的母线）
    for i in 1:nb_ac
        bus_id = Int(jpc.busAC[i, BUS_I])
        if !(bus_id in all_valid_island_buses) && !(bus_id in isolated)
            push!(isolated, bus_id)
        end
    end
    
    # 检测是否在Vdc恒定模式下存在蓄电池
    for i in 1:length(jpc_list)
        if !isempty(jpc_list[i].storage)
            # 使用逻辑或的向量化操作 .||
            Vdc_converter = findall((jpc_list[i].converter[:,CONV_MODE] .== 4) .|| (jpc_list[i].converter[:,CONV_MODE] .== 5))
            if !isempty(Vdc_converter)
                @error "在岛屿 $(i) 中发现Vdc恒定模式下的蓄电池，可能导致功率流计算错误。请检查JPC数据。"
            end
        end

        slack_bus_indices = findall(jpc_list[i].busAC[:, BUS_TYPE] .== REF)
        if length(slack_bus_indices) > 1
            # 检查是否有多个参考节点
                @warn "在岛屿 $(i) 中发现多个参考节点，可能导致功率流计算错误。请检查JPC数据。"
        end
    end

    return jpc_list, isolated
end


