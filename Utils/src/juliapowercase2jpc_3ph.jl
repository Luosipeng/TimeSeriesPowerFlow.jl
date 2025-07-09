function JuliaPowerCase2Jpc_3ph(case::JuliaPowerCase)
    # 1. 合并虚拟节点
    case = PowerFlow.merge_virtual_nodes(case)

    # 2. 创建JPC_3ph对象
    jpc_3ph = PowerFlow.JPC_3ph()

    # 3. 设置基本参数
    jpc_3ph.baseMVA = case.baseMVA
    jpc_3ph.basef = case.basef

    # 4. 设置节点数据
    JPC_3ph_buses_process(case, jpc_3ph)

    # 5. 设置线路数据
    JPC_3ph_branches_process(case, jpc_3ph)
    # 6. 设置发电机数据
    JPC_3ph_gens_process(case, jpc_3ph)
    # 7. 设置负荷数据
    JPC_3ph_loads_process(case, jpc_3ph)

    
    jpc_3ph, _, _ = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 1)
    jpc_3ph, gs_eg, bs_eg = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 2)
    jpc_3ph, _ , _ = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 0)

    return jpc_3ph, gs_eg, bs_eg

end

function JPC_3ph_buses_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 获取节点数据并深拷贝防止误操作
    buses = deepcopy(case.busesAC)
    
    # 创建一个空矩阵，行数为节点数，列数为13
    num_buses = length(buses)
    bus_matrix = zeros(num_buses, 13)
    
    for (i, bus) in enumerate(buses)
        # 设置电压初始值（根据序号）
        vm = 0.0
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
    jpc_3ph.busAC_0 = bus_matrix
    jpc_3ph.busAC_2 = bus_matrix

    bus_matrix_1 = zeros(num_buses, 13)
    for (i, bus) in enumerate(buses)
        # 设置电压初始值（根据序号）
        vm = 1.0
        va = 0.0
        
        # 填充矩阵的每一行
        bus_matrix_1[i, :] = [
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
    jpc_3ph.busAC_1 = bus_matrix_1
    
    return jpc_3ph
end

function JPC_3ph_branches_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 计算线路参数
    calculate_3ph_line_parameters(case, jpc_3ph)
    # 计算变压器参数
    calculate_3ph_transformer2w_parameters(case, jpc_3ph)
    # 计算零序支路参数
    jpc_3ph = calculate_branch_JPC_zero(case, jpc_3ph)

    return jpc_3ph
end

function calculate_3ph_line_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 处理线路数据，转换为JPC格式
    nbr = length(case.branchesAC)
    branch_1 = zeros(nbr, 14)  # 正序分量
    branch_2 = zeros(nbr, 14)  # 负序分量
    lines = case.branchesAC

    for (i, line) in enumerate(lines)
        # 获取起始和终止母线编号
        from_bus_idx = line.from_bus
        to_bus_idx = line.to_bus
        
        # 获取起始母线的基准电压(kV)
        basekv = jpc_3ph.busAC_1[from_bus_idx, BASE_KV]  # 使用三相系统的母线数据
        
        # 计算基准阻抗
        baseR = (basekv^2) / case.baseMVA
        
        # 考虑并联线路的情况
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # 计算正序阻抗参数
        r1_pu = line.length_km * line.r_ohm_per_km / baseR / parallel
        x1_pu = line.length_km * line.x_ohm_per_km / baseR / parallel
        
        # 计算负序阻抗参数（通常与正序相同或略有不同）
        r2_pu = r1_pu  # 对于大多数线路，负序电阻等于正序电阻
        x2_pu = x1_pu  # 对于大多数线路，负序电抗等于正序电抗
        
        # 如果线路有特定的负序参数，可以使用以下代码：
        if hasfield(typeof(line), :r2_ohm_per_km)
            r2_pu = line.length_km * line.r2_ohm_per_km / baseR / parallel
        end
        if hasfield(typeof(line), :x2_ohm_per_km)
            x2_pu = line.length_km * line.x2_ohm_per_km / baseR / parallel
        end
        
        # 计算并联电纳(p.u.) - 正序和负序通常相同
        b_pu = 2 * π * case.basef * line.length_km * line.c_nf_per_km * 1e-9 * baseR * parallel
        
        # 计算并联电导(p.u.)
        g_pu = 0.0
        if hasfield(typeof(line), :g_us_per_km)
            g_pu = line.g_us_per_km * 1e-6 * baseR * line.length_km * parallel
        end
        
        # 填充正序分量矩阵 (branchAC_1)
        branch_1[i, F_BUS] = from_bus_idx
        branch_1[i, T_BUS] = to_bus_idx
        branch_1[i, BR_R] = r1_pu
        branch_1[i, BR_X] = x1_pu
        branch_1[i, BR_B] = b_pu
        
        # 填充负序分量矩阵 (branchAC_2)
        branch_2[i, F_BUS] = from_bus_idx
        branch_2[i, T_BUS] = to_bus_idx
        branch_2[i, BR_R] = r2_pu
        branch_2[i, BR_X] = x2_pu
        branch_2[i, BR_B] = b_pu
        
        # 设置额定容量（正序和负序相同）
        rate_a = 100.0  # 默认值
        if hasfield(typeof(line), :max_i_ka)
            rate_a = line.max_i_ka * basekv * sqrt(3)  # 额定容量(MVA)
        end
        branch_1[i, RATE_A] = rate_a
        branch_2[i, RATE_A] = rate_a
        
        # 设置支路状态（正序和负序相同）
        status = line.in_service ? 1.0 : 0.0
        branch_1[i, BR_STATUS] = status
        branch_2[i, BR_STATUS] = status
        
        # 设置相角限制（正序和负序相同）
        branch_1[i, ANGMIN] = -360.0
        branch_1[i, ANGMAX] = 360.0
        branch_2[i, ANGMIN] = -360.0
        branch_2[i, ANGMAX] = 360.0
    end

    # 为三相系统赋值
    jpc_3ph.branchAC_1 = branch_1  # 正序分量
    jpc_3ph.branchAC_2 = branch_2  # 负序分量
end

function calculate_3ph_transformer2w_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 处理变压器数据，转换为JPC格式
    transformers = case.transformers_2w_etap
    nbr = length(transformers)
    
    if nbr == 0
        return  # 如果没有变压器，直接返回
    end
    
    # 创建变压器分支矩阵 - 正序和负序
    branch_1 = zeros(nbr, 14)  # 正序分量
    branch_2 = zeros(nbr, 14)  # 负序分量
    
    for (i, transformer) in enumerate(transformers)
        # 获取高压侧和低压侧母线编号
        hv_bus_idx = transformer.hv_bus
        lv_bus_idx = transformer.lv_bus
        
        # 获取高压侧母线的基准电压(kV)
        hv_basekv = jpc_3ph.busAC_1[hv_bus_idx, BASE_KV]  # 使用三相系统的母线数据
        
        # 计算阻抗参数
        # 变压器阻抗百分比转换为标幺值
        z_percent = transformer.z_percent
        x_r_ratio = transformer.x_r
        
        # 计算电阻和电抗（考虑基准功率转换）
        s_ratio = transformer.sn_mva / case.baseMVA
        z_pu = z_percent / 100.0 / s_ratio  # 转换到系统基准（百分比转小数）
        
        # 计算正序阻抗参数
        r1_pu = z_pu / sqrt(1 + x_r_ratio^2)
        x1_pu = r1_pu * x_r_ratio
        
        # 计算负序阻抗参数
        # 对于变压器，负序阻抗通常等于正序阻抗
        r2_pu = r1_pu
        x2_pu = x1_pu
        
        # 如果变压器有特定的负序参数，可以使用以下代码：
        if hasfield(typeof(transformer), :z2_percent)
            z2_pu = transformer.z2_percent / 100.0 / s_ratio
            if hasfield(typeof(transformer), :x2_r)
                x2_r_ratio = transformer.x2_r
            else
                x2_r_ratio = x_r_ratio  # 使用正序的X/R比
            end
            r2_pu = z2_pu / sqrt(1 + x2_r_ratio^2)
            x2_pu = r2_pu * x2_r_ratio
        end
        
        # 考虑并联变压器
        parallel = transformer.parallel
        if parallel > 1
            r1_pu = r1_pu / parallel
            x1_pu = x1_pu / parallel
            r2_pu = r2_pu / parallel
            x2_pu = x2_pu / parallel
        end
        
        # 填充正序分支矩阵 (branchAC_1)
        branch_1[i, F_BUS] = hv_bus_idx
        branch_1[i, T_BUS] = lv_bus_idx
        branch_1[i, BR_R] = r1_pu
        branch_1[i, BR_X] = x1_pu
        branch_1[i, BR_B] = 0.0  # 变压器通常没有并联电纳
        
        # 填充负序分支矩阵 (branchAC_2)
        branch_2[i, F_BUS] = hv_bus_idx
        branch_2[i, T_BUS] = lv_bus_idx
        branch_2[i, BR_R] = r2_pu
        branch_2[i, BR_X] = x2_pu
        branch_2[i, BR_B] = 0.0  # 变压器通常没有并联电纳
        
        # 设置变比和相移（正序和负序相同）
        tap_ratio = 1.0
        shift_angle = 0.0
        
        # 如果变压器有变比信息
        if hasfield(typeof(transformer), :tap_ratio)
            tap_ratio = transformer.tap_ratio
        end
        if hasfield(typeof(transformer), :shift_angle)
            shift_angle = transformer.shift_angle
        end
        
        branch_1[i, TAP] = tap_ratio
        branch_1[i, SHIFT] = shift_angle
        branch_2[i, TAP] = tap_ratio
        branch_2[i, SHIFT] = shift_angle
        
        # 设置额定容量（正序和负序相同）
        rate_a = transformer.sn_mva
        branch_1[i, RATE_A] = rate_a
        branch_2[i, RATE_A] = rate_a
        
        # 设置支路状态（正序和负序相同）
        status = transformer.in_service ? 1.0 : 0.0
        branch_1[i, BR_STATUS] = status
        branch_2[i, BR_STATUS] = status
        
        # 设置相角限制（正序和负序相同）
        branch_1[i, ANGMIN] = -360.0
        branch_1[i, ANGMAX] = 360.0
        branch_2[i, ANGMIN] = -360.0
        branch_2[i, ANGMAX] = 360.0
    end
    
    # 将变压器分支数据添加到三相JPC结构体
    if isempty(jpc_3ph.branchAC_1)
        jpc_3ph.branchAC_1 = branch_1
        jpc_3ph.branchAC_2 = branch_2
    else
        jpc_3ph.branchAC_1 = [jpc_3ph.branchAC_1; branch_1]
        jpc_3ph.branchAC_2 = [jpc_3ph.branchAC_2; branch_2]
    end
end

function calculate_branch_JPC_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 初始化分支数据矩阵
    nbr = length(case.branchesAC)
    branch = zeros(nbr, 14)
    
    # 设置默认值
    branch[:, RATE_A] .= 250.0
    branch[:, RATE_B] .= 250.0
    branch[:, RATE_C] .= 250.0
    branch[:, TAP] .= 1.0
    branch[:, SHIFT] .= 0.0
    branch[:, BR_STATUS] .= 1.0
    branch[:, ANGMIN] .= -360.0
    branch[:, ANGMAX] .= 360.0
    
    # 添加线路零序阻抗
    add_line_sc_impedance_zero(case, jpc_3ph, branch)
    
    # 添加变压器零序阻抗
    branch = add_trafo_sc_impedance_zero(case, jpc_3ph, branch)
    
    # 将计算结果存入JPC对象
    jpc_3ph.branchAC_0 = branch
    
    return jpc_3ph
end

function add_line_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)
    # 检查是否有线路数据
    if isempty(case.branchesAC)
        return
    end
    
    # 处理所有线路（不包括变压器）
    line_indices = findall(l -> !hasfield(typeof(l), :is_transformer) || !l.is_transformer, case.branchesAC)
    
    if isempty(line_indices)
        return
    end
    
    for i in line_indices
        line = case.branchesAC[i]
        
        # 获取起始和终止母线编号
        fb = line.from_bus
        tb = line.to_bus
        
        # 获取长度和并联数据
        length_km = line.length_km
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # 计算基准阻抗（注意除以3，与原代码保持一致）
        base_kv = jpc_3ph.busAC_0[fb, BASE_KV]
        baseR = (base_kv^2) / (3 * case.baseMVA)
        
        # 检查零序参数是否存在
        if hasfield(typeof(line), :r0_ohm_per_km) && hasfield(typeof(line), :x0_ohm_per_km) && hasfield(typeof(line), :c0_nf_per_km)
            # 计算零序阻抗
            r0_pu = length_km * line.r0_ohm_per_km / baseR / parallel
            x0_pu = length_km * line.x0_ohm_per_km / baseR / parallel
            
            # 计算零序并联电纳
            b0_pu = 2 * π * case.basef * length_km * line.c0_nf_per_km * 1e-9 * baseR * parallel
        else
            # 如果没有零序数据，使用正序数据的近似值
            r0_pu = length_km * line.r_ohm_per_km * 3.0 / baseR / parallel
            x0_pu = length_km * line.x_ohm_per_km * 3.0 / baseR / parallel
            b0_pu = 0.0
        end
        
        # 填充branch矩阵
        branch[i, F_BUS] = fb
        branch[i, T_BUS] = tb
        branch[i, BR_R] = r0_pu
        branch[i, BR_X] = x0_pu
        branch[i, BR_B] = b0_pu
        branch[i, BR_STATUS] = line.in_service ? 1.0 : 0.0
    end
    return branch  # 返回更新后的branch矩阵
end

function add_trafo_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)
    # 处理二绕组ETAP变压器的零序阻抗
    transformers = case.transformers_2w_etap
    
    if isempty(transformers)
        return branch  # 如果没有变压器，直接返回原始branch
    end
    
    # 获取线路数量，用于确定变压器在branch矩阵中的起始索引
    n_lines = count(l -> !hasfield(typeof(l), :is_transformer) || !l.is_transformer, case.branchesAC)
    
    # 检查branch矩阵是否有足够的行
    n_transformers = length(transformers)
    current_rows = size(branch, 1)
    
    if current_rows < n_lines + n_transformers
        # 如果branch矩阵行数不足，扩展矩阵
        additional_rows = n_lines + n_transformers - current_rows
        branch = vcat(branch, zeros(additional_rows, size(branch, 2)))
    end
    
    # 创建存储母线、电导和电纳的数组
    buses_all = Int[]
    gs_all = Float64[]
    bs_all = Float64[]
    
    # 处理每个变压器
    for (i, transformer) in enumerate(transformers)
        branch_idx = n_lines + i  # branch矩阵中的索引
        
        # 获取高压侧和低压侧母线编号
        hv_bus = transformer.hv_bus
        lv_bus = transformer.lv_bus
        
        # 获取向量组
        vector_group = transformer.vector_group
        
        # 跳过不支持的向量组
        if lowercase(vector_group) in ["yy", "yd", "dy", "dd"]
            continue
        end
        
        # 设置branch矩阵的基本参数
        branch[branch_idx, F_BUS] = hv_bus
        branch[branch_idx, T_BUS] = lv_bus
        
        # 默认设置变压器不参与零序网络
        branch[branch_idx, BR_STATUS] = 0.0
        
        # 获取零序阻抗参数（如果没有零序参数，使用正序参数）
        z0_percent = hasfield(typeof(transformer), :z0_percent) ? transformer.z0_percent : transformer.z_percent
        x0_r0 = hasfield(typeof(transformer), :x0_r0) ? transformer.x0_r0 : transformer.x_r
        
        # 获取额定功率(MVA)和并联数
        sn_trafo_mva = transformer.sn_mva
        parallel = transformer.parallel
        
        # 计算零序阻抗（与正序类似，但使用零序参数）
        s_ratio = sn_trafo_mva / case.baseMVA
        z0_pu = z0_percent / 100.0 / s_ratio  # 转换到系统基准
        
        r0_pu = z0_pu / sqrt(1 + x0_r0^2)
        x0_pu = r0_pu * x0_r0
        
        # 考虑并联变压器
        if parallel > 1
            r0_pu = r0_pu / parallel
            x0_pu = x0_pu / parallel
        end
        
        # 根据变压器接线类型处理零序网络
        if lowercase(vector_group) == "ynyn"
            # YNyn型变压器在零序网络中表现为普通支路
            branch[branch_idx, BR_R] = r0_pu
            branch[branch_idx, BR_X] = x0_pu
            branch[branch_idx, BR_B] = 0.0
            branch[branch_idx, BR_STATUS] = transformer.in_service ? 1.0 : 0.0
            
        elseif lowercase(vector_group) == "dyn"
            # Dyn型变压器在零序网络中，高压侧不参与，低压侧接地
            # 将等效阻抗作为分流导纳添加到低压侧母线
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, lv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "ynd"
            # YNd型变压器在零序网络中，低压侧不参与，高压侧接地
            # 将等效阻抗作为分流导纳添加到高压侧母线
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, hv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "yyn"
            # YYn型变压器在零序网络中，将等效阻抗作为分流导纳添加到低压侧母线
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, lv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "yny"
            # YNy型变压器在零序网络中，将等效阻抗作为分流导纳添加到高压侧母线
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, hv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "yzn"
            # YZn型变压器在零序网络中，需要考虑特殊的阻抗关系
            # 通常将等效阻抗乘以系数后作为分流导纳添加到低压侧母线
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                # 系数1.1547 = sqrt(3)/sqrt(2)，与变压器连接方式有关
                push!(buses_all, lv_bus)
                push!(gs_all, 1.1547 * real(y0))
                push!(bs_all, 1.1547 * imag(y0))
            end
            
        else
            # 不支持的向量组
            @warn "变压器向量组 $(vector_group) 不支持或未实现，变压器索引: $i"
        end
        
        # 设置变压器变比和相移（与正序相同）
        # 获取变压器额定电压
        vn_trafo_hv = transformer.vn_hv_kv
        vn_trafo_lv = transformer.vn_lv_kv
        
        # 获取母线额定电压
        vn_bus_hv = jpc_3ph.busAC_0[hv_bus, BASE_KV]
        vn_bus_lv = jpc_3ph.busAC_0[lv_bus, BASE_KV]
        
        # 计算变比
        ratio = (vn_trafo_hv / vn_bus_hv) / (vn_trafo_lv / vn_bus_lv)
        
        # 获取相移角度
        shift = hasfield(typeof(transformer), :shift_degree) ? transformer.shift_degree : 0.0
        
        # 设置变压器变比和相移
        branch[branch_idx, TAP] = ratio
        branch[branch_idx, SHIFT] = shift
        
        # 设置额定容量和角度限制（与正序相同）
        branch[branch_idx, RATE_A] = case.baseMVA 
        branch[branch_idx, ANGMIN] = -360.0
        branch[branch_idx, ANGMAX] = 360.0
    end
    
    # 将电导和电纳添加到母线矩阵
    # 合并相同母线的值
    if !isempty(buses_all)
        buses_unique = unique(buses_all)
        for bus in buses_unique
            indices = findall(x -> x == bus, buses_all)
            gs_sum = sum(gs_all[indices])
            bs_sum = sum(bs_all[indices])
            
            jpc_3ph.busAC_0[bus, GS] += gs_sum
            jpc_3ph.busAC_0[bus, BS] += bs_sum
        end
    end
    
    return branch  # 返回更新后的branch矩阵
end

function JPC_3ph_gens_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
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
        
        # 更新母线类型为参考节点(REF/平衡节点) - 对所有序分量都设置
        jpc_3ph.busAC_1[bus_idx, 2] = 3  # 3表示REF节点
        jpc_3ph.busAC_2[bus_idx, 2] = 3  # 3表示REF节点
        jpc_3ph.busAC_0[bus_idx, 2] = 3  # 3表示REF节点
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
        
        # 如果母线尚未设置为参考节点，则设置为PV节点 - 对所有序分量都设置
        if jpc_3ph.busAC_1[bus_idx, 2] != 3  # 3表示REF节点
            jpc_3ph.busAC_1[bus_idx, 2] = 2  # 2表示PV节点
        end
        if jpc_3ph.busAC_2[bus_idx, 2] != 3
            jpc_3ph.busAC_2[bus_idx, 2] = 2
        end
        if jpc_3ph.busAC_0[bus_idx, 2] != 3
            jpc_3ph.busAC_0[bus_idx, 2] = 2
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
        if sgen.controllable
            if jpc_3ph.busAC_1[bus_idx, 2] == 1  # 1表示PQ节点
                jpc_3ph.busAC_1[bus_idx, 2] = 2  # 2表示PV节点
            end
            if jpc_3ph.busAC_2[bus_idx, 2] == 1
                jpc_3ph.busAC_2[bus_idx, 2] = 2
            end
            if jpc_3ph.busAC_0[bus_idx, 2] == 1
                jpc_3ph.busAC_0[bus_idx, 2] = 2
            end
        end
    end
    
    # 移除未使用的行(对应未投运的发电设备)
    active_rows = findall(x -> x > 0, gen_data[:, 8])  # 第8列是GEN_STATUS
    gen_data = gen_data[active_rows, :]
    
    # 将发电机数据存储到三相JPC结构体的所有序分量
    jpc_3ph.genAC_1 = copy(gen_data)  # 正序分量
    jpc_3ph.genAC_2 = copy(gen_data)  # 负序分量
    jpc_3ph.genAC_0 = copy(gen_data)  # 零序分量
    
    # 确保至少有一个平衡节点 - 对所有序分量都设置
    if !any(jpc_3ph.busAC_1[:, 2] .== 3) && size(gen_data, 1) > 0  # 3表示REF节点
        # 如果没有平衡节点，选择第一个发电机所在母线作为平衡节点
        first_gen_bus = Int(gen_data[1, 1])
        jpc_3ph.busAC_1[first_gen_bus, 2] = 3  # 3表示REF节点
        jpc_3ph.busAC_2[first_gen_bus, 2] = 3
        jpc_3ph.busAC_0[first_gen_bus, 2] = 3
    end
end

function JPC_3ph_loads_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # 处理负荷数据，转换为三相JPC格式并更新busAC的PD和QD
    
    # 过滤出投运的负荷
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # 如果没有投运的负荷，直接返回
    if isempty(in_service_loads)
        return
    end
    
    # 创建一个空矩阵，行数为负荷数，列数为8
    num_loads = length(in_service_loads)
    load_matrix = zeros(num_loads, 8)
    
    # 创建三个字典，用于累加连接到同一母线的负荷（三个序分量）
    bus_load_sum_1 = Dict{Int, Vector{Float64}}()  # 正序
    bus_load_sum_2 = Dict{Int, Vector{Float64}}()  # 负序
    bus_load_sum_0 = Dict{Int, Vector{Float64}}()  # 零序
    
    for (i, load) in enumerate(in_service_loads)
        # 计算实际的有功和无功负荷（考虑缩放因子）
        actual_p_mw = load.p_mw * load.scaling
        actual_q_mvar = load.q_mvar * load.scaling
        
        # 根据负荷类型分配到三相序分量
        if load.type == "wye"
            # Y型连接：三相平衡负荷，主要为正序分量
            p_1 = actual_p_mw    # 正序分量（每相功率）
            q_1 = actual_q_mvar
            p_2 = 0.0                  # 负序分量（平衡负荷为零）
            q_2 = 0.0
            p_0 = 0.0                  # 零序分量（平衡负荷为零）
            q_0 = 0.0
        elseif load.type == "delta"
            # Δ型连接：线电压，零序分量为零
            p_1 = actual_p_mw    # 正序分量
            q_1 = actual_q_mvar
            p_2 = 0.0                  # 负序分量（平衡负荷为零）
            q_2 = 0.0
            p_0 = 0.0                  # 零序分量（Δ型连接无零序电流）
            q_0 = 0.0
        else
            # 默认处理为Y型连接
            p_1 = actual_p_mw
            q_1 = actual_q_mvar
            p_2 = 0.0
            q_2 = 0.0
            p_0 = 0.0
            q_0 = 0.0
        end
        
        # 填充基础负荷矩阵（使用正序分量）
        load_matrix[i, :] = [
            i,                         # 负荷编号
            load.bus,                  # 负荷连接的母线编号
            1.0,                       # 负荷状态(1=投运)
            p_1,                       # 正序有功负荷(MW)
            q_1,                       # 正序无功负荷(MVAr)
            load.const_z_percent/100,  # 恒阻抗负荷百分比
            load.const_i_percent/100,  # 恒电流负荷百分比
            load.const_p_percent/100   # 恒功率负荷百分比
        ]
        
        # 累加连接到同一母线的负荷（正序分量）
        bus_idx = load.bus
        if haskey(bus_load_sum_1, bus_idx)
            bus_load_sum_1[bus_idx][1] += p_1
            bus_load_sum_1[bus_idx][2] += q_1
        else
            bus_load_sum_1[bus_idx] = [p_1, q_1]
        end
        
        # 累加连接到同一母线的负荷（负序分量）
        if haskey(bus_load_sum_2, bus_idx)
            bus_load_sum_2[bus_idx][1] += p_2
            bus_load_sum_2[bus_idx][2] += q_2
        else
            bus_load_sum_2[bus_idx] = [p_2, q_2]
        end
        
        # 累加连接到同一母线的负荷（零序分量）
        if haskey(bus_load_sum_0, bus_idx)
            bus_load_sum_0[bus_idx][1] += p_0
            bus_load_sum_0[bus_idx][2] += q_0
        else
            bus_load_sum_0[bus_idx] = [p_0, q_0]
        end
    end
    
    # 创建三个序分量的负荷矩阵
    load_matrix_1 = copy(load_matrix)  # 正序分量负荷矩阵
    load_matrix_2 = copy(load_matrix)  # 负序分量负荷矩阵
    load_matrix_0 = copy(load_matrix)  # 零序分量负荷矩阵
    
    # 更新负序和零序负荷矩阵的功率值
    for (i, load) in enumerate(in_service_loads)
        # 负序分量矩阵更新
        load_matrix_2[i, 4] = 0.0  # 平衡负荷的负序有功功率为0
        load_matrix_2[i, 5] = 0.0  # 平衡负荷的负序无功功率为0
        
        # 零序分量矩阵更新
        load_matrix_0[i, 4] = 0.0  # 平衡负荷的零序有功功率为0
        load_matrix_0[i, 5] = 0.0  # 平衡负荷的零序无功功率为0
    end
    
    # 将负荷数据存储到三相JPC结构体的三个序分量
    jpc_3ph.loadAC_1 = load_matrix_1  # 正序分量
    jpc_3ph.loadAC_2 = load_matrix_2  # 负序分量
    jpc_3ph.loadAC_0 = load_matrix_0  # 零序分量
    
    # 更新三个序分量的busAC矩阵中的PD和QD字段
    
    # 更新正序busAC_1
    for (bus_idx, load_values) in bus_load_sum_1
        # 找到对应的母线行
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_1[:, 1])
        
        if !isnothing(bus_row)
            # 更新PD(第3列)和QD(第4列)
            jpc_3ph.busAC_1[bus_row, PD] = load_values[1]  # PD - 有功负荷(MW)
            jpc_3ph.busAC_1[bus_row, QD] = load_values[2]  # QD - 无功负荷(MVAr)
        end
    end
    
    # 更新负序busAC_2
    for (bus_idx, load_values) in bus_load_sum_2
        # 找到对应的母线行
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_2[:, 1])
        
        if !isnothing(bus_row)
            # 更新PD(第3列)和QD(第4列)
            jpc_3ph.busAC_2[bus_row, PD] = load_values[1]  # PD - 有功负荷(MW)
            jpc_3ph.busAC_2[bus_row, QD] = load_values[2]  # QD - 无功负荷(MVAr)
        end
    end
    
    # 更新零序busAC_0
    for (bus_idx, load_values) in bus_load_sum_0
        # 找到对应的母线行
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_0[:, 1])
        
        if !isnothing(bus_row)
            # 更新PD(第3列)和QD(第4列)
            jpc_3ph.busAC_0[bus_row, PD] = load_values[1]  # PD - 有功负荷(MW)
            jpc_3ph.busAC_0[bus_row, QD] = load_values[2]  # QD - 无功负荷(MVAr)
        end
    end
end

function JPC_3ph_add_grid_external_sc_impedance(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, sequence::Int)    
    # 如果没有外部电网，返回空数组和原始jpc_3ph
    if isempty(case.ext_grids)
        return jpc_3ph, Float64[], Float64[]
    end
    
    # 只在负序(sequence=2)或零序(sequence=0)时处理外部电网阻抗
    if sequence == 2 || sequence == 0
        # 收集所有外部电网的数据
        external_buses = Int[]
        Y_grid_real = Float64[]
        Y_grid_imag = Float64[]
        
        # 遍历所有外部电网
        for ext_grid in case.ext_grids
            # 获取外部电网连接的母线ID
            external_bus = ext_grid.bus
            
            # 计算外部电网的短路阻抗
            c = 1.1
            s_sc = ext_grid.s_sc_max_mva / jpc_3ph.baseMVA
            rx = ext_grid.rx_max
            z_grid = c / (s_sc / 3)
            x_grid = z_grid / sqrt(1 + rx^2)
            r_grid = x_grid * rx
            
            # 计算导纳
            Y_grid = 1 / (r_grid + 1im * x_grid)
            
            push!(external_buses, external_bus)
            push!(Y_grid_real, real(Y_grid))
            push!(Y_grid_imag, imag(Y_grid))
        end
        
        # 对相同母线的导纳值进行汇总
        buses, gs, bs = PowerFlow.sum_by_group(external_buses, Y_grid_real, Y_grid_imag)
        
        # 根据序分量更新对应的母线数据
        if sequence == 2
            # 更新负序busAC_2中的母线数据
            for (i, bus_id) in enumerate(buses)
                # 找到对应母线在jpc_3ph.busAC_2中的索引
                bus_idx = findfirst(x -> x == bus_id, jpc_3ph.busAC_2[:, 1])
                
                if !isnothing(bus_idx)
                    # 更新母线的电导和电纳值
                    jpc_3ph.busAC_2[bus_idx, GS] = gs[i] * jpc_3ph.baseMVA
                    jpc_3ph.busAC_2[bus_idx, BS] = bs[i] * jpc_3ph.baseMVA
                end
            end
        elseif sequence == 0
            # 更新零序busAC_0中的母线数据
            for (i, bus_id) in enumerate(buses)
                # 找到对应母线在jpc_3ph.busAC_0中的索引
                bus_idx = findfirst(x -> x == bus_id, jpc_3ph.busAC_0[:, 1])
                
                if !isnothing(bus_idx)
                    # 更新母线的电导和电纳值
                    jpc_3ph.busAC_0[bus_idx, GS] = gs[i] * jpc_3ph.baseMVA
                    jpc_3ph.busAC_0[bus_idx, BS] = bs[i] * jpc_3ph.baseMVA
                end
            end
        end
        
        # 返回更新后的jpc_3ph和计算出的电导、电纳值
        gs_eg = gs .* jpc_3ph.baseMVA
        bs_eg = bs .* jpc_3ph.baseMVA
        return jpc_3ph, gs_eg, bs_eg
    else
        # 正序(sequence=1)时不处理外部电网阻抗，直接返回
        return jpc_3ph, Float64[], Float64[]
    end
end
