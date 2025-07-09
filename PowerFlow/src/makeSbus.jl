function makeSbus(baseMVA, bus, gen, Vm, load::Matrix{Float64}, pvarray; dc=false, Sg=nothing, return_derivative=false)


    nb = size(bus, 1)
    pw_1=zeros(size(bus,1),1)
    pw_2=zeros(size(bus,1),1)
    pw_3=zeros(size(bus,1),1)
    pw_1[Int64.(load[:,LOAD_CND])]=load[:,LOADP_PERCENT]
    pw_2[Int64.(load[:,LOAD_CND])]=load[:,LOADI_PERCENT]
    pw_3[Int64.(load[:,LOAD_CND])]=load[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = makeSdzip(baseMVA, bus,pw_1,pw_2,pw_3)

     # 初始化光伏功率和导数
    Spv = zeros(Complex{Float64}, nb)
    dSpv_dVm = spzeros(nb, nb)
    
    # 如果有光伏阵列且不是直流潮流计算，计算其功率和导数
    if !isnothing(pvarray)
        Spv, dSpv_dVm = calculate_pv_power(pvarray, bus, Vm, baseMVA)
    end

    if return_derivative
        if isempty(Vm)
            dSbus_dVm = spzeros(nb, nb)
        else
            # 负荷导数 + 光伏导数
            dSbus_dVm = -(spdiagm(0 => Sd.i + 2 .* Vm .* Sd.z)) + dSpv_dVm
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = sparse(gbus, 1:ngon, 1, nb, ngon)  # connection matrix
        # element i, j is 1 if gen on(j) at bus i is ON
        if Sg !== nothing
            Sbusg = Cg * Sg[on]
        else
            Sbusg = Cg * (gen[on, PG] .+ 1im * gen[on, QG]) / baseMVA
        end

        if dc
            Vm = ones(nb,1)
        end
        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm .+ Sd.z .* Vm.^2

       # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads + PV power)
        Sbus = Sbusg - Sbusd + Spv
        return Sbus
    end
end

function makeSbus(baseMVA, bus, gen, Vm, Sg=nothing, return_derivative=false)
    # Define named indices into bus, gen matrices
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = idx_bus();
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
     PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
     RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
      COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN,GEN_AREA) = idx_gen();

    nb = size(bus, 1)

    # Get load parameters
    Sd = makeSdzip(baseMVA, bus)

    if return_derivative
        if isempty(Vm)
            dSbus_dVm = spzeros(nb, nb)
        else
            dSbus_dVm = -(spdiagm(0 => Sd.i + 2 .* Vm .* Sd.z))
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = sparse(gbus, 1:ngon, 1, nb, ngon)  # connection matrix
        # element i, j is 1 if gen on(j) at bus i is ON
        if !isnothing(Sg)
            Sbusg = Cg * Sg[on]
        else
            Sbusg = Cg * (gen[on, PG] .+ 1im * gen[on, QG]) / baseMVA
        end

        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm .+ Sd.z .* Vm.^2

        # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads)
        Sbus = Sbusg - Sbusd
        return Sbus
    end
end

function calculate_pv_power(pvarray, bus, Vm, baseMVA)
    
    nb = size(Vm, 1)
    Spv = zeros(Complex{Float64}, nb)
    dSpv_dVm = spzeros(nb, nb)
    
    # 筛选在服务中的光伏装置
    active_pv = pvarray[pvarray[:, PV_IN_SERVICE] .> 0, :]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # 创建母线号到索引的映射
    bus_num_to_idx = Dict{Int, Int}()
    for i in eachindex(bus[:,1])
        bus_num_to_idx[Int(bus[i, BUS_I])] = i
    end
    
    # 获取光伏连接的母线索引
    bus_nums = Int.(active_pv[:, PV_BUS])
    valid_pv = falses(size(active_pv, 1))
    bus_indices = zeros(Int, size(active_pv, 1))
    
    for i in eachindex(bus_nums)
        if haskey(bus_num_to_idx, bus_nums[i])
            valid_pv[i] = true
            bus_indices[i] = bus_num_to_idx[bus_nums[i]]
        end
    end
    
    # 只保留有效的光伏装置
    active_pv = active_pv[valid_pv, :]
    bus_indices = bus_indices[valid_pv]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # 获取母线基准电压 (kV)
    base_kvs = bus[bus_indices, BASE_KV]
    
    # 获取当前母线电压 (标幺值)
    V_buses = Vm[bus_indices]
    
    # 修正幂函数模型参数
    a = 10.000000
    b = 0.547596
    c = 0.023812
    
    # 提取光伏参数
    Vocs = active_pv[:, PV_VOC]      # 开路电压 (V)
    Iscs = active_pv[:, PV_ISC]      # 短路电流 (A)
    areas = active_pv[:, PV_AREA]    # 面积或其他缩放因子
    Vmpps = active_pv[:, PV_VMPP]    # 最大功率点电压 (V)，需要确保这个列存在
    
    # 将标幺值电压转换为实际电压 (V)
    voltage_ratios = 1
    V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # 转换为伏特
    
    # 初始化电流和导数数组
    I_arrays = zeros(size(active_pv, 1))
    dI_dVs = zeros(size(active_pv, 1))
    
    # 创建掩码进行向量化操作
    neg_v_mask = V_arrays .< 0
    valid_v_mask = (V_arrays .>= 0) .& (V_arrays .<= Vocs)
    over_v_mask = V_arrays .> Vocs
    
    # 计算比例
    v_ratios = zeros(size(V_arrays))
    vmpp_ratios = zeros(size(V_arrays))
    v_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vocs[valid_v_mask]
    vmpp_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vmpps[valid_v_mask] .- 1
    
    # 使用修正幂函数模型计算电流
    for i in 1:length(V_arrays)
        if valid_v_mask[i]
            # 基础幂函数部分
            base_term = (1 - v_ratios[i]^a)^b
            
            # 修正项
            correction_term = (1 - c * vmpp_ratios[i]^2)
            
            # 计算电流
            I_arrays[i] = Iscs[i] * base_term * correction_term
            
            # 确保电流不为负
            I_arrays[i] = max(0, I_arrays[i])
        end
    end
    
    # 计算导数（对于有效电压范围内的值）
    for i in eachindex(V_arrays)
        if valid_v_mask[i] && V_arrays[i] > 0
            # 基础幂函数部分的导数
            d_base_term_dv = -a * b * (1 - v_ratios[i]^a)^(b-1) * v_ratios[i]^(a-1) / Vocs[i]
            
            # 修正项的导数
            d_correction_term_dv = -c * 2 * vmpp_ratios[i] / Vmpps[i]
            
            # 使用乘积法则计算总导数
            base_term = (1 - v_ratios[i]^a)^b
            correction_term = (1 - c * vmpp_ratios[i]^2)
            
            dI_dVs[i] = Iscs[i] * (
                d_base_term_dv * correction_term + 
                base_term * d_correction_term_dv
            )
            
            # 如果电流为0，导数也应该为0
            if I_arrays[i] <= 0
                dI_dVs[i] = 0
            end
        end
    end
    
    # 计算功率 (W)
    P_arrays = V_arrays .* I_arrays
    
    # 换算到标幺值 (考虑baseMVA)
    P_pus = P_arrays ./ (baseMVA * 1e6)  # 从W转换为标幺值
    
    # 计算功率对电压的导数
    dP_dVs = I_arrays .+ V_arrays .* dI_dVs
    dP_dV_pus = (dP_dVs .* base_kvs .* 1000 .* voltage_ratios) ./ (baseMVA * 1e6)
    
    # 更新功率注入和导数矩阵
    for i in eachindex(bus_indices)
        Spv[bus_indices[i]] += P_pus[i] + 0im  # 只有有功功率
        dSpv_dVm[bus_indices[i], bus_indices[i]] += dP_dV_pus[i]
    end
    
    return Spv, dSpv_dVm
end


# function calculate_pv_power(pvarray, bus, Vm, baseMVA)
    
#     nb = size(Vm, 1)
#     Spv = zeros(Complex{Float64}, nb)
#     dSpv_dVm = spzeros(nb, nb)
    
#     # 筛选在服务中的光伏装置
#     active_pv = pvarray[pvarray[:, PV_IN_SERVICE] .> 0, :]
    
#     if isempty(active_pv)
#         return Spv, dSpv_dVm
#     end
    
#     # 创建母线号到索引的映射
#     bus_num_to_idx = Dict{Int, Int}()
#     for i in 1:size(bus, 1)
#         bus_num_to_idx[Int(bus[i, BUS_I])] = i
#     end
    
#     # 获取光伏连接的母线索引
#     bus_nums = Int.(active_pv[:, PV_BUS])
#     valid_pv = falses(size(active_pv, 1))
#     bus_indices = zeros(Int, size(active_pv, 1))
    
#     for i in 1:length(bus_nums)
#         if haskey(bus_num_to_idx, bus_nums[i])
#             valid_pv[i] = true
#             bus_indices[i] = bus_num_to_idx[bus_nums[i]]
#         end
#     end
    
#     # 只保留有效的光伏装置
#     active_pv = active_pv[valid_pv, :]
#     bus_indices = bus_indices[valid_pv]
    
#     if isempty(active_pv)
#         return Spv, dSpv_dVm
#     end
    
#     # 获取母线基准电压 (kV)
#     base_kvs = bus[bus_indices, BASE_KV]
    
#     # 获取当前母线电压 (标幺值)
#     V_buses = Vm[bus_indices]
    
#     # 经验公式参数
#     a = 14.5
    
#     # 提取光伏参数
#     Vocs = active_pv[:, PV_VOC]      # 开路电压 (V)
#     Iscs = active_pv[:, PV_ISC]      # 短路电流 (A)
#     areas = active_pv[:, PV_AREA]    # 面积或其他缩放因子
    
#     # 将标幺值电压转换为实际电压 (V)
#     voltage_ratios = 1
#     V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # 转换为伏特
    
#     # 初始化电流和导数数组
#     I_arrays = zeros(size(active_pv, 1))
#     dI_dVs = zeros(size(active_pv, 1))
    
#     # 创建掩码进行向量化操作
#     neg_v_mask = V_arrays .< 0
#     valid_v_mask = (V_arrays .>= 0) .& (V_arrays .<= Vocs)
#     over_v_mask = V_arrays .> Vocs
    
#     # 计算比例
#     ratios = zeros(size(V_arrays))
#     ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vocs[valid_v_mask]
    
#     # 使用经验公式计算电流（向量化操作）
#     I_arrays[valid_v_mask] = Iscs[valid_v_mask] .* areas[valid_v_mask] .* (1 .- ratios[valid_v_mask].^a)
    
#     # 计算导数（向量化操作）
#     # 避免0处的导数问题
#     nonzero_v = (V_arrays .> 0) .& valid_v_mask
#     dI_dVs[nonzero_v] = -Iscs[nonzero_v] .* areas[nonzero_v] .* a .* ratios[nonzero_v].^(a-1) ./ Vocs[nonzero_v]
    
#     # 计算功率 (W)
#     P_arrays = V_arrays .* I_arrays
    
#     # 换算到标幺值 (考虑baseMVA)
#     P_pus = P_arrays ./ (baseMVA * 1e6)  # 从W转换为标幺值
    
#     # 计算功率对电压的导数
#     dP_dVs = I_arrays .+ V_arrays .* dI_dVs
#     dP_dV_pus = (dP_dVs .* base_kvs .* 1000 .* voltage_ratios) ./ (baseMVA * 1e6)
    
#     # 更新功率注入和导数矩阵
#     for i in 1:length(bus_indices)
#         Spv[bus_indices[i]] += P_pus[i] + 0im  # 只有有功功率
#         dSpv_dVm[bus_indices[i], bus_indices[i]] += dP_dV_pus[i]
#     end
    
#     return Spv, dSpv_dVm
# end