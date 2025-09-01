"""
    makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu; dc=false, Sg=nothing, return_derivative=false)

Build the vector of complex bus power injections using GPU acceleration.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `gen_gpu`: Generator data matrix on GPU with columns representing generator parameters
- `gen`: Generator data matrix on CPU with columns representing generator parameters
- `Vm_gpu`: Vector of bus voltage magnitudes on GPU
- `load_gpu`: Load data matrix on GPU with columns representing load parameters

# Keyword Arguments
- `dc`: Boolean indicating whether to use DC power flow assumptions (default: false)
- `Sg`: Optional pre-computed generator complex power injections (default: nothing)
- `return_derivative`: Boolean indicating whether to return derivative of Sbus with respect to Vm (default: false)

# Returns
- If `return_derivative=false`: Vector of complex bus power injections (Sbus)
- If `return_derivative=true`: Sparse matrix of partial derivatives of power injections with respect to voltage magnitude (dSbus_dVm)

# Description
This function computes the vector of complex bus power injections (Sbus) for power flow analysis using GPU acceleration.
It accounts for ZIP load models (constant power, constant current, and constant impedance components) and generator injections.

When `return_derivative=true`, it returns the partial derivatives of the power injections with respect to voltage magnitude,
which is useful for power flow Jacobian calculations.

# Notes
- All power values are converted to per-unit on system MVA base
- The function handles ZIP load models with percentages specified in load_gpu
- Generator status is considered when computing injections
- When dc=true, voltage magnitudes are set to 1.0 p.u.

# Constants Used (assumed to be defined elsewhere)
- LOAD_CND: Column index for load bus number in load_gpu matrix
- LOADP_PERCENT: Column index for constant power percentage in load_gpu matrix
- LOADI_PERCENT: Column index for constant current percentage in load_gpu matrix
- LOADZ_PERCENT: Column index for constant impedance percentage in load_gpu matrix
- GEN_STATUS: Column index for generator status in gen matrix
- GEN_BUS: Column index for generator bus number in gen matrix
- PG: Column index for real power output in gen_gpu matrix
- QG: Column index for reactive power output in gen_gpu matrix
"""
function makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu, pvarray_gpu; dc=false, Sg=nothing, return_derivative=false)

    nb = size(bus_gpu, 1)
    pw_1=CUDA.zeros(size(bus_gpu,1),1)
    pw_2=CUDA.zeros(size(bus_gpu,1),1)
    pw_3=CUDA.zeros(size(bus_gpu,1),1)
    pw_1[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADP_PERCENT]
    pw_2[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADI_PERCENT]
    pw_3[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = PowerFlow.makeSdzip_gpu(baseMVA, bus_gpu,pw_1,pw_2,pw_3)
    dSpv_dVm_gpu = CUDA.CUSPARSE.spzeros(Complex{Float64}, nb)

    # If there are PV arrays and not using DC power flow, calculate their power and derivatives
    if !isnothing(pvarray_gpu)
        Spv, dSpv_dVm = PowerFlow.calculate_pv_power_gpu(pvarray_gpu, bus_gpu, Vm_gpu, baseMVA)
    end

    if return_derivative
        if isempty(Vm_gpu)
            dSbus_dVm = PowerFlow.CUDA.spzeros(nb, nb)
        else
            dSpv_dVm_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(dSpv_dVm)
            diag_elements = Sd.i + 2 .* Vm_gpu .* Sd.z
            dSbus_dVm = -PowerFlow.Diagonal(diag_elements) + dSpv_dVm_gpu
        end
        return dSbus_dVm
    else
        # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = CUSPARSE.CuSparseMatrixCSR(sparse(Int64.(gbus), collect(1:ngon), ones(ngon), nb, ngon))

        # element i, j is 1 if gen on(j) at bus i is ON
        if Sg !== nothing
            Sbusg = Cg * Sg[on]
        else
            # Step 1: Create generator complex power vector
            Sg = gen_gpu[on, PG] .+ 1im * gen_gpu[on, QG]

            # Step 2: Create result vector (bus injection power)
            Sbusg = CUDA.zeros(ComplexF64, size(Cg, 1))

            # Step 3: Use CUSPARSE.mv! function to perform matrix-vector multiplication
            # Add extra character parameter 'O' to indicate operation type
            CUDA.CUSPARSE.mv!('N', one(ComplexF64), Cg, Sg, zero(ComplexF64), Sbusg, 'O')

            # Step 4: Divide by base power baseMVA for per-unit normalization
            Sbusg = Sbusg ./ baseMVA
        end

        if dc
            Vm = PowerFlow.CUDA.ones(nb,1)
        end
        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm_gpu .+ Sd.z .* Vm_gpu.^2

        # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads)
        Sbus = Sbusg - Sbusd + Spv
        return Sbus
    end
end

function calculate_pv_power_gpu(pvarray_gpu, bus_gpu, Vm_gpu, baseMVA)
    nb = size(Vm_gpu, 1)
    Spv = CUDA.zeros(Complex{Float64}, nb)
    dSpv_dVm = CUDA.CUSPARSE.spzeros(nb, nb)
    
    # 过滤在服务中的PV阵列
    active_pv = pvarray_gpu[pvarray_gpu[:, PV_IN_SERVICE] .> 0, :]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # 创建从总线号到索引的映射
    valid_pv, bus_indices = process_pv_connections(bus_gpu, active_pv, BUS_I, PV_BUS)
    
    # 只保留有效的PV阵列
    active_pv = active_pv[valid_pv, :]
    bus_indices = bus_indices[valid_pv]
    
    if isempty(active_pv)
        return Spv, dSpv_dVm
    end
    
    # 将必要的数据传回CPU进行处理
    bus_cpu = Array(bus_gpu)
    Vm_cpu = Array(Vm_gpu)
    active_pv_cpu = Array(active_pv)
    
    # 获取总线的基准电压 (kV)
    base_kvs = bus_cpu[bus_indices, BASE_KV]
    
    # 获取当前总线电压 (标幺值)
    V_buses = Vm_cpu[bus_indices]
    
    # 修改的功率函数模型参数
    a = 10.000000
    b = 0.547596
    c = 0.023812
    
    # 提取PV参数
    Vocs = active_pv_cpu[:, PV_VOC]      # 开路电压 (V)
    Iscs = active_pv_cpu[:, PV_ISC]      # 短路电流 (A)
    areas = active_pv_cpu[:, PV_AREA]    # 面积或其他缩放因子
    Vmpps = active_pv_cpu[:, PV_VMPP]    # 最大功率点电压 (V)
    
    # 将标幺值电压转换为实际电压 (V)
    voltage_ratios = 1
    V_arrays = V_buses .* base_kvs .* 1000 .* voltage_ratios  # 转换为伏特
    
    # 初始化电流和导数数组
    I_arrays = zeros(size(active_pv_cpu, 1))
    dI_dVs = zeros(size(active_pv_cpu, 1))
    
    # 创建用于向量化操作的掩码
    neg_v_mask = real.(V_arrays) .< 0
    valid_v_mask = (real.(V_arrays) .>= 0) .& (real.(V_arrays) .<= Vocs)
    over_v_mask = real.(V_arrays) .> Vocs
    
    # 计算比率
    v_ratios = zeros(size(V_arrays))
    vmpp_ratios = zeros(size(V_arrays))
    v_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vocs[valid_v_mask]
    vmpp_ratios[valid_v_mask] = V_arrays[valid_v_mask] ./ Vmpps[valid_v_mask] .- 1
    
    # 使用修改的功率函数模型计算电流
    for i in 1:length(V_arrays)
        if valid_v_mask[i]
            # 基础功率函数项
            base_term = (1 - v_ratios[i]^a)^b
            
            # 修正项
            correction_term = (1 - c * vmpp_ratios[i]^2)
            
            # 计算电流
            I_arrays[i] = Iscs[i] * base_term * correction_term
            
            # 确保电流非负
            I_arrays[i] = max(0, I_arrays[i])
        end
    end
    
    # 计算导数（对于有效电压范围）
    for i in eachindex(V_arrays)
        if valid_v_mask[i] && real(V_arrays[i]) > 0
            # 基础功率函数项的导数
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
            
            # 如果电流为0，导数也应为0
            if I_arrays[i] <= 0
                dI_dVs[i] = 0
            end
        end
    end
    
    # 计算功率 (W)
    P_arrays = V_arrays .* I_arrays
    
    # 转换为标幺值（考虑baseMVA）
    P_pus = P_arrays ./ (baseMVA * 1e6)  # 从W转换为标幺值
    
    # 计算功率对电压的导数
    dP_dVs = I_arrays .+ V_arrays .* dI_dVs
    dP_dV_pus = (dP_dVs .* base_kvs .* 1000 .* voltage_ratios) ./ (baseMVA * 1e6)
    
    # 更新功率注入和导数矩阵
    Spv_cpu = Array(Spv)
    for i in eachindex(bus_indices)
        Spv_cpu[bus_indices[i]] += P_pus[i] + 0im  # 仅实功率
    end
    
    # 将结果传回GPU
    copyto!(Spv, Spv_cpu)
    
    # 处理稀疏导数矩阵
    I_indices = bus_indices
    J_indices = bus_indices
    V_values = dP_dV_pus
    
    # 创建新的稀疏矩阵
    dSpv_dVm = CUDA.CUSPARSE.sparse(I_indices, J_indices, V_values, nb, nb)
    
    return Spv, dSpv_dVm
end
# 创建映射矩阵
function create_bus_mapping(bus_gpu, BUS_I)
    # 将bus数据传回CPU处理映射关系
    bus_ids = Array(Int.(bus_gpu[:, BUS_I]))
    n = length(bus_ids)
    
    # 创建映射字典
    bus_to_idx = Dict{Int, Int}()
    for i in 1:n
        bus_to_idx[bus_ids[i]] = i
    end
    
    return bus_to_idx
end

# 查找总线索引 - 完全在CPU上处理
function find_bus_indices(bus_to_idx, active_pv, PV_BUS)
    # 将PV总线号传回CPU
    bus_nums = Array(Int.(active_pv[:, PV_BUS]))
    n = length(bus_nums)
    
    # 创建结果数组
    valid_pv = falses(n)
    bus_indices = zeros(Int, n)
    
    # 处理每个总线号
    for i in 1:n
        bus_id = bus_nums[i]
        if haskey(bus_to_idx, bus_id)
            valid_pv[i] = true
            bus_indices[i] = bus_to_idx[bus_id]
        end
    end
    
    return valid_pv, bus_indices
end

# 处理PV连接
function process_pv_connections(bus_gpu, active_pv, BUS_I, PV_BUS)
    # 创建映射字典
    bus_to_idx = create_bus_mapping(bus_gpu, BUS_I)
    
    # 查找总线索引
    valid_pv, bus_indices = find_bus_indices(bus_to_idx, active_pv, PV_BUS)
    
    return valid_pv, bus_indices
end