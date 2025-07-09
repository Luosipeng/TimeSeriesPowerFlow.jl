function makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm_gpu, load_gpu; dc=false, Sg=nothing, return_derivative=false)
    # Define named indices into bus, gen matrices
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
        PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
        RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
         COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN,GEN_AREA) = PowerFlow.idx_gen();
    (LOAD_I,LOAD_CND,LOAD_STATUS,LOAD_PD,LOAD_QD,LOADZ_PERCENT,LOADI_PERCENT,
    LOADP_PERCENT)=PowerFlow.idx_ld()

    nb = size(bus_gpu, 1)
    pw_1=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_2=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_3=PowerFlow.CUDA.zeros(size(bus_gpu,1),1)
    pw_1[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADP_PERCENT]
    pw_2[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADI_PERCENT]
    pw_3[Int64.(load_gpu[:,LOAD_CND])]=load_gpu[:,LOADZ_PERCENT]
    # Get load parameters
    Sd = makeSdzip_gpu(baseMVA, bus_gpu,pw_1,pw_2,pw_3)

    if return_derivative
        if isempty(Vm_gpu)
            dSbus_dVm = PowerFlow.CUDA.spzeros(nb, nb)
        else
            diag_elements = Sd.i + 2 .* Vm_gpu .* Sd.z
            dSbus_dVm = -PowerFlow.Diagonal(diag_elements)
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
            # 步骤1: 创建发电机复功率向量
            Sg = gen_gpu[on, PG] .+ 1im * gen_gpu[on, QG]

            # 步骤2: 创建结果向量 (母线注入功率)
            Sbusg = CUDA.zeros(ComplexF64, size(Cg, 1))

            # 步骤3: 使用 CUSPARSE.mv! 函数执行矩阵-向量乘法
            # 添加额外的字符参数 'O' 表示操作类型
            CUDA.CUSPARSE.mv!('N', one(ComplexF64), Cg, Sg, zero(ComplexF64), Sbusg, 'O')

            # 步骤4: 除以基准功率 baseMVA 进行标幺化
            Sbusg = Sbusg ./ baseMVA
        end

        if dc
            Vm = PowerFlow.CUDA.ones(nb,1)
        end
        # Compute per-bus loads in p.u.
        Sbusd = Sd.p .+ Sd.i .* Vm_gpu .+ Sd.z .* Vm_gpu.^2

        # Form net complex bus power injection vector
        # (power injected by generators + power injected by loads)
        Sbus = Sbusg - Sbusd
        return Sbus
    end
end
