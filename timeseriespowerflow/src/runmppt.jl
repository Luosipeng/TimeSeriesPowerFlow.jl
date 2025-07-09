
using  JuMP
using  Ipopt
using  LinearAlgebra

function runmppt(jpc)
    # Extract the DC parameters from the jpc structure
    branchDC = jpc.branchDC 
    busDC = jpc.busDC
    converter = jpc.converter
    loadDC = jpc.loadDC
    genDC = jpc.genDC
    pvarray = jpc.pv

    #ext2int
    genDC  = genDC[genDC[:, GEN_STATUS] .!= 0, :]
    busDC  = busDC[busDC[:, BUS_TYPE] .!= 0, :]
    loadDC = loadDC[loadDC[:, LOAD_STATUS] .!= 0, :]
    branchDC = branchDC[branchDC[:, BR_STATUS] .!= 0, :]
    converter = converter[converter[:, 3] .!= 0, :]
    if size(pvarray, 1) > 0
        pvarray = pvarray[pvarray[:, PV_IN_SERVICE] .!= 0, :]
    end

    # remove zero load and generation
    genDC = genDC[genDC[:, 8] .!= 0, :]
    branchDC = branchDC[branchDC[:, 11] .!= 0, :]
    # create map of external bus numbers to bus indices
    i2e = Int.(busDC[:, BUS_I])  # 确保i2e是整数类型
    e2i = sparsevec(zeros(Int, Int(maximum(i2e))))
    e2i[Int.(i2e)] = 1:size(busDC, 1)
    # renumber buses consecutively
    busDC[:, BUS_I] = e2i[busDC[:, BUS_I]]
    genDC[:, GEN_BUS] = e2i[genDC[:, GEN_BUS]]
    branchDC[:, F_BUS] = e2i[branchDC[:, F_BUS]]
    branchDC[:, T_BUS] = e2i[branchDC[:, T_BUS]]
    loadDC[:, LOAD_CND] = e2i[loadDC[:, LOAD_CND]]
    converter[:,2] = e2i[converter[:,2]]
    if size(pvarray, 1) > 0
        pvarray[:, PV_BUS] = e2i[pvarray[:, PV_BUS]]
    end

    # Initialize
    n_branch = size(branchDC, 1)
    n_bus = size(busDC, 1)
    n_converter = size(converter, 1)
    n_load = size(loadDC, 1)
    n_gen = size(genDC, 1)
    Cd = zeros(n_bus, n_load)

    slack_bus = findall(busDC[:, 2] .== 2)

    for i in 1:n_load
        Cd[Int(loadDC[i, 2]), i] = 1
    end

    Cg = zeros(n_bus, n_gen)
    for i in 1:n_gen
        Cg[Int(genDC[i, 1]), i] = 1
    end

    Crec = zeros(n_bus, n_converter)
    for i in 1:n_converter
        Crec[Int(converter[i, 2]), i] = 1
    end

    Pd = loadDC[:, 4]/jpc.baseMVA
    PG = genDC[:, 2]/jpc.baseMVA

    # Calculate the addmitance matrix
    Y , Yf, Yt = PowerFlow.makeYbus(jpc.baseMVA, busDC, branchDC)

    # 创建优化模型
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 5)

    # 定义变量
    @variable(model, pdc_mw[1:n_converter])  # 换流器的功率
    @variable(model, vm[1:n_bus] >= 0)  # 节点电压

    @variable(model, slack[1:n_bus])
    @constraint(model, vm .* (Y * vm) + Cd * Pd - Cg * PG + Crec * pdc_mw .== slack)
    @constraint(model, -1e-3 .<= slack .<= 1e-3)


    # 电压约束
    @constraint(model, vm .>= 0.85)
    @constraint(model, vm .<= 1.05)
    @constraint(model, vm[slack_bus] .== 1.0)  # 平衡节点电压约束

    @objective(model, Max, sum(pdc_mw))  # 目标函数：最大化功率传输

    # 设置初始值
    for i in 1:n_bus
    set_start_value(vm[i], 1.0)
    end
    
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-6)
    set_optimizer_attribute(model, "acceptable_tol", 1e-4)
    set_optimizer_attribute(model, "acceptable_iter", 15)
    set_optimizer_attribute(model, "hessian_approximation", "limited-memory")

    # 求解模型
    optimize!(model)
    # 检查求解状态
    if has_values(model) && (termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL])
        # 获取结果
        pdc_mw = value.(pdc_mw)*jpc.baseMVA
        println("换流器功率传输结果：")
        println(value.(vm))
        
        return pdc_mw
    else
        error("MPPT optimization failed")
    end
end