using HiGHS
using JuMP 
using PowerFlow
using SparseArrays
using Random
using LinearAlgebra
using Printf
using Base.Threads

struct NetworkVariableIndices
    # Flow variables
    Fij_start::Int
    Fg_start::Int 
    Fm_start::Int
    
    # Binary variables
    α_start::Int
    β_start::Int
    
    # Power variables
    pls_start::Int
    qls_start::Int
    pij_start::Int
    qij_start::Int
    pg_start::Int
    qg_start::Int
    pm_start::Int
    qm_start::Int
    
    # Voltage variables
    v_start::Int
    
    # Total number of variables
    n_vars::Int
    
    function NetworkVariableIndices(nl::Int, ng::Int, nmg::Int, nd::Int, nb::Int)
        # Define sequential indexing
        Fij_start = 1
        Fg_start = Fij_start + nl
        Fm_start = Fg_start + ng
        α_start = Fm_start + nmg
        β_start = α_start + nl
        pls_start = β_start + nl
        qls_start = pls_start + nd
        pij_start = qls_start + nd
        qij_start = pij_start + nl
        pg_start = qij_start + nl
        qg_start = pg_start + ng
        pm_start = qg_start + ng
        qm_start = pm_start + nmg
        v_start = qm_start + nmg
        n_vars = v_start + nb - 1
        
        new(
            Fij_start, Fg_start, Fm_start,
            α_start, β_start,
            pls_start, qls_start,
            pij_start, qij_start,
            pg_start, qg_start,
            pm_start, qm_start,
            v_start,
            n_vars
        )
    end
end

# Helper functions to get variable ranges
function get_range(start::Int, length::Int)
    return start:start+length-1
end

function get_variable_values(x::Vector, idx::NetworkVariableIndices, var_name::Symbol, dim::Int)
    start_idx = getfield(idx, Symbol(var_name, "_start"))
    return x[get_range(start_idx, dim)]
end

function network_reconfiguration(mpc::Dict{String, Any}, mg_buses::Vector{Int}, tie_lines::Vector{Int}, ζ::Vector{Float64}, Betamg::Vector{Int}; verbose::Bool=false)
    # Get PowerFlow indices
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = idx_bus();
    
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, 
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, 
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX) = PowerFlow.idx_brch()
    
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
     PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
     RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
      COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN,GEN_AREA) = PowerFlow.idx_gen();

    # Extract network data
    bus, gen, branch = mpc["bus"], mpc["gen"], mpc["branch"]
    baseMVA = mpc["baseMVA"]
    
    # System dimensions
    nb = size(bus,1)  # Number of buses
    nl = size(branch,1)  # Number of branches
    ng = size(gen,1)  # Number of generators
    load_buses = findall(bus[:,PD].>0)
    nd = length(load_buses)  # Number of demand buses
    nmg = length(mg_buses)  # Number of microgrids
    
    # Connection matrices
    Cg = sparse(1:ng, gen[:,GEN_BUS], ones(ng), ng, nb)
    Cd = sparse(1:nd, load_buses, ones(nd), nd, nb)
    Cmg = sparse(1:nmg, mg_buses, ones(nmg), nmg, nb)
    f = branch[:,F_BUS]
    t = branch[:,T_BUS]
    i = [1:nl; 1:nl]
    Cft = sparse(i, [f; t], [ones(nl); -ones(nl)], nl, nb)

    # Parameters
    Pd = bus[load_buses,PD]/baseMVA
    Qd = bus[load_buses,QD]/baseMVA
    Fd = ones(nd)
    
    # Generator parameters
    Pgmax = sum(Pd)*ones(ng)
    Pgmin = zeros(ng)
    Qgmax = sum(Qd)*ones(ng)
    Qgmin = zeros(ng)
    Fgmax = ones(ng)*nb
    Fgmin = zeros(ng)
    
    # Microgrid parameters
    Pmgmax = sum(Pd)*ones(nmg)
    Pmgmin = zeros(nmg)
    Qmgmax = sum(Qd)*ones(nmg)
    Qmgmin = zeros(nmg)
    Fmgmax = ones(nmg)*nb
    Fmgmin = zeros(nmg)

    # Branch parameters
    R = branch[:, BR_R]
    X = branch[:, BR_X]
    Smax = branch[:,RATE_A]
    for i = 1:nl
        if Smax[i] == 0
            Smax[i] = sum(Pd)*baseMVA
        end
    end
    
    # Bus parameters
    VMAX = 1.1^2*ones(nb)
    VMIN = 0.9^2*ones(nb)
    bigM = 100000

    # Create optimization model
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # Decision variables
    @variable(model, -nb <= Fij[i = 1:nl] <= nb)
    @variable(model, 0 <= Fg[i=1:ng] <= nb)
    @variable(model, 0 <= Fm[i=1:nmg] <= nb)
    @variable(model, α[i=1:nl], Bin)
    @variable(model, β[i=1:nl], Bin)
    @variable(model, 0 <= pls[i=1:nd] .<= Pd[i])
    @variable(model, 0 <= qls[i=1:nd] .<= Qd[i])
    @variable(model, -Smax[i] <= pij[i=1:nl] <= Smax[i])
    @variable(model, -Smax[i] <= qij[i=1:nl] <= Smax[i])
    @variable(model, 0 <= pg[i=1:ng] <= Pgmax[i])
    @variable(model, -Qgmax[i] <= qg[i=1:ng] <= Qgmax[i])
    @variable(model, 0 <= pm[i=1:nmg] <= Pmgmax[i])
    @variable(model, -Qmgmax[i] <= qm[i=1:nmg] <= Qmgmax[i])
    @variable(model, VMIN[i] <= v[i=1:nb] <= VMAX[i])

    # Objective function
    @objective(model, Min, sum(pls) + sum(β[tie_lines]));

    # Constraints
    # (1) fictious flow conservation
    @constraint(model, Cft' * Fij .== Cg'*Fg - Cd'*Fd + Cmg'*Fm);
    # (2) Fij <= Smax.*α
    @constraint(model, Fij.<= nb.*α);
    # (3) Fij >= -Smax.*α
    @constraint(model, Fij.>= -nb.*α);
    # (4) Fmg .<= Betamg.*Fmgmax
    @constraint(model, Fm .<= Fmgmax);
    # (5) active power balance
    @constraint(model, Cft'*pij .== Cg'*pg - Cd'*Pd + Cmg'*pm + Cd'*pls);
    # (6) reactive power balance
    @constraint(model, Cft'*qij .== Cg'*qg - Cd'*Qd+Cmg'*qm + Cd'*qls);
    # (7) KVL: Cft*v-2*diagm(R)*pij-2*diagm(X)*qij .<= (1-α)*bigM
    @constraint(model, Cft*v .- 2*diagm(R)*pij .- 2*diagm(X)*qij .<= (1 .- β)*bigM);
    # (8) KVL: Cft*v-2*diagm(R)*pij-2*diagm(X)*qij .>= -(1-α)*bigM
    @constraint(model, Cft*v .- 2*diagm(R)*pij .- 2*diagm(X)*qij .>= -(1 .- β)*bigM);
    # (9) pij .<= Smax.*β
    @constraint(model, pij.<=Smax.*β);
    # (10) pij .>= -Smax.*β
    @constraint(model, pij.>=-Smax.*β);
    # (11) qij .<= Smax.*β
    @constraint(model, qij.<=Smax.*β);
    # (12) qij .>= -Smax.*β
    @constraint(model, qij.>=-Smax.*β);
    # (13) -\sqrt{2} S_{ij}^{\max} \leq p_{ij} \pm q_{ij} \leq \sqrt{2} S_{ij}^{\max}
    @constraint(model, -sqrt(2)*Smax .<= pij .+ qij .<=sqrt(2)*Smax);
    @constraint(model, -sqrt(2)*Smax .<= pij .- qij .<=sqrt(2)*Smax);
    # (14) pls .== (Pd./Qd).*qls
    @constraint(model, pls.==(Pd./Qd).*qls);
    # (15) pmg .<= Betamg.*Pmgmax
    @constraint(model, pm.<=Betamg.*Pmgmax);
    # (16) pmg .>= Betamg.*Pmgmin
    @constraint(model, pm.>=Betamg.*Pmgmin);
    # (17) qmg .<= Betamg.*Qmgmax
    @constraint(model, qm.<=Betamg.*Qmgmax);
    # (18) qmg .>= Betamg.*Qmgmin
    @constraint(model, qm.>=Betamg.*Qmgmin);
    # (19) β = α .* ζ
    @constraint(model, β .>= α .- (1 .- ζ));
    @constraint(model, β .<= α);
    @constraint(model, β .<= ζ);
    # (20) topology constraints, \sum_{ij\in\mathcal{L}} \alpha_{ij} = |\mathcal{N}| - |\mathcal{G}| - \sum_{m\in\mathcal{M}} \beta_{m}
    @constraint(model, sum(α) == nb - ng - sum(Betamg));
    # (21) limit the voltage magnitude at the generator and mg bues
    @constraint(model, v[Int.(gen[:,GEN_BUS])].>= VMIN[Int.(gen[:,GEN_BUS])]);
    @constraint(model, v[Int.(gen[:,GEN_BUS])].<= VMAX[Int.(gen[:,GEN_BUS])]);

    @constraint(model, v[Int.(mg_buses)].>= VMIN[Int.(mg_buses)]);
    @constraint(model, v[Int.(mg_buses)].<= VMAX[Int.(mg_buses)]);

    # Solve model
    optimize!(model)

    # Extract results
    results = Dict(
        "Fij" => value.(Fij),
        "Fg" => value.(Fg),
        "Fm" => value.(Fm),
        "α" => value.(α),
        "β" => value.(β),
        "pls" => value.(pls),
        "qls" => value.(qls),
        "pij" => value.(pij),
        "qij" => value.(qij),
        "pg" => value.(pg),
        "qg" => value.(qg),
        "pm" => value.(pm),
        "qm" => value.(qm),
        "v" => value.(v),
        "objective" => objective_value(model),
        "status" => termination_status(model)
    )

    # Update mpc
    mpc_out = deepcopy(mpc)
    
    # Update bus data
    mpc_out["bus"][load_buses, PD] .-= results["pls"]*baseMVA
    mpc_out["bus"][load_buses, QD] .-= results["qls"]*baseMVA
    
    # Update generator data including microgrids
    gen_new = vcat(mpc_out["gen"], zeros(nmg, size(mpc_out["gen"],2)))
    for i = 1:nmg
        gen_new[ng+i, GEN_BUS] = mg_buses[i]
        gen_new[ng+i, PG] = results["pm"][i]*baseMVA
        gen_new[ng+i, QG] = results["qm"][i]*baseMVA
        gen_new[ng+i, VG] = sqrt(results["v"][mg_buses[i]])
        gen_new[ng+i, MBASE] = baseMVA
        gen_new[ng+i, GEN_STATUS] = 1
        gen_new[ng+i, PMAX] = Pmgmax[i]*baseMVA
        gen_new[ng+i, PMIN] = Pmgmin[i]*baseMVA
        gen_new[ng+i, QMAX] = Qmgmax[i]*baseMVA
        gen_new[ng+i, QMIN] = Qmgmin[i]*baseMVA
        mpc_out["bus"][mg_buses[i], BUS_TYPE] = REF
    end
    mpc_out["gen"] = gen_new

    # Update branch status
    for i = 1:nl
        mpc_out["branch"][i, BR_STATUS] = results["β"][i] > 0 ? 1 : 0
    end

    if verbose
        print_report(results, mpc_out, baseMVA, Pd, load_buses, Smax, Betamg, nb, ng)
    end

    return mpc_out, results
end

# Helper function to print detailed results
function print_report(results, mpc, baseMVA, Pd, load_buses, Smax, Betamg, nb, ng)
    println("\n=== 优化结果校核报告 ===")
    println("------------------------")
    # 1. 优化求解状态校核
    pls = results["pls"]
    pm = results["pm"]
    qm = results["qm"]
    v = results["v"]
    pg = results["pg"]
    println("\n1. 优化求解状态：")
    println("优化是否成功：", results["status"])
    println("目标函数值（负荷切除量）：", sum(pls)*baseMVA)

    # # 2. 电压约束校核
    println("\n2. 电压约束校核：")
    voltage_violations = findall(x -> x < 0.9^2 || x > 1.1^2, v)
    if isempty(voltage_violations)
        println("✓ 所有节点电压在允许范围内")
    else
        println("⚠ 存在电压越限节点：")
        for bus_idx in voltage_violations
            @printf("  节点 %d: V = %.4f p.u.\n", bus_idx, sqrt(v[bus_idx]))
        end
    end
    # 3. 功率平衡校核
    println("\n3. 功率平衡校核：")
    # 总发电功率
    total_gen = sum(pg)*baseMVA
    # 总微电网功率
    total_mg = sum(pm)*baseMVA
    # 总负荷功率
    total_load = sum(Pd .- pls)*baseMVA
    # 总网损
    total_loss = max(total_gen + total_mg - total_load, 0)

    # 4. 线路容量约束校核
    println("\n4. 线路容量约束校核：")
    pij = results["pij"]
    qij = results["qij"]
    line_violations = findall(x -> x > 0, (pij.^2 + qij.^2) .- Smax.^2)
    if isempty(line_violations)
        println("✓ 所有线路功率在额定容量范围内")
    else
        println("⚠ 存在越限线路：")
        for line_idx in line_violations
            apparent_power = sqrt(pij[line_idx]^2 + qij[line_idx]^2)
            @printf("  线路 %d: S = %.4f p.u., Smax = %.4f p.u.\n", 
                    line_idx, apparent_power, Smax[line_idx])
        end
    end

    # 5. 拓扑约束校核
    println("\n5. 拓扑约束校核：")
    α = results["α"]
    topology_constraint = sum(α) - (nb - ng - sum(Betamg))
    if abs(topology_constraint) < 1e-6
        println("✓ 拓扑约束满足")
    else
        println("⚠ 拓扑约束不满足，误差：", topology_constraint)
    end

    # 6. 负荷切除分析
    println("\n6. 负荷切除分析：")
    total_load_shedding = sum(pls)
    load_shedding_percentage = 100 * total_load_shedding / sum(Pd)
    @printf("总负荷切除量: %.4f p.u.\n", total_load_shedding)
    @printf("负荷切除比例: %.2f%%\n", load_shedding_percentage)

    # 7. 微电网参与情况
    println("\n7. 微电网参与情况：")
    for i in 1:nmg
        @printf("微电网 %d:\n", mg_buses[i])
        @printf("  有功功率: %.4f p.u.\n", pm[i])
        @printf("  无功功率: %.4f p.u.\n", qm[i])
    end

end

# Example usage
mpc = PowerFlow.case33bw()  # A default case
opt = PowerFlow.options()
mg_buses = [10; 14; 33]  # Microgrid buses
nmg = length(mg_buses)
Betamg = Int.(ones(nmg))
K = 5
branch = mpc["branch"]
(F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
     ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = idx_brch()

always_on_lines = findall(branch[:, BR_STATUS] .== 1);
failed_lines = randperm(length(always_on_lines))[1:K];
tie_lines  =  findall(branch[:, BR_STATUS] .== 0);
ζ = ones(size(branch, 1))
ζ[failed_lines] .= 0

mpc_out, results = network_reconfiguration(mpc, mg_buses, tie_lines, ζ, Betamg; verbose = false)

# Solve the power flow
# mpc_list, isolated = PowerFlow.extract_islands(mpc)
println("使用 $(Threads.nthreads()) 个线程进行计算")
mpc_list, isolated = PowerFlow.dc_preprocess(mpc_out, opt)

n_islands = length(mpc_list)
println("共提取出 $(n_islands) 个孤岛")

println("开始多线程潮流计算...")
t_start = time()

# 创建数组来存储每个岛屿的计算时间和线程ID
results_array = Vector{Any}(undef, n_islands)
island_times = zeros(n_islands)
thread_ids = zeros(Int, n_islands)

@threads for i in 1:n_islands
    local_start = time()
    thread_ids[i] = threadid()  # 记录当前线程ID
    results_array[i] = PowerFlow.runpf(mpc_list[i], opt)
    island_times[i] = time() - local_start
end

t_end = time()
elapsed = t_end - t_start

# 构造类似@timed返回的结果
results = (value=results_array, time=elapsed)
println("计算完成，总耗时: $(results.time) 秒")

# 输出每个岛屿的计算时间和线程ID
println("\n每个岛屿的计算详情:")
for i in 1:n_islands
    println("岛屿 $i: 线程 $(thread_ids[i]), 耗时 $(island_times[i]) 秒")
end

PowerFlow.process_result(results, isolated, "powerflow_report.txt")
# @time results_pf = PowerFlow.runpf(mpc_out, opt)

# println("\n=== Power Flow Results ===")
# println("------------------------")
# println("Success: ", results_pf["success"])