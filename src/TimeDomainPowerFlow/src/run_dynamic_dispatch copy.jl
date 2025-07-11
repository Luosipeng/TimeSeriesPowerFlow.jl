include("../../Solvers/types.jl")
include("../../Solvers/interior_point_method.jl")
include("../../Solvers/mips.jl")
using  LinearAlgebra

"""
    run_dynamic_dispatch_ipopt(new_jpc, Cld_ac, Cld_dc, loadAC_PD, loadAC_QD, loadDC_PD, genAC_PG, 
                              Cgen_ac, Cconv_ac, Cconv_dc, η_rec, η_inv, Cpv_ac, Cpv_dc, 
                              pv_ac_p_mw_ratio, pv_ac_p_mw, pv_max_p_mw, pv_max_p_mw_ratio, 
                              Cstorage_ac, ess_initial_soc, ess_max_soc, ess_min_soc, 
                              ess_power_capacity_mw, ess_energy_capacity_mwh, ess_efficiency, 
                              day_price_line, num_hours=24)

Solve the dynamic economic dispatch problem for a hybrid AC-DC power system using the interior point method.

# Arguments
- `new_jpc`: Renumbered hybrid power system data structure
- `Cld_ac`: AC load connection matrix
- `Cld_dc`: DC load connection matrix
- `loadAC_PD`: AC load active power demand (MW) over time
- `loadAC_QD`: AC load reactive power demand (MVar) over time
- `loadDC_PD`: DC load power demand (MW) over time
- `genAC_PG`: AC generator active power output (MW)
- `Cgen_ac`: AC generator connection matrix
- `Cconv_ac`: AC side converter connection matrix
- `Cconv_dc`: DC side converter connection matrix
- `η_rec`: Rectifier efficiency (AC to DC conversion)
- `η_inv`: Inverter efficiency (DC to AC conversion)
- `Cpv_ac`: AC PV system connection matrix
- `Cpv_dc`: DC PV system connection matrix
- `pv_ac_p_mw_ratio`: AC PV power output ratio over time
- `pv_ac_p_mw`: AC PV system maximum power output (MW)
- `pv_max_p_mw`: Maximum power output of each PV system (MW)
- `pv_max_p_mw_ratio`: PV power output ratio over time
- `Cstorage_ac`: Energy storage system connection matrix
- `ess_initial_soc`: Initial state of charge for energy storage systems
- `ess_max_soc`: Maximum state of charge for energy storage systems
- `ess_min_soc`: Minimum state of charge for energy storage systems
- `ess_power_capacity_mw`: Power capacity of energy storage systems (MW)
- `ess_energy_capacity_mwh`: Energy capacity of energy storage systems (MWh)
- `ess_efficiency`: Energy storage charging/discharging efficiency
- `day_price_line`: Electricity price data over time
- `num_hours`: Number of hours in the optimization horizon (default: 24)

# Description
This function formulates and solves the dynamic economic dispatch problem for a hybrid AC-DC power system
using an interior point optimization method. The objective is to minimize the total generation cost while
satisfying power balance constraints, converter operation constraints, and energy storage operation constraints.

The optimization variables include:
- Branch power flows (Pij)
- Generator power outputs (Pgen)
- Converter power flows (Pij_inv, Pij_rec)
- PV system power outputs (P_pv_mw)
- Energy storage state of charge (soc)
- Energy storage charging/discharging power (ess_charge, ess_discharge)
- Energy storage operation mode (ess_mode)

Key constraints include:
- Power balance at each node for each time period
- Converter mutual exclusivity (cannot operate in both inverter and rectifier modes simultaneously)
- Energy storage state of charge evolution
- Energy storage charging/discharging mutual exclusivity
- Initial and final state of charge requirements

The function returns detailed optimization results including the optimal solution, objective value,
convergence status, and constraint violations.
"""

function run_dynamic_dispatch_ipopt(new_jpc,
        Cld_ac, Cld_dc, 
        loadAC_PD, loadAC_QD,
        loadDC_PD,
        genAC_PG,
        Cgen_ac, 
        Cconv_ac, Cconv_dc, 
        η_rec, η_inv,
        Cpv_ac,
        Cpv_dc,
        pv_ac_p_mw_ratio,
        pv_ac_p_mw,
        pv_max_p_mw,
        pv_max_p_mw_ratio,
        Cstorage_ac,
        ess_initial_soc,
        ess_max_soc,
        ess_min_soc,
        ess_power_capacity_mw,
        ess_energy_capacity_mwh,
        ess_efficiency,
        day_price_line,
        num_hours=24)

    # 提取连接信息
    branchAC = new_jpc.branchAC
    branchDC = new_jpc.branchDC
    converter = new_jpc.converter

    rac = branchAC[:,BR_R]
    rdc = branchDC[:,BR_R]
    rconv = zeros(size(converter, 1))
    
    # 提取母线数量
    n_nodes = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)
    A, branch_data = TimeSeriesPowerFlow.build_incidence_matrix_td(n_nodes, branchAC, branchDC, converter)

    # 初始化
    nbr = size(A, 1)  # 支路数
    nc = size(converter, 1)  # 变流器数
    ng = size(genAC_PG, 1)  # 交流发电机数
    npv = size(new_jpc.pv, 1)  # 光伏设备数
    ns = size(new_jpc.storage, 1)  # 储能设备数

    # 搜寻支路索引
    dcbranch_indices = findall(x -> x[3] == 2, branch_data)
    converter_indices = findall(x -> x[3] == 3, branch_data)
    acbranch_indices = findall(x -> x[3] == 1, branch_data)

    r = zeros(nbr)  # 初始化支路电阻
    r[acbranch_indices] = rac  # 交流支路电阻
    r[dcbranch_indices] = rdc  # 直流支路电阻
    r[converter_indices] = rconv  # 变流器电阻

    # 定义优化问题的变量数量和约束数量
    n_vars = nbr*num_hours + ng*num_hours + 2*nc*num_hours + npv*num_hours + ns*num_hours*3 + ns*num_hours
    
    # 变量索引映射
    indices = create_reorganized_indices(nbr, ng, nc, npv, ns, num_hours)
    
    # 变量下界和上界
    lb = fill(-Inf, n_vars)
    ub = fill(Inf, n_vars)
    
    # 设置支路功率流的界限 (Pij)
    # 这里可以根据支路容量设置限制，暂时不设置
    
    # 设置发电机功率的界限 (Pgen)
    # 假设发电机有最小和最大输出限制
    for i in 1:ng
        for t in 1:num_hours
            idx = indices["Pgen"][(i, t)]
            lb[idx] = -Inf  # 允许负值（如果需要）
            ub[idx] = Inf   # 无上限
        end
    end
    
    # 设置变流器功率的界限 (Pij_inv, Pij_rec)
    for i in 1:nc
        for t in 1:num_hours
            inv_idx = indices["Pij_inv"][(i, t)]
            rec_idx = indices["Pij_rec"][(i, t)]
            lb[inv_idx] = 0.0  # 逆变器功率非负
            lb[rec_idx] = 0.0  # 整流器功率非负
        end
    end
    
    # 设置光伏功率的界限 (P_pv_mw)
    for i in 1:npv
        for t in 1:num_hours
            idx = indices["P_pv_mw"][(i, t)]
            lb[idx] = 0.0  # 光伏功率非负
            ub[idx] = pv_max_p_mw[i] * pv_max_p_mw_ratio[i, t]  # 最大光伏功率
        end
    end
    
    # 设置储能SOC的界限 (soc)
    for i in 1:ns
        for t in 1:num_hours
            idx = indices["soc"][(i, t)]
            lb[idx] = ess_min_soc[i]  # 最小SOC
            ub[idx] = ess_max_soc[i]  # 最大SOC
        end
    end
    
    # 设置储能充放电功率的界限 (ess_charge, ess_discharge)
    for i in 1:ns
        for t in 1:num_hours
            charge_idx = indices["ess_charge"][(i, t)]
            discharge_idx = indices["ess_discharge"][(i, t)]
            mode_idx = indices["ess_mode"][(i, t)]
            
            lb[charge_idx] = 0.0  # 充电功率非负
            ub[charge_idx] = ess_power_capacity_mw[i]  # 最大充电功率
            
            lb[discharge_idx] = 0.0  # 放电功率非负
            ub[discharge_idx] = ess_power_capacity_mw[i]  # 最大放电功率
            
            lb[mode_idx] = 0.0  # 模式变量下界
            ub[mode_idx] = 1.0  # 模式变量上界
        end
    end
    
    # 设置换流器功率约束 (Pij[converter_indices] = 0)
    for i in converter_indices
        for t in 1:num_hours
            idx = indices["Pij"][(i, t)]
            lb[idx] = 0.0
            ub[idx] = 0.0
        end
    end
    
    # # 初始点
    # x0 = zeros(n_vars)
    
    # # 智能初始化
    # # 发电机初始功率
    # for i in 1:ng
    #     for t in 1:num_hours
    #         idx = indices["Pgen"][(i, t)]
    #         x0[idx] = 0.5 * (lb[idx] + (ub[idx] == Inf ? lb[idx] + 100.0 : ub[idx]))  # 在界限内取中点
    #     end
    # end
    
    # # 光伏初始功率
    # for i in 1:npv
    #     for t in 1:num_hours
    #         idx = indices["P_pv_mw"][(i, t)]
    #         x0[idx] = 0.5 * pv_max_p_mw[i] * pv_max_p_mw_ratio[i, t]  # 初始为最大功率的一半
    #     end
    # end
    
    # # 储能初始SOC
    # for i in 1:ns
    #     for t in 1:num_hours
    #         idx = indices["soc"][(i, t)]
    #         x0[idx] = ess_initial_soc[i]  # 初始SOC
    #     end
    # end
    x0 = set_warm_start_initial_point()
    
    # 目标函数: 最小化发电成本
    function f_obj(x)
        cost = 0.0
        for i in 1:ng
            for t in 1:num_hours
                idx = indices["Pgen"][(i, t)]
                cost += day_price_line[t, 2] * x[idx]
            end
        end
        return cost
    end
    
    # 目标函数梯度
    function ∇f_obj(x)
        grad = zeros(n_vars)
        for i in 1:ng
            for t in 1:num_hours
                idx = indices["Pgen"][(i, t)]
                grad[idx] = day_price_line[t, 2]
            end
        end
        return grad
    end
    
    # 目标函数Hessian (零矩阵，因为目标函数是线性的)
    function ∇2f_obj(x)
        return spzeros(n_vars, n_vars)
    end
    
    # 等式约束: 功率平衡约束
    function g(x)
        constraints = Float64[]
        
        # 对每个时间步
        for t in 1:num_hours
            # 提取当前时间步的变量
            Pij_t = [x[indices["Pij"][(i, t)]] for i in 1:nbr]
            Pgen_t = [x[indices["Pgen"][(i, t)]] for i in 1:ng]
            Pij_inv_t = [x[indices["Pij_inv"][(i, t)]] for i in 1:nc]
            Pij_rec_t = [x[indices["Pij_rec"][(i, t)]] for i in 1:nc]
            P_pv_mw_t = [x[indices["P_pv_mw"][(i, t)]] for i in 1:npv]
            ess_charge_t = [x[indices["ess_charge"][(i, t)]] for i in 1:ns]
            ess_discharge_t = [x[indices["ess_discharge"][(i, t)]] for i in 1:ns]
            
            # 功率平衡约束
            power_balance = A' * Pij_t + 
                           Cld_ac * loadAC_PD[:, t] + 
                           Cld_dc * loadDC_PD[:, t] - 
                           Cgen_ac * Pgen_t - 
                           Cpv_ac * (pv_ac_p_mw .* pv_ac_p_mw_ratio[:, t]) + 
                           Cconv_ac * (Pij_inv_t - η_rec .* Pij_rec_t) + 
                           Cconv_dc * (Pij_rec_t - η_inv .* Pij_inv_t) - 
                           Cpv_dc * P_pv_mw_t + 
                           Cstorage_ac * (ess_charge_t - ess_discharge_t)
            
            append!(constraints, power_balance)
            
            # 换流器互斥约束: Pij_inv .* Pij_rec = 0
            for i in 1:nc
                inv_idx = indices["Pij_inv"][(i, t)]
                rec_idx = indices["Pij_rec"][(i, t)]
                push!(constraints, x[inv_idx] * x[rec_idx])
            end
            if t == 1
                # 初始SOC约束
                for i in 1:ns
                    soc_idx_1 = indices["soc"][(i, 1)]
                    charge_idx_1 = indices["ess_charge"][(i, 1)]
                    discharge_idx_1 = indices["ess_discharge"][(i, 1)]
                    
                    soc_constraint = x[soc_idx_1] - ess_initial_soc[i] + 
                                    x[discharge_idx_1] / ess_efficiency[i] / ess_energy_capacity_mwh[i] - 
                                    x[charge_idx_1] * ess_efficiency[i] / ess_energy_capacity_mwh[i]
                    push!(constraints, soc_constraint)
                end
            else
                # SOC随时间的演变
                for i in 1:ns
                    soc_idx_t = indices["soc"][(i, t)]
                    soc_idx_prev = indices["soc"][(i, t-1)]
                    charge_idx_t = indices["ess_charge"][(i, t)]
                    discharge_idx_t = indices["ess_discharge"][(i, t)]
                    
                    soc_constraint = x[soc_idx_t] - x[soc_idx_prev] + 
                                    x[discharge_idx_t] / ess_efficiency[i] / ess_energy_capacity_mwh[i] - 
                                    x[charge_idx_t] * ess_efficiency[i] / ess_energy_capacity_mwh[i]
                    push!(constraints, soc_constraint)
                end
            end
        end
        # 循环约束：最终SOC等于初始SOC
        for i in 1:ns
            soc_idx_final = indices["soc"][(i, num_hours)]
            final_soc_constraint = x[soc_idx_final] - ess_initial_soc[i]
            push!(constraints, final_soc_constraint)
        end
        
        return constraints
    end
    
    # 修改等式约束雅可比矩阵的计算函数
    function ∇g(x)
        # 首先计算约束向量，以确保维度一致
        constraints = g(x)
        n_constraints = length(constraints)
        n_vars = length(x)
        
        # 确定各类约束的数量
        n_nodes = size(A, 2)  # 节点数
        n_power_balance = n_nodes * num_hours  # 功率平衡约束
        n_converter_mutual = nc * num_hours    # 换流器互斥约束
        n_soc_initial = ns                     # 初始SOC约束
        n_soc_evolution = ns * (num_hours - 1) # SOC演变约束
        n_soc_final = ns                       # 最终SOC约束
        
        # 验证约束总数
        expected_constraints = n_power_balance + n_converter_mutual + n_soc_initial + n_soc_evolution + n_soc_final
        if n_constraints != expected_constraints
            error("约束函数返回的约束数量 $(n_constraints) 与预期数量 $(expected_constraints) 不匹配")
        end
        
        Jacobian = zeros(0, n_vars)  # 初始化雅可比矩阵
        
        # 按时间步处理功率平衡和换流器约束
        for t in 1:num_hours
            Jac_balance = zeros(n_nodes, n_vars)
            Jac_offset = zeros(nc, n_vars)
            
            # 1. 功率平衡约束的雅可比矩阵
            # 对于每个节点和每个变量类型，正确填充雅可比矩阵
            
            # Pij的贡献 (A'的转置)
            for br in 1:nbr
                col_idx = indices["Pij"][(br, t)]
                Jac_balance[:, col_idx] = A[br, :]
            end
            
            # Pgen的贡献 (-Cgen_ac)
            for g in 1:ng
                col_idx = indices["Pgen"][(g, t)]
                Jac_balance[:, col_idx] = -Cgen_ac[:, g]
            end
            
            # Pij_inv的贡献 (Cconv_ac - η_inv.*Cconv_dc)
            for c in 1:nc
                col_idx = indices["Pij_inv"][(c, t)]
                Jac_balance[:, col_idx] = Cconv_ac[:, c] - η_inv[c] * Cconv_dc[:, c]
            end
            
            # Pij_rec的贡献 (Cconv_dc - η_rec.*Cconv_ac)
            for c in 1:nc
                col_idx = indices["Pij_rec"][(c, t)]
                Jac_balance[:, col_idx] = Cconv_dc[:, c] - η_rec[c] * Cconv_ac[:, c]
            end
            
            # P_pv_mw的贡献 (-Cpv_dc)
            for pv in 1:npv
                col_idx = indices["P_pv_mw"][(pv, t)]
                Jac_balance[:, col_idx] = -Cpv_dc[:, pv]
            end
            
            # ess_charge的贡献 (Cstorage_ac)
            for s in 1:ns
                col_idx = indices["ess_charge"][(s, t)]
                Jac_balance[:, col_idx] = Cstorage_ac[:, s]
            end
            
            # ess_discharge的贡献 (-Cstorage_ac)
            for s in 1:ns
                col_idx = indices["ess_discharge"][(s, t)]
                Jac_balance[:, col_idx] = -Cstorage_ac[:, s]
            end
            
            # 2. 换流器互斥约束的雅可比矩阵
            for i in 1:nc
                Jac_offset[i, indices["Pij_inv"][(i, t)]] = x[indices["Pij_rec"][(i, t)]]
                Jac_offset[i, indices["Pij_rec"][(i, t)]] = x[indices["Pij_inv"][(i, t)]]
            end
            
            Jacobian = [Jacobian; Jac_balance; Jac_offset]
            
            # 3. SOC初始约束的雅可比矩阵
            if t == 1
                Jac_soc_initial = zeros(ns, n_vars)
                for s in 1:ns
                    Jac_soc_initial[s, indices["soc"][(s, 1)]] = 1.0
                    Jac_soc_initial[s, indices["ess_discharge"][(s, 1)]] = 1.0 / (ess_efficiency[s] * ess_energy_capacity_mwh[s])
                    Jac_soc_initial[s, indices["ess_charge"][(s, 1)]] = -ess_efficiency[s] / ess_energy_capacity_mwh[s]
                end
                Jacobian = [Jacobian; Jac_soc_initial]
            end
            
            # 4. SOC演变约束的雅可比矩阵
            if t > 1
                Jac_soc_evolution = zeros(ns, n_vars)
                for s in 1:ns
                    Jac_soc_evolution[s, indices["soc"][(s, t)]] = 1.0
                    Jac_soc_evolution[s, indices["soc"][(s, t-1)]] = -1.0
                    Jac_soc_evolution[s, indices["ess_discharge"][(s, t)]] = 1.0 / (ess_efficiency[s] * ess_energy_capacity_mwh[s])
                    Jac_soc_evolution[s, indices["ess_charge"][(s, t)]] = -ess_efficiency[s] / ess_energy_capacity_mwh[s]
                end
                Jacobian = [Jacobian; Jac_soc_evolution]
            end
        end
        
        # 5. 储能最终SOC约束的雅可比矩阵
        Jac_final_soc = zeros(ns, n_vars)
        for s in 1:ns
            Jac_final_soc[s, indices["soc"][(s, num_hours)]] = 1.0
        end
        Jacobian = [Jacobian; Jac_final_soc]
        
        return Jacobian'
    end
    
    # 不等式约束: 储能充放电互斥约束
    function h(x)
        constraints = Float64[]
        
        # 储能充放电互斥约束
        for t in 1:num_hours
            for i in 1:ns
                discharge_idx = indices["ess_discharge"][(i, t)]
                mode_idx = indices["ess_mode"][(i, t)]
                charge_idx = indices["ess_charge"][(i, t)]
            
                # ess_charge <= ess_power_capacity_mw * ess_mode
                push!(constraints, x[charge_idx] - ess_power_capacity_mw[i] * x[mode_idx])

                # ess_discharge <= ess_power_capacity_mw * (1 - ess_mode)
                push!(constraints, x[discharge_idx] - ess_power_capacity_mw[i] * (1 - x[mode_idx]))
            end
        end
        
        return constraints
    end
    
    # 修改不等式约束雅可比矩阵的计算函数
    function ∇h(x)
        # 首先计算约束向量，以确保维度一致
        constraints = h(x)
        n_constraints = length(constraints)
        n_vars = length(x)
        
        # 确定不等式约束数量
        n_ess_discharge_constraints = ns * num_hours  # ess_discharge <= power_capacity * (1 - mode)
        n_ess_charge_constraints = ns * num_hours     # ess_charge <= power_capacity * mode
        
        # 验证约束总数
        expected_constraints = n_ess_discharge_constraints + n_ess_charge_constraints
        if n_constraints != expected_constraints
            error("不等式约束函数返回的约束数量 $(n_constraints) 与预期数量 $(expected_constraints) 不匹配")
        end
        
        Jacobian = zeros(0, n_vars)  # 初始化雅可比矩阵
        for t in 1:num_hours
            Jac_charge = zeros(ns, n_vars)
            Jac_discharge = zeros(ns, n_vars)
            for s in 1:ns
                # 1. 储能充电约束的雅可比矩阵
                Jac_charge[s, indices["ess_charge"][(s, t)]] = 1.0  # ∂/∂charge = 1
                Jac_charge[s, indices["ess_mode"][(s, t)]] = -ess_power_capacity_mw[s]  # ∂/∂mode = -power_capacity
                # 2. 储能放电约束的雅可比矩阵
                Jac_discharge[s, indices["ess_discharge"][(s, t)]] = 1.0  # ∂/∂discharge = 1
                Jac_discharge[s, indices["ess_mode"][(s, t)]] = ess_power_capacity_mw[s]  # ∂/∂mode = power_capacity
            end
            Jacobian = [Jacobian; vcat(Jac_charge, Jac_discharge)]
        end
        
       
        return Jacobian'
    end

    
    function Lxx(x, λ, μ)
        
        hessians = zeros(n_vars, n_vars)  # 初始化Hessian矩阵
        
        # 计算各类约束在λ中的起始索引
        n_power_balance = n_nodes * num_hours
        n_converter_mutual = nc * num_hours
        n_soc_initial = ns
        n_soc_evolution = ns * (num_hours - 1)
        
        # 换流器互斥约束在λ中的起始索引
        converter_constraint_start = n_power_balance + 1
        
        for t in 1:num_hours
            for c in 1:nc
                # 当前换流器互斥约束在λ中的索引
                constraint_idx = converter_constraint_start + (t-1) * nc + (c-1)
                
                # 检查索引是否有效
                if constraint_idx <= length(λ)
                    λ_val = λ[constraint_idx]
                    
                    # 添加Hessian贡献
                    inv_idx = indices["Pij_inv"][(c, t)]
                    rec_idx = indices["Pij_rec"][(c, t)]
                    
                    hessians[inv_idx, rec_idx] += λ_val * 1.0
                    hessians[rec_idx, inv_idx] += λ_val * 1.0
                end
            end
        end
        
        return hessians
    end

    
    # 拉格朗日函数梯度
    function Lx(x, λ, μ)
        # 计算目标函数梯度
        grad_f = ∇f_obj(x)
        
        # 计算等式约束的贡献
        constraints_g = g(x)
        J_g = ∇g(x)  # 现在这是 (变量数 × 约束数) 矩阵
        
        # 计算不等式约束的贡献
        constraints_h = h(x)
        J_h = ∇h(x)  # 现在这是 (变量数 × 约束数) 矩阵
        
        # 确保λ和μ的维度正确
        if length(λ) != length(constraints_g)
            error("拉格朗日乘子λ的维度 $(length(λ)) 与等式约束数量 $(length(constraints_g)) 不匹配")
        end
        
        if length(μ) != length(constraints_h)
            error("拉格朗日乘子μ的维度 $(length(μ)) 与不等式约束数量 $(length(constraints_h)) 不匹配")
        end
        
        # 由于雅可比矩阵已经是转置的，直接相乘即可
        eq_contribution = J_g * λ      # (变量数 × 约束数) × (约束数 × 1) = (变量数 × 1)
        ineq_contribution = J_h * μ    # (变量数 × 约束数) × (约束数 × 1) = (变量数 × 1)
        
        return grad_f + eq_contribution + ineq_contribution
    end
    
    # 在求解前添加调试信息
    println("变量数量: ", n_vars)
    println("等式约束数量: ", length(g(x0)))
    println("不等式约束数量: ", length(h(x0)))
    println("等式约束雅可比矩阵维度: ", size(∇g(x0)))
    println("不等式约束雅可比矩阵维度: ", size(∇h(x0)))
    
    # 设置NonConvexOPT问题结构
    nonlinear = NonConvexOPT(
        f_obj,      # 目标函数
        ∇f_obj,     # 目标函数梯度
        ∇2f_obj,    # 目标函数Hessian
        g,          # 等式约束 g(x) = 0
        ∇g,         # 等式约束梯度
        h,          # 不等式约束 h(x) <= 0
        ∇h,         # 不等式约束梯度
        Lx,         # 拉格朗日函数梯度
        Lxx,        # 拉格朗日函数Hessian
        x0          # 初始点
    )
    
    # 求解问题
    println("开始求解内部点法...")
    
    ipopt_result = interior_point_method(nonlinear, IPM(1e-6, 100, 1e-6, 1e-6, true, 0.99995))
    mips_result = mips(nonlinear, IPM(1e-6, 100, 1e-6, 1e-6, true, 0.99995))
    # Print detailed comparison
    println("\n" * "="^60)
    println("OPTIMIZATION RESULTS COMPARISON")
    println("="^60)

    println("\nIPM Results:")
    # println("  Solution: ", result.x)
    println("  Objective: ", ipopt_result.obj)
    println("  Exit flag: ", ipopt_result.eflag)
    println("  Iterations: ", ipopt_result.iterations)
    println("  Constraint violations (g): ", norm(g(ipopt_result.x), Inf))
    println("  Constraint violations (h): ", maximum(max.(h(ipopt_result.x), 0.0)))
    println("  History ", ipopt_result.hist.x_record)
    println("  History ", ipopt_result.hist.obj_record)
    println("  History ", ipopt_result.hist.feascond_record)
    println("  History ", ipopt_result.hist.gradcond_record)
    println("  History ", ipopt_result.hist.compcond_record)
    println("  History ", ipopt_result.hist.costcond_record)
end


function create_reorganized_indices(nbr, ng, nc, npv, ns, num_hours)
  # 计算每个时间步内的变量数量
  vars_per_timestep = nbr + ng + nc + nc + npv + ns + ns + ns + ns
  # 对应: Pij + Pgen + Pij_inv + Pij_rec + P_pv_mw + soc + ess_charge + ess_discharge + ess_mode
  
  # 总变量数
  n_vars = vars_per_timestep * num_hours
  
  # 初始化索引字典
  indices = Dict(
      "n_vars" => n_vars,
      "vars_per_timestep" => vars_per_timestep,
      "Pij" => Dict{Tuple{Int, Int}, Int}(),
      "Pgen" => Dict{Tuple{Int, Int}, Int}(),
      "Pij_inv" => Dict{Tuple{Int, Int}, Int}(),
      "Pij_rec" => Dict{Tuple{Int, Int}, Int}(),
      "P_pv_mw" => Dict{Tuple{Int, Int}, Int}(),
      "soc" => Dict{Tuple{Int, Int}, Int}(),
      "ess_charge" => Dict{Tuple{Int, Int}, Int}(),
      "ess_discharge" => Dict{Tuple{Int, Int}, Int}(),
      "ess_mode" => Dict{Tuple{Int, Int}, Int}()
  )
  
  # 为每个时间步和变量类型计算索引
  for t in 1:num_hours
      # 当前时间步的起始索引
      time_offset = (t - 1) * vars_per_timestep
      
      # 在每个时间步内的偏移量
      current_offset = 0
      
      # Pij: 支路功率流 (索引 1 到 nbr)
      for i in 1:nbr
          indices["Pij"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += nbr
      
      # Pgen: 发电机功率 (索引 nbr+1 到 nbr+ng)
      for i in 1:ng
          indices["Pgen"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += ng
      
      # Pij_inv: 变流器逆变功率 (索引 nbr+ng+1 到 nbr+ng+nc)
      for i in 1:nc
          indices["Pij_inv"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += nc
      
      # Pij_rec: 变流器整流功率 (索引 nbr+ng+nc+1 到 nbr+ng+2*nc)
      for i in 1:nc
          indices["Pij_rec"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += nc
      
      # P_pv_mw: 光伏功率 (索引 nbr+ng+2*nc+1 到 nbr+ng+2*nc+npv)
      for i in 1:npv
          indices["P_pv_mw"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += npv
      
      # soc: 储能SOC (索引 nbr+ng+2*nc+npv+1 到 nbr+ng+2*nc+npv+ns)
      for i in 1:ns
          indices["soc"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += ns
      
      # ess_charge: 储能充电功率 (索引 nbr+ng+2*nc+npv+ns+1 到 nbr+ng+2*nc+npv+2*ns)
      for i in 1:ns
          indices["ess_charge"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += ns
      
      # ess_discharge: 储能放电功率 (索引 nbr+ng+2*nc+npv+2*ns+1 到 nbr+ng+2*nc+npv+3*ns)
      for i in 1:ns
          indices["ess_discharge"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += ns
      
      # ess_mode: 储能模式 (索引 nbr+ng+2*nc+npv+3*ns+1 到 nbr+ng+2*nc+npv+4*ns)
      for i in 1:ns
          indices["ess_mode"][(i, t)] = time_offset + current_offset + i
      end
      current_offset += ns
  end
  
  return indices
end

function set_warm_start_initial_point()
    x0 = zeros(n_vars)
    
    # 之前求解的结果数据
    Pgen_values = [2.9217006660653846, 2.6807916001104717, 2.437251907027197, 2.3290843435344986, 2.2309101319499924, 2.087252365238023, 1.9185509207746714, 0.6072760699579337, 0.28966877331729396, -1.0535419649316617, -2.035015429514435, -1.357153554417465, -1.199062246607309, -0.7577326879584394, -1.4033330097503667, 0.13815507247823008, 1.4420888807825594, 2.448719174718103, 2.8606868355603914, 2.5141104044459257, 3.1658644216804808, 3.436383235771799, 3.6250743342245397, 3.1736250899906677]
    
    Pij_inv_values = [0.13952021459561664, 0.13331314864070562, 0.13099580555743043, 0.1317727420647328, 0.13589688048022514, 0.1449552137682554, 0.16433170985490347, 0.22291798294616783, 1.3450352996444596e-7, 1.4133057077104394e-7, 1.4793779689584358e-7, 5.3744276660479e-7, 6.209249234411527e-7, 7.372137928550491e-7, 9.163378672092868e-7, 1.246276466373501e-6, 2.2032027918439455e-6, 1.4120062912550012e-34, 2.4746510517341903e-34, 6.72492418504592e-35, -6.498998388541293e-36, -2.9302212613466825e-35, 0.2627280827547721, 0.20023458852090126]
    
    Pij_rec_values = [7.612408666839252e-41, -6.234934287501052e-41, 9.32596260252787e-41, 1.8371857553140555e-40, -5.774942113044937e-44, 3.995236672657833e-41, 2.0086219888834206e-40, -3.813744405067803e-40, -1.262286725834044e-34, -1.2766039478985325e-34, -1.260858084617512e-35, 9.48130958251469e-33, -3.6918791210406676e-32, -3.964025202712717e-35, -4.481376616680721e-33, 1.0110349393720196e-32, -4.670511096110994e-33, 1.307501849166111e-6, 0.6074974065659732, 0.6074973855820428, 1.033099208186765e-6, 1.5729977420428322e-6, -4.0245099092177075e-40, -8.263269681081897e-36]
    
    ess_charge_values = [0.12556963648770467, 0.11998327712962298, 0.11789766835520195, 0.11859691121159577, 0.12230863578462117, 0.13046113574198565, 0.1478999822165434, 0.20062762799102946, 2.635980496649934e-6, 2.4157445309411417e-6, 2.250854853256369e-6, 2.3188712338418776e-6, 2.2483615783651785e-6, 2.2088282589864653e-6, 2.2178971002478016e-6, 2.3385950036514587e-6, 2.9585741314982483e-6, 6.635478898637536e-7, 4.663653611951111e-7, 4.6636536119623193e-7, 6.749245632905391e-7, 6.534875609452414e-7, 0.23645591040957023, 0.18021176560040847]
    
    ess_discharge_values = [1.443351649700378e-6, 1.443352987925805e-6, 1.4433535145532933e-6, 1.443353336253817e-6, 1.443352418558582e-6, 1.4433505557827594e-6, 1.4433471302796916e-6, 1.4433394784332278e-6, 2.514927319681933e-6, 2.2885470172472015e-6, 2.1177108360501096e-6, 1.8351727438975664e-6, 1.689529147268141e-6, 1.545335845416921e-6, 1.3931930197594436e-6, 1.2169461839153075e-6, 9.756916188386968e-7, 1.971049739029865e-6, 0.6074978729313344, 0.607497851947404, 1.7080237714773043e-6, 2.2264853029880737e-6, 6.359302753041705e-7, 6.359315973246111e-7]
    
    soc_values = [0.37534071274325265, 0.447329609870665, 0.5180671417330346, 0.5892242193093726, 0.6626083316302056, 0.7408839439268372, 0.8296228641107408, 0.9499983717650041, 0.9499980904441764, 0.9499978446708822, 0.9499976265091008, 0.9499976584446234, 0.9499977559584983, 0.949997936562235, 0.9499982353056656, 0.9499987370210501, 0.9499997894309965, 0.9499987275228867, 0.5000005829485223, 0.05000245391781017, 0.05000159366975446, 0.05000033651391844, 0.19187341170019737, 0.3]
    
    # 设置发电机功率
    for i in 1:ng
        for t in 1:num_hours
            idx = indices["Pgen"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(Pgen_values)
                x0[idx] = Pgen_values[linear_idx]
            end
        end
    end
    
    # 设置变流器逆变功率
    for i in 1:nc
        for t in 1:num_hours
            idx = indices["Pij_inv"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(Pij_inv_values)
                x0[idx] = Pij_inv_values[linear_idx]
            end
        end
    end
    
    # 设置变流器整流功率
    for i in 1:nc
        for t in 1:num_hours
            idx = indices["Pij_rec"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(Pij_rec_values)
                x0[idx] = Pij_rec_values[linear_idx]
            end
        end
    end
    
    # 设置储能充电功率
    for i in 1:ns
        for t in 1:num_hours
            idx = indices["ess_charge"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(ess_charge_values)
                x0[idx] = ess_charge_values[linear_idx]
            end
        end
    end
    
    # 设置储能放电功率
    for i in 1:ns
        for t in 1:num_hours
            idx = indices["ess_discharge"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(ess_discharge_values)
                x0[idx] = ess_discharge_values[linear_idx]
            end
        end
    end
    
    # 设置储能SOC
    for i in 1:ns
        for t in 1:num_hours
            idx = indices["soc"][(i, t)]
            linear_idx = (i-1) * num_hours + t
            if linear_idx <= length(soc_values)
                x0[idx] = soc_values[linear_idx]
            end
        end
    end
    
    # 其他变量保持默认初始化
    # 支路功率 - 初始为0或小值
    for i in 1:nbr
        for t in 1:num_hours
            idx = indices["Pij"][(i, t)]
            x0[idx] = 0.0
        end
    end
    
    # 光伏功率 - 根据可用功率设置
    for i in 1:npv
        for t in 1:num_hours
            idx = indices["P_pv_mw"][(i, t)]
            x0[idx] = 0.5 * pv_max_p_mw[i] * pv_max_p_mw_ratio[i, t]
        end
    end
    
    # 储能模式变量
    for i in 1:ns
        for t in 1:num_hours
            idx = indices["ess_mode"][(i, t)]
            # 根据充放电状态设置模式
            charge_idx = indices["ess_charge"][(i, t)]
            discharge_idx = indices["ess_discharge"][(i, t)]
            
            if x0[charge_idx] > x0[discharge_idx]
                x0[idx] = 0.8  # 充电模式
            elseif x0[discharge_idx] > x0[charge_idx]
                x0[idx] = 0.2  # 放电模式
            else
                x0[idx] = 0.5  # 中性模式
            end
        end
    end
    
    return x0
end