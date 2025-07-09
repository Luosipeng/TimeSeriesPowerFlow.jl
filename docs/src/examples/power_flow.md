# PowerFlow 模块示例

本文档提供了 PowerFlow 模块的详细使用示例，展示了各种潮流计算算法的应用、高级选项配置以及结果分析方法。这些示例旨在帮助用户深入理解 PowerFlow 模块的核心功能。

## 基础潮流计算

### 牛顿-拉夫森法潮流计算

牛顿-拉夫森法是最常用的潮流计算方法，适用于大多数电力系统：

```julia
using PowerFlow

# 加载标准测试系统
case = load_case("case14.m")

# 使用牛顿-拉夫森法计算潮流
opt = Dict("pf_alg" => "NR", "verbose" => 1)
results = runpf(case, opt)

# 显示收敛过程
println("迭代次数: $(results.iterations)")
println("最终误差: $(results.max_mismatch)")

# 查看计算结果
println("系统总发电: $(sum(results.gen[:, PG])) MW")
println("系统总负荷: $(sum(results.bus[:, PD])) MW")
println("系统总损耗: $(sum(results.branch[:, PL])) MW")
```

### 快速解耦潮流计算

快速解耦法在某些情况下可以提供更快的计算速度：

```julia
# 使用快速解耦法计算潮流
opt = Dict("pf_alg" => "FDXB", "verbose" => 1)
fdxb_results = runpf(case, opt)

# 比较与牛顿-拉夫森法的结果
println("牛顿-拉夫森法迭代次数: $(results.iterations)")
println("快速解耦法迭代次数: $(fdxb_results.iterations)")

# 比较计算精度
nr_losses = sum(results.branch[:, PL])
fdxb_losses = sum(fdxb_results.branch[:, PL])
println("牛顿-拉夫森法计算的系统损耗: $nr_losses MW")
println("快速解耦法计算的系统损耗: $fdxb_losses MW")
println("相对误差: $(abs(nr_losses - fdxb_losses) / nr_losses * 100)%")
```

### 高斯-赛德尔法潮流计算

高斯-赛德尔法适用于某些特定类型的系统：

```julia
# 使用高斯-赛德尔法计算潮流
opt = Dict("pf_alg" => "GS", "verbose" => 1, "max_it" => 100)
gs_results = runpf(case, opt)

# 比较计算结果
println("高斯-赛德尔法迭代次数: $(gs_results.iterations)")
println("高斯-赛德尔法计算的系统损耗: $(sum(gs_results.branch[:, PL])) MW")
```

## 直流潮流计算

### 基本直流潮流

直流潮流提供了快速的近似解：

```julia
# 执行直流潮流计算
dc_results = rundcpf(case)

# 查看计算结果
println("直流潮流计算的系统总发电: $(sum(dc_results.gen[:, PG])) MW")
println("直流潮流计算的系统总负荷: $(sum(dc_results.bus[:, PD])) MW")

# 比较交流和直流潮流的计算速度
using BenchmarkTools

println("交流潮流计算性能:")
@btime runpf($case);

println("直流潮流计算性能:")
@btime rundcpf($case);
```

### 直流潮流的线性敏感度分析

直流潮流可用于快速进行线性敏感度分析：

```julia
# 计算功率传输分布因子(PTDF)
ptdf_matrix = calc_ptdf(case)

# 查看特定线路对特定注入的敏感度
line_idx = 5
bus_idx = 3
println("线路 $(case.branch[line_idx, F_BUS])-$(case.branch[line_idx, T_BUS]) 对母线 $bus_idx 注入的敏感度: $(ptdf_matrix[line_idx, bus_idx])")

# 计算线路中断因子(LODF)
lodf_matrix = calc_lodf(case)

# 查看特定线路中断对其他线路的影响
outage_line = 1
affected_line = 5
println("线路 $(case.branch[outage_line, F_BUS])-$(case.branch[outage_line, T_BUS]) 中断对线路 $(case.branch[affected_line, F_BUS])-$(case.branch[affected_line, T_BUS]) 的影响系数: $(lodf_matrix[affected_line, outage_line])")
```

## 不平衡潮流计算

### 三相不平衡潮流

处理三相不平衡系统：

```julia
# 加载三相不平衡系统案例
unbalanced_case = load_case("case_unbalanced.m")

# 执行不平衡潮流计算
opt = Dict("verbose" => 1)
upf_results = runupf(unbalanced_case, opt)

# 查看各相电压
for i in 1:size(upf_results.bus, 1)
    bus_id = upf_results.bus[i, BUS_I]
    println("母线 $bus_id:")
    println("  A相电压: $(upf_results.bus[i, VM_A]) ∠ $(upf_results.bus[i, VA_A])°")
    println("  B相电压: $(upf_results.bus[i, VM_B]) ∠ $(upf_results.bus[i, VA_B])°")
    println("  C相电压: $(upf_results.bus[i, VM_C]) ∠ $(upf_results.bus[i, VA_C])°")
end

# 计算电压不平衡度
for i in 1:size(upf_results.bus, 1)
    vm_a = upf_results.bus[i, VM_A]
    vm_b = upf_results.bus[i, VM_B]
    vm_c = upf_results.bus[i, VM_C]
    
    vm_avg = (vm_a + vm_b + vm_c) / 3
    vm_dev = max(abs(vm_a - vm_avg), abs(vm_b - vm_avg), abs(vm_c - vm_avg))
    vuf = vm_dev / vm_avg * 100  # 电压不平衡度(%)
    
    println("母线 $(upf_results.bus[i, BUS_I]) 电压不平衡度: $vuf%")
end
```

## 高级选项与自定义

### 自定义收敛标准

设置自定义的收敛标准：

```julia
# 设置严格的收敛标准
strict_opt = Dict(
    "pf_tol" => 1e-8,        # 更严格的收敛容差
    "max_it" => 50,          # 增加最大迭代次数
    "verbose" => 2           # 更详细的输出
)

# 执行潮流计算
strict_results = runpf(case, strict_opt)

println("使用严格收敛标准:")
println("  迭代次数: $(strict_results.iterations)")
println("  最终误差: $(strict_results.max_mismatch)")
```

### 自定义起始点

设置自定义的计算起始点：

```julia
# 保存原始案例数据
original_case = deepcopy(case)

# 修改起始电压值
for i in 1:size(case.bus, 1)
    # 设置平坦起始点
    case.bus[i, VM] = 1.0
    case.bus[i, VA] = 0.0
end

# 使用平坦起始点计算潮流
flat_results = runpf(case)

# 使用"热启动"起始点
warm_start_case = deepcopy(original_case)
# 将上一次的结果作为起始点
for i in 1:size(warm_start_case.bus, 1)
    bus_idx = findfirst(x -> x == warm_start_case.bus[i, BUS_I], flat_results.bus[:, BUS_I])
    if bus_idx !== nothing
        warm_start_case.bus[i, VM] = flat_results.bus[bus_idx, VM]
        warm_start_case.bus[i, VA] = flat_results.bus[bus_idx, VA]
    end
end

# 使用热启动点计算潮流
warm_results = runpf(warm_start_case)

# 比较迭代次数
println("平坦起始点迭代次数: $(flat_results.iterations)")
println("热启动起始点迭代次数: $(warm_results.iterations)")
```

### 处理收敛问题

处理难以收敛的情况：

```julia
# 创建一个难以收敛的案例
difficult_case = deepcopy(case)
# 增加负荷，减少电压支撑
difficult_case.bus[:, PD] *= 1.5
# 降低发电机电压设定点
difficult_case.gen[:, VG] *= 0.95

# 尝试标准计算
try
    standard_results = runpf(difficult_case)
    println("标准方法成功收敛")
catch e
    println("标准方法未能收敛: $e")
    
    # 使用改进的选项
    robust_opt = Dict(
        "pf_alg" => "NR",
        "max_it" => 100,       # 增加迭代次数
        "pf_tol" => 1e-4,      # 放宽收敛标准
        "enforce_q_lims" => 0, # 暂时不考虑无功限制
        "verbose" => 2
    )
    
    try
        robust_results = runpf(difficult_case, robust_opt)
        println("改进方法成功收敛，迭代次数: $(robust_results.iterations)")
    catch e2
        println("改进方法仍未收敛: $e2")
        
        # 使用连续潮流法
        continuation_opt = Dict(
            "pf_alg" => "CONT",  # 连续潮流法
            "verbose" => 2
        )
        
        try
            cont_results = runpf(difficult_case, continuation_opt)
            println("连续潮流法成功收敛")
        catch e3
            println("连续潮流法未能收敛: $e3")
        end
    end
end
```

## GPU加速计算

### 使用GPU加速大规模系统计算

对于大规模系统，可以使用GPU加速计算：

```julia
# 检查是否有可用的GPU
using CUDA
if CUDA.functional()
    # 加载大规模系统
    large_case = load_case("case1888rte.m")
    
    # 使用CPU计算
    println("使用CPU计算...")
    @time cpu_results = runpf(large_case)
    
    # 使用GPU加速计算
    println("使用GPU加速计算...")
    gpu_opt = Dict("use_gpu" => true)
    @time gpu_results = runpf(large_case, gpu_opt)
    
    # 比较结果
    cpu_losses = sum(cpu_results.branch[:, PL])
    gpu_losses = sum(gpu_results.branch[:, PL])
    println("CPU计算的系统损耗: $cpu_losses MW")
    println("GPU计算的系统损耗: $gpu_losses MW")
    println("相对误差: $(abs(cpu_losses - gpu_losses) / cpu_losses * 100)%")
else
    println("没有可用的GPU或CUDA环境")
end
```

## 混合交直流系统

### 混合交直流潮流计算

计算包含交流和直流部分的混合系统：

```julia
# 加载混合交直流系统案例
hybrid_case = load_case("case_hybrid.m")

# 执行混合交直流潮流计算
hybrid_results = runpf(hybrid_case)

# 查看交流部分结果
println("交流系统结果:")
println("  交流母线数量: $(sum(hybrid_results.bus[:, BUS_TYPE] .!= BUS_TYPE_DC))")
println("  交流系统总发电: $(sum(hybrid_results.gen[hybrid_results.gen[:, GEN_STATUS] .== 1, PG])) MW")
println("  交流系统总负荷: $(sum(hybrid_results.bus[hybrid_results.bus[:, BUS_TYPE] .!= BUS_TYPE_DC, PD])) MW")
println("  交流系统总损耗: $(sum(hybrid_results.branch[hybrid_results.branch[:, BR_TYPE] .== BR_TYPE_AC, PL])) MW")

# 查看直流部分结果
println("直流系统结果:")
println("  直流母线数量: $(sum(hybrid_results.bus[:, BUS_TYPE] .== BUS_TYPE_DC))")
println("  直流系统总负荷: $(sum(hybrid_results.bus[hybrid_results.bus[:, BUS_TYPE] .== BUS_TYPE_DC, PD])) MW")
println("  直流系统总损耗: $(sum(hybrid_results.branch[hybrid_results.branch[:, BR_TYPE] .== BR_TYPE_DC, PL])) MW")

# 查看换流器结果
println("换流器结果:")
for i in 1:size(hybrid_results.converter, 1)
    conv_id = hybrid_results.converter[i, CONV_ID]
    ac_bus = hybrid_results.converter[i, CONV_AC_BUS]
    dc_bus = hybrid_results.converter[i, CONV_DC_BUS]
    p_ac = hybrid_results.converter[i, CONV_P_AC]
    p_dc = hybrid_results.converter[i, CONV_P_DC]
    loss = hybrid_results.converter[i, CONV_LOSS]
    
    println("  换流器 $conv_id (AC母线 $ac_bus - DC母线 $dc_bus):")
    println("    交流侧功率: $p_ac MW")
    println("    直流侧功率: $p_dc MW")
    println("    换流损耗: $loss MW")
    println("    效率: $((abs(p_dc) / abs(p_ac)) * 100)%")
end
```

## 结果分析与处理

### 详细的结果分析

深入分析潮流计算结果：

```julia
# 执行潮流计算
results = runpf(case)

# 分析母线电压分布
vm = results.bus[:, VM]
va = results.bus[:, VA]

println("电压统计:")
println("  最高电压: $(maximum(vm)) pu，位于母线 $(results.bus[argmax(vm), BUS_I])")
println("  最低电压: $(minimum(vm)) pu，位于母线 $(results.bus[argmin(vm), BUS_I])")
println("  平均电压: $(mean(vm)) pu")
println("  标准差: $(std(vm)) pu")

# 分析支路功率流
branch_flow = results.branch[:, PF]
branch_from = results.branch[:, F_BUS]
branch_to = results.branch[:, T_BUS]
branch_rate = results.branch[:, RATE_A]

# 计算支路负载率
loading_rate = abs.(branch_flow) ./ branch_rate * 100
branch_rate_valid = branch_rate .> 0  # 只考虑有额定值的支路

if any(branch_rate_valid)
    max_loading_idx = argmax(loading_rate .* branch_rate_valid)
    println("支路负载率统计(对有额定值的支路):")
    println("  最高负载率: $(loading_rate[max_loading_idx])%，位于支路 $(branch_from[max_loading_idx])-$(branch_to[max_loading_idx])")
    println("  平均负载率: $(mean(loading_rate[branch_rate_valid]))%")
end

# 分析发电机出力
gen_pg = results.gen[:, PG]
gen_qg = results.gen[:, QG]
gen_bus = results.gen[:, GEN_BUS]
gen_pmax = results.gen[:, PMAX]
gen_pmin = results.gen[:, PMIN]
gen_qmax = results.gen[:, QMAX]
gen_qmin = results.gen[:, QMIN]

# 检查发电机是否达到限制
p_at_max = abs.(gen_pg - gen_pmax) .< 1e-3
p_at_min = abs.(gen_pg - gen_pmin) .< 1e-3
q_at_max = abs.(gen_qg - gen_qmax) .< 1e-3
q_at_min = abs.(gen_qg - gen_qmin) .< 1e-3

println("发电机限制状态:")
println("  有功达到上限的发电机数量: $(sum(p_at_max))")
println("  有功达到下限的发电机数量: $(sum(p_at_min))")
println("  无功达到上限的发电机数量: $(sum(q_at_max))")
println("  无功达到下限的发电机数量: $(sum(q_at_min))")

# 分析系统损耗
total_gen = sum(gen_pg)
total_load = sum(results.bus[:, PD])
total_loss = total_gen - total_load

println("系统损耗分析:")
println("  总发电: $total_gen MW")
println("  总负荷: $total_load MW")
println("  总损耗: $total_loss MW")
println("  损耗率: $(total_loss / total_gen * 100)%")
```

### 结果导出与可视化

导出和可视化计算结果：

```julia
using DataFrames, CSV, Plots

# 将结果转换为DataFrame格式
# 母线结果
bus_df = DataFrame(
    Bus_ID = results.bus[:, BUS_I],
    Type = results.bus[:, BUS_TYPE],
    Voltage_Magnitude = results.bus[:, VM],
    Voltage_Angle = results.bus[:, VA],
    Load_MW = results.bus[:, PD],
    Load_MVAr = results.bus[:, QD]
)

# 支路结果
branch_df = DataFrame(
    From_Bus = results.branch[:, F_BUS],
    To_Bus = results.branch[:, T_BUS],
    Power_Flow_MW = results.branch[:, PF],
    Power_Flow_MVAr = results.branch[:, QF],
    Loss_MW = results.branch[:, PL],
    Loss_MVAr = results.branch[:, QL]
)

# 发电机结果
gen_df = DataFrame(
    Bus_ID = results.gen[:, GEN_BUS],
    Active_Power_MW = results.gen[:, PG],
    Reactive_Power_MVAr = results.gen[:, QG],
    Voltage_Setpoint = results.gen[:, VG]
)

# 导出到CSV文件
CSV.write("bus_results.csv", bus_df)
CSV.write("branch_results.csv", branch_df)
CSV.write("gen_results.csv", gen_df)

# 创建电压分布图
p1 = bar(bus_df.Bus_ID, bus_df.Voltage_Magnitude,
    title="母线电压分布",
    xlabel="母线编号",
    ylabel="电压幅值(pu)",
    legend=false)

# 创建功率流分布图
p2 = bar(1:size(branch_df, 1), abs.(branch_df.Power_Flow_MW),
    title="支路功率流分布",
    xlabel="支路编号",
    ylabel="功率流大小(MW)",
    legend=false)

# 创建损耗分布图
p3 = bar(1:size(branch_df, 1), branch_df.Loss_MW,
    title="支路损耗分布",
    xlabel="支路编号",
    ylabel="有功损耗(MW)",
    legend=false)

# 组合图表
plot(p1, p2, p3, layout=(3,1), size=(800, 1200))
savefig("power_flow_results.png")
```

## 高级应用

### 或然性分析

进行简单的或然性分析：

```julia
# 定义要分析的线路中断情况
outage_branches = [1, 5, 10]

# 保存原始案例
base_case = deepcopy(case)
base_results = runpf(base_case)

println("基础情况:")
println("  总发电: $(sum(base_results.gen[:, PG])) MW")
println("  总损耗: $(sum(base_results.branch[:, PL])) MW")
println("  最低电压: $(minimum(base_results.bus[:, VM])) pu")

# 分析每个中断情况
for br_idx in outage_branches
    # 创建中断情况
    outage_case = deepcopy(base_case)
    from_bus = outage_case.branch[br_idx, F_BUS]
    to_bus = outage_case.branch[br_idx, T_BUS]
    outage_case.branch[br_idx, BR_STATUS] = 0  # 断开线路
    
    println("\n分析线路 $from_bus-$to_bus 中断情况:")
    
    # 执行潮流计算
    try
        outage_results = runpf(outage_case)
        
        # 分析结果
        total_gen = sum(outage_results.gen[:, PG])
        total_loss = sum(outage_results.branch[:, PL])
        min_voltage = minimum(outage_results.bus[:, VM])
        
        # 找出最大负载率
        loading_rate = abs.(outage_results.branch[:, PF]) ./ outage_results.branch[:, RATE_A] * 100
        branch_rate_valid = outage_results.branch[:, RATE_A] .> 0
        if any(branch_rate_valid)
            max_loading = maximum(loading_rate .* branch_rate_valid)
            max_loading_idx = argmax(loading_rate .* branch_rate_valid)
            overloaded = loading_rate .> 100 .& branch_rate_valid
            
            println("  总发电: $total_gen MW")
            println("  总损耗: $total_loss MW")
            println("  最低电压: $min_voltage pu")
            println("  最大支路负载率: $max_loading%")
            println("  过载支路数量: $(sum(overloaded))")
            
            if sum(overloaded) > 0
                for i in findall(overloaded)
                    from = outage_results.branch[i, F_BUS]
                    to = outage_results.branch[i, T_BUS]
                    load_rate = loading_rate[i]
                    println("    过载支路 $from-$to: $load_rate%")
                end
            end
        end
    catch e
        println("  潮流计算未收敛: $e")
    end
end
```

### 灵敏度分析

进行简单的灵敏度分析：

```julia
# 执行基础潮流计算
base_results = runpf(case)
base_losses = sum(base_results.branch[:, PL])

# 定义要分析的参数和变化范围
load_factors = 0.8:0.05:1.2

# 存储结果
loss_results = zeros(length(load_factors))
voltage_results = zeros(length(load_factors))

# 进行灵敏度分析
for (i, factor) in enumerate(load_factors)
    # 修改负荷
    test_case = deepcopy(case)
    test_case.bus[:, PD] = case.bus[:, PD] * factor
    test_case.bus[:, QD] = case.bus[:, QD] * factor
    
    # 执行潮流计算
    try
        test_results = runpf(test_case)
        
        # 记录结果
        loss_results[i] = sum(test_results.branch[:, PL])
        voltage_results[i] = minimum(test_results.bus[:, VM])
    catch e
        println("负荷因子 $factor 时潮流计算未收敛: $e")
        loss_results[i] = NaN
        voltage_results[i] = NaN
    end
end

# 计算灵敏度
valid_indices = .!isnan.(loss_results)
if sum(valid_indices) >= 2
    # 计算损耗对负荷的灵敏度
    loss_sensitivity = (loss_results[valid_indices] .- base_losses) ./ (load_factors[valid_indices] .- 1.0) ./ base_losses
    
    println("损耗对负荷变化的灵敏度:")
    for (f, s) in zip(load_factors[valid_indices], loss_sensitivity)
        println("  负荷因子 $f: 灵敏度 $s")
    end
    
    # 绘制灵敏度曲线
    p1 = plot(load_factors, loss_results,
        title="系统损耗对负荷变化的灵敏度",
        xlabel="负荷因子",
        ylabel="系统损耗(MW)",
        marker=:circle,
        label="系统损耗")
    
    p2 = plot(load_factors, voltage_results,
        title="最低电压对负荷变化的灵敏度",
        xlabel="负荷因子",
        ylabel="最低电压(pu)",
        marker=:circle,
        label="最低电压")
    
    plot(p1, p2, layout=(2,1), size=(800, 800))
    savefig("sensitivity_analysis.png")
end
```

## 结论

本文档展示了 PowerFlow 模块的各种功能和应用场景，从基础的潮流计算到高级的系统分析。通过这些示例，用户可以了解如何:

1. 使用不同的潮流计算算法
2. 配置高级计算选项
3. 处理混合交直流系统
4. 分析和可视化计算结果
5. 进行或然性和灵敏度分析

PowerFlow 模块提供了强大而灵活的功能，可以满足从教学到研究再到工程实践的各种需求。通过合理配置和使用这些功能，用户可以高效地分析和解决电力系统中的各种问题。