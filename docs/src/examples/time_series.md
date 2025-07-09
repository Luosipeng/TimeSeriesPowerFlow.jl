# TimeSeriesPowerFlow 模块示例

本文档提供了 TimeSeriesPowerFlow 模块的详细使用示例，展示了时间序列潮流计算、动态调度、可再生能源集成以及结果可视化的应用方法。这些示例旨在帮助用户深入理解 TimeSeriesPowerFlow 模块的核心功能。

## 基础时间序列潮流计算

### 加载数据与基本计算

执行基本的时间序列潮流计算：

```julia
using PowerFlow
using TimeSeriesPowerFlow
using Plots

# 加载标准测试系统
case = load_case("case14.m")

# 准备24小时负荷数据
# 这里假设load_data.xlsx包含了24小时的负荷曲线数据
load_file = "load_data.xlsx"

# 执行时间序列潮流计算
results = runtdpf(case, load_file)

# 检查计算结果
println("时间序列潮流计算结果:")
println("  时间点数量: $(length(results.time_points))")
println("  成功计算的时间点数量: $(sum(results.success))")

# 查看第一个时间点的结果
if results.success[1]
    first_result = results.bus_results[1]
    println("第一个时间点:")
    println("  总负荷: $(sum(first_result.bus[:, PD])) MW")
    println("  总发电: $(sum(first_result.gen[:, PG])) MW")
    println("  总损耗: $(sum(first_result.branch[:, PL])) MW")
end
```

### 结果可视化

可视化时间序列潮流计算结果：

```julia
# 提取时间序列数据
time_points = results.time_points
total_load = zeros(length(time_points))
total_gen = zeros(length(time_points))
total_loss = zeros(length(time_points))
min_voltage = zeros(length(time_points))

for t in 1:length(time_points)
    if results.success[t]
        result = results.bus_results[t]
        total_load[t] = sum(result.bus[:, PD])
        total_gen[t] = sum(result.gen[:, PG])
        total_loss[t] = sum(result.branch[:, PL])
        min_voltage[t] = minimum(result.bus[:, VM])
    else
        # 对未收敛的时间点使用NaN
        total_load[t] = NaN
        total_gen[t] = NaN
        total_loss[t] = NaN
        min_voltage[t] = NaN
    end
end

# 绘制负荷曲线
p1 = plot(time_points, total_load,
    title="24小时负荷曲线",
    xlabel="时间(小时)",
    ylabel="总负荷(MW)",
    marker=:circle,
    label="总负荷",
    legend=:topleft)

# 绘制发电曲线
p2 = plot(time_points, total_gen,
    title="24小时发电曲线",
    xlabel="时间(小时)",
    ylabel="总发电(MW)",
    marker=:circle,
    label="总发电",
    legend=:topleft)

# 绘制损耗曲线
p3 = plot(time_points, total_loss,
    title="24小时系统损耗曲线",
    xlabel="时间(小时)",
    ylabel="总损耗(MW)",
    marker=:circle,
    label="总损耗",
    legend=:topleft)

# 绘制最低电压曲线
p4 = plot(time_points, min_voltage,
    title="24小时最低电压曲线",
    xlabel="时间(小时)",
    ylabel="最低电压(pu)",
    marker=:circle,
    label="最低电压",
    legend=:bottomleft)

# 组合图表
plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
savefig("time_series_results.png")
```

### 使用内置可视化函数

使用模块提供的内置可视化函数：

```julia
# 使用内置函数绘制电压时间序列
plot_time_day = 1  # 分析第一天的数据
plot_voltage_time_series(results, case, plot_time_day)

# 绘制系统损耗时间序列
plot_losses_time_series(results, case, plot_time_day)

# 绘制特定母线的电压时间序列
bus_id = 5  # 选择特定母线
plot_bus_voltage_time_series(results, case, bus_id, plot_time_day)

# 绘制特定支路的功率流时间序列
branch_idx = 3  # 选择特定支路
plot_branch_flow_time_series(results, case, branch_idx, plot_time_day)
```

## 包含可再生能源的时间序列分析

### 光伏发电集成

模拟包含光伏发电的时间序列潮流：

```julia
# 准备光照强度数据
# 假设irradiance_data.xlsx包含了24小时的光照强度数据
irradiance_file = "irradiance_data.xlsx"

# 执行包含光照数据的时间序列潮流计算
pv_results = runtdpf(case, load_file, irradiance_file)

# 分析光伏发电输出
# 假设系统中的PV发电机编号为[3, 5]
pv_gen_indices = [3, 5]

# 提取光伏发电时间序列
pv_generation = zeros(length(pv_results.time_points), length(pv_gen_indices))

for t in 1:length(pv_results.time_points)
    if pv_results.success[t]
        for (i, gen_idx) in enumerate(pv_gen_indices)
            pv_generation[t, i] = pv_results.bus_results[t].gen[gen_idx, PG]
        end
    end
end

# 绘制光伏发电曲线
plot(pv_results.time_points, pv_generation,
    title="光伏发电时间序列",
    xlabel="时间(小时)",
    ylabel="发电功率(MW)",
    label=["PV发电机 $(pv_gen_indices[i])" for i in 1:length(pv_gen_indices)],
    marker=:circle,
    legend=:topleft)
savefig("pv_generation.png")

# 分析光伏渗透率对系统的影响
total_gen_with_pv = zeros(length(pv_results.time_points))
total_loss_with_pv = zeros(length(pv_results.time_points))

for t in 1:length(pv_results.time_points)
    if pv_results.success[t]
        total_gen_with_pv[t] = sum(pv_results.bus_results[t].gen[:, PG])
        total_loss_with_pv[t] = sum(pv_results.bus_results[t].branch[:, PL])
    end
end

# 计算光伏发电占比
pv_penetration = sum(pv_generation, dims=2) ./ total_gen_with_pv * 100

# 绘制光伏渗透率曲线
plot(pv_results.time_points, pv_penetration,
    title="光伏渗透率时间序列",
    xlabel="时间(小时)",
    ylabel="光伏渗透率(%)",
    marker=:circle,
    legend=false)
savefig("pv_penetration.png")
```

### 最大功率点跟踪

模拟光伏发电的最大功率点跟踪：

```julia
# 设置MPPT选项
mppt_opt = Dict(
    "pv_mode" => "MPPT",
    "verbose" => 1
)

# 执行带有MPPT的时间序列潮流计算
mppt_results = runtdpf(case, load_file, irradiance_file, nothing, mppt_opt)

# 比较MPPT模式与默认模式的光伏发电
mppt_generation = zeros(length(mppt_results.time_points), length(pv_gen_indices))

for t in 1:length(mppt_results.time_points)
    if mppt_results.success[t]
        for (i, gen_idx) in enumerate(pv_gen_indices)
            mppt_generation[t, i] = mppt_results.bus_results[t].gen[gen_idx, PG]
        end
    end
end

# 绘制比较图
for i in 1:length(pv_gen_indices)
    plot(pv_results.time_points, [pv_generation[:, i] mppt_generation[:, i]],
        title="PV发电机 $(pv_gen_indices[i]) 发电比较",
        xlabel="时间(小时)",
        ylabel="发电功率(MW)",
        label=["默认模式" "MPPT模式"],
        marker=:circle,
        legend=:topleft)
    savefig("pv_mppt_comparison_$(pv_gen_indices[i]).png")
end
```

## 动态经济调度

### 基本经济调度

执行基本的动态经济调度：

```julia
# 准备电价数据
# 假设price_data.xlsx包含了24小时的电价数据
price_file = "price_data.xlsx"

# 设置动态调度选项
dispatch_opt = Dict(
    "mode" => "dynamic_dispatch",
    "verbose" => 1
)

# 执行动态经济调度计算
dispatch_results = runtdpf(case, load_file, irradiance_file, price_file, dispatch_opt)

# 分析调度结果
# 提取各发电机的出力曲线
gen_dispatch = zeros(length(dispatch_results.time_points), size(case.gen, 1))
prices = zeros(length(dispatch_results.time_points))

for t in 1:length(dispatch_results.time_points)
    if dispatch_results.success[t]
        for g in 1:size(case.gen, 1)
            gen_dispatch[t, g] = dispatch_results.bus_results[t].gen[g, PG]
        end
        prices[t] = dispatch_results.prices[t]
    end
end

# 绘制发电机出力曲线
plot(dispatch_results.time_points, gen_dispatch,
    title="发电机经济调度结果",
    xlabel="时间(小时)",
    ylabel="发电功率(MW)",
    label=["发电机 $g" for g in 1:size(case.gen, 1)],
    legend=:outerright)
savefig("economic_dispatch.png")

# 绘制电价曲线
plot(dispatch_results.time_points, prices,
    title="24小时电价曲线",
    xlabel="时间(小时)",
    ylabel="电价($/MWh)",
    marker=:circle,
    legend=false)
savefig("price_curve.png")

# 计算总发电成本
total_cost = sum(prices .* sum(gen_dispatch, dims=2))
println("总发电成本: $total_cost $")
```

### 考虑网络约束的经济调度

执行考虑网络约束的经济调度：

```julia
# 设置考虑网络约束的动态调度选项
network_dispatch_opt = Dict(
    "mode" => "network_constrained_dispatch",
    "enforce_branch_limits" => 1,
    "verbose" => 1
)

# 执行考虑网络约束的动态经济调度
network_results = runtdpf(case, load_file, irradiance_file, price_file, network_dispatch_opt)

# 分析网络约束下的调度结果
network_gen_dispatch = zeros(length(network_results.time_points), size(case.gen, 1))

for t in 1:length(network_results.time_points)
    if network_results.success[t]
        for g in 1:size(case.gen, 1)
            network_gen_dispatch[t, g] = network_results.bus_results[t].gen[g, PG]
        end
    end
end

# 比较有无网络约束的调度结果
for g in 1:size(case.gen, 1)
    plot(dispatch_results.time_points, [gen_dispatch[:, g] network_gen_dispatch[:, g]],
        title="发电机 $g 调度比较",
        xlabel="时间(小时)",
        ylabel="发电功率(MW)",
        label=["无网络约束" "有网络约束"],
        marker=:circle,
        legend=:topleft)
    savefig("dispatch_comparison_gen_$g.png")
end

# 分析网络约束的影响
# 检查是否有支路达到限制
branch_loading = zeros(length(network_results.time_points), size(case.branch, 1))
branch_limits = case.branch[:, RATE_A]

for t in 1:length(network_results.time_points)
    if network_results.success[t]
        for br in 1:size(case.branch, 1)
            if branch_limits[br] > 0  # 只考虑有限制的支路
                branch_loading[t, br] = abs(network_results.bus_results[t].branch[br, PF]) / branch_limits[br] * 100
            end
        end
    end
end

# 找出负载率最高的支路
max_loading = maximum(branch_loading, dims=1)[:]
critical_branches = findall(max_loading .> 90)  # 负载率超过90%的支路

if !isempty(critical_branches)
    println("关键支路(最大负载率超过90%):")
    for br in critical_branches
        println("  支路 $(case.branch[br, F_BUS])-$(case.branch[br, T_BUS]): $(max_loading[br])%")
    end
    
    # 绘制关键支路的负载率曲线
    plot(network_results.time_points, branch_loading[:, critical_branches],
        title="关键支路负载率",
        xlabel="时间(小时)",
        ylabel="负载率(%)",
        label=["支路 $(case.branch[br, F_BUS])-$(case.branch[br, T_BUS])" for br in critical_branches],
        legend=:outerright)
    savefig("critical_branch_loading.png")
end
```

## 线性分布潮流

### 基本线性分布潮流计算

执行线性分布潮流计算：

```julia
# 设置线性分布潮流选项
lindist_opt = Dict(
    "pf_alg" => "LINDISTFLOW",
    "verbose" => 1
)

# 执行线性分布潮流计算
lindist_results = runlindistflow(case, load_file)

# 比较线性分布潮流与完整潮流计算的结果
# 提取两种方法的电压结果
voltage_full = zeros(length(results.time_points), size(case.bus, 1))
voltage_lindist = zeros(length(lindist_results.time_points), size(case.bus, 1))

for t in 1:length(results.time_points)
    if results.success[t] && lindist_results.success[t]
        for b in 1:size(case.bus, 1)
            voltage_full[t, b] = results.bus_results[t].bus[b, VM]
            voltage_lindist[t, b] = lindist_results.bus_results[t].bus[b, VM]
        end
    end
end

# 计算电压误差
voltage_error = abs.(voltage_full - voltage_lindist)
max_error = maximum(voltage_error)
avg_error = mean(voltage_error)

println("线性分布潮流与完整潮流计算的电压比较:")
println("  最大误差: $max_error pu")
println("  平均误差: $avg_error pu")

# 绘制特定母线的电压比较
bus_to_plot = 5  # 选择一个母线进行比较
plot(results.time_points, [voltage_full[:, bus_to_plot] voltage_lindist[:, bus_to_plot]],
    title="母线 $(case.bus[bus_to_plot, BUS_I]) 电压比较",
    xlabel="时间(小时)",
    ylabel="电压幅值(pu)",
    label=["完整潮流" "线性分布潮流"],
    marker=:circle,
    legend=:bottomleft)
savefig("voltage_comparison_bus_$(case.bus[bus_to_plot, BUS_I]).png")

# 比较计算时间
using BenchmarkTools

println("完整潮流计算性能:")
@btime runtdpf($case, $load_file);

println("线性分布潮流计算性能:")
@btime runlindistflow($case, $load_file);
```

### 配电网应用

在配电网中应用线性分布潮流：

```julia
# 加载配电网测试系统
dist_case = load_case("case_distribution.m")

# 执行线性分布潮流计算
dist_lindist_results = runlindistflow(dist_case, load_file)

# 分析配电网中的电压分布
# 提取末端节点的电压
# 假设末端节点是编号最大的几个节点
num_end_nodes = 5
num_buses = size(dist_case.bus, 1)
end_nodes = (num_buses-num_end_nodes+1):num_buses

# 提取末端节点的电压时间序列
end_voltage = zeros(length(dist_lindist_results.time_points), length(end_nodes))

for t in 1:length(dist_lindist_results.time_points)
    if dist_lindist_results.success[t]
        for (i, node) in enumerate(end_nodes)
            end_voltage[t, i] = dist_lindist_results.bus_results[t].bus[node, VM]
        end
    end
end

# 绘制末端节点电压曲线
plot(dist_lindist_results.time_points, end_voltage,
    title="配电网末端节点电压",
    xlabel="时间(小时)",
    ylabel="电压幅值(pu)",
    label=["节点 $(dist_case.bus[node, BUS_I])" for node in end_nodes],
    legend=:bottomleft)
savefig("distribution_end_voltage.png")

# 检查电压越限情况
voltage_violations = (end_voltage .< 0.95) .| (end_voltage .> 1.05)
has_violations = any(voltage_violations)

if has_violations
    println("检测到电压越限情况:")
    for (i, node) in enumerate(end_nodes)
        node_violations = sum(voltage_violations[:, i])
        if node_violations > 0
            println("  节点 $(dist_case.bus[node, BUS_I]) 有 $node_violations 个时间点电压越限")
        end
    end
end
```

## 违例分析与可视化

### 电压违例分析

分析时间序列中的电压违例：

```julia
# 设置电压限制
v_min = 0.95
v_max = 1.05

# 分析电压违例
voltage_violations = analyze_voltage_violations(results, case, v_min, v_max)

# 显示违例摘要
println("电压违例摘要:")
println("  低电压违例数量: $(voltage_violations.low_count)")
println("  高电压违例数量: $(voltage_violations.high_count)")

# 如果有违例，显示详细信息
if voltage_violations.low_count > 0 || voltage_violations.high_count > 0
    println("违例详情:")
    
    # 低电压违例
    if voltage_violations.low_count > 0
        println("  低电压违例:")
        for violation in voltage_violations.low_violations
            println("    时间点 $(violation.time), 母线 $(violation.bus_id): $(violation.value) pu")
        end
    end
    
    # 高电压违例
    if voltage_violations.high_count > 0
        println("  高电压违例:")
        for violation in voltage_violations.high_violations
            println("    时间点 $(violation.time), 母线 $(violation.bus_id): $(violation.value) pu")
        end
    end
    
    # 绘制违例图表
    plot_voltage_violations(results, case, v_min, v_max)
end
```

### 支路功率流违例分析

分析时间序列中的支路功率流违例：

```julia
# 设置违例阈值(%)
flow_threshold = 90

# 分析支路功率流违例
flow_violations = analyze_flow_violations(results, case, flow_threshold)

# 显示违例摘要
println("支路功率流违例摘要:")
println("  违例数量: $(flow_violations.count)")

# 如果有违例，显示详细信息
if flow_violations.count > 0
    println("违例详情:")
    for violation in flow_violations.violations
        from_bus = case.branch[violation.branch_idx, F_BUS]
        to_bus = case.branch[violation.branch_idx, T_BUS]
        println("  时间点 $(violation.time), 支路 $from_bus-$to_bus: $(violation.value)% of rating")
    end
    
    # 绘制违例图表
    plot_flow_violations(results, case, flow_threshold)
end
```

### 时间序列动画

创建时间序列结果的动画：

```julia
using Plots

# 创建电压分布动画
anim = @animate for t in 1:length(results.time_points)
    if results.success[t]
        result = results.bus_results[t]
        bar(result.bus[:, BUS_I], result.bus[:, VM],
            title="电压分布 - 时间点 $(results.time_points[t])",
            xlabel="母线编号",
            ylabel="电压幅值(pu)",
            ylims=(0.9, 1.1),
            legend=false)
    end
end

gif(anim, "voltage_animation.gif", fps=2)

# 创建系统负荷和损耗动画
anim2 = @animate for t in 1:length(results.time_points)
    if results.success[t]
        result = results.bus_results[t]
        total_load = sum(result.bus[:, PD])
        total_loss = sum(result.branch[:, PL])
        
        bar([1, 2], [total_load, total_loss],
            title="系统负荷和损耗 - 时间点 $(results.time_points[t])",
            xlabel="",
            ylabel="功率(MW)",
            xticks=(1:2, ["总负荷", "总损耗"]),
            legend=false)
    end
end

gif(anim2, "load_loss_animation.gif", fps=2)
```

## 高级应用

### 储能系统集成

模拟包含储能系统的时间序列潮流：

```julia
# 加载包含储能系统的案例
storage_case = load_case("case_with_storage.m")

# 设置储能调度选项
storage_opt = Dict(
    "storage_control" => "optimal",
    "verbose" => 1
)

# 执行包含储能的时间序列潮流计算
storage_results = runtdpf(storage_case, load_file, irradiance_file, price_file, storage_opt)

# 分析储能系统的充放电行为
# 假设储能系统连接在特定母线上
storage_bus_id = 10
storage_p = zeros(length(storage_results.time_points))
storage_soc = zeros(length(storage_results.time_points))

for t in 1:length(storage_results.time_points)
    if storage_results.success[t]
        # 提取储能功率(正值表示放电，负值表示充电)
        storage_p[t] = storage_results.storage_results[t][storage_bus_id, STORAGE_P]
        # 提取储能SOC
        storage_soc[t] = storage_results.storage_results[t][storage_bus_id, STORAGE_SOC]
    end
end

# 绘制储能充放电曲线
p1 = plot(storage_results.time_points, storage_p,
    title="储能充放电功率",
    xlabel="时间(小时)",
    ylabel="功率(MW)",
    marker=:circle,
    label="储能功率(正值为放电)")

# 绘制储能SOC曲线
p2 = plot(storage_results.time_points, storage_soc * 100,
    title="储能SOC",
    xlabel="时间(小时)",
    ylabel="SOC(%)",
    marker=:circle,
    label="储能SOC")

# 组合图表
plot(p1, p2, layout=(2,1), size=(800, 600))
savefig("storage_behavior.png")

# 分析储能对系统的影响
# 比较有无储能的系统损耗
storage_loss = zeros(length(storage_results.time_points))
for t in 1:length(storage_results.time_points)
    if storage_results.success[t]
        storage_loss[t] = sum(storage_results.bus_results[t].branch[:, PL])
    end
end

# 绘制损耗比较曲线
plot(results.time_points, [total_loss storage_loss],
    title="系统损耗比较",
    xlabel="时间(小时)",
    ylabel="系统损耗(MW)",
    label=["无储能" "有储能"],
    marker=:circle,
    legend=:topleft)
savefig("loss_comparison_with_storage.png")
```

### 需求响应模拟

模拟需求响应对系统的影响：

```julia
# 设置需求响应选项
dr_opt = Dict(
    "demand_response" => true,
    "dr_price_threshold" => 50,  # 当电价超过此值时触发需求响应
    "dr_reduction" => 0.1,       # 需求响应减少10%的负荷
    "verbose" => 1
)

# 执行包含需求响应的时间序列潮流计算
dr_results = runtdpf(case, load_file, irradiance_file, price_file, dr_opt)

# 分析需求响应的效果
dr_load = zeros(length(dr_results.time_points))
for t in 1:length(dr_results.time_points)
    if dr_results.success[t]
        dr_load[t] = sum(dr_results.bus_results[t].bus[:, PD])
    end
end

# 绘制需求响应效果
plot(results.time_points, [total_load dr_load],
    title="需求响应效果",
    xlabel="时间(小时)",
    ylabel="总负荷(MW)",
    label=["基准负荷" "需求响应后负荷"],
    marker=:circle,
    legend=:topleft)
savefig("demand_response_effect.png")

# 计算需求响应的经济效益
if !isempty(dr_results.prices)
    base_cost = sum(total_load .* dr_results.prices)
    dr_cost = sum(dr_load .* dr_results.prices)
    savings = base_cost - dr_cost
    
    println("需求响应经济效益分析:")
    println("  基准情况总成本: $base_cost $")
    println("  需求响应后总成本: $dr_cost $")
    println("  节约成本: $savings $ ($(savings/base_cost*100)%)")
end
```

