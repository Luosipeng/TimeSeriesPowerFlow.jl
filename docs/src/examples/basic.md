# Distribution PowerFlow 库基础示例

本文档提供了 PowerFlow 库的基础使用示例，涵盖了从简单的潮流计算到时间序列分析的典型应用场景。这些示例旨在帮助新用户快速上手并理解库的核心功能。

## 安装与导入

首先，我们需要安装并导入 PowerFlow 库及其相关模块：

```julia
# 安装库（仅需执行一次）
using Pkg
Pkg.add("PowerFlow")

# 导入所需模块
using PowerFlow
using TimeSeriesPowerFlow
using Utils
using ComponentModel
```

## 基础潮流计算

### 加载案例数据

PowerFlow 支持多种格式的电力系统案例数据，包括 MATPOWER 格式的 `.m` 文件：

```julia
# 加载标准 IEEE 14 节点测试系统
case = load_case("case14.m")

# 查看案例基本信息
println("系统包含 $(length(case.bus)) 个母线和 $(length(case.branch)) 条支路")
```

### 执行交流潮流计算

最基本的功能是执行交流潮流计算：

```julia
# 使用默认参数执行交流潮流计算
results = runpf(case)

# 检查计算是否收敛
if results.success
    println("潮流计算成功收敛！")
else
    println("潮流计算未收敛，请检查系统参数")
end

# 查看计算结果摘要
println("系统总负荷: $(sum(results.bus[:, PD])) MW")
println("系统总损耗: $(sum(results.branch[:, PL])) MW")
```

### 自定义计算选项

可以通过选项字典自定义潮流计算的参数：

```julia
# 设置自定义计算选项
opt = Dict(
    "pf_alg" => "NR",        # 使用牛顿-拉夫森法
    "pf_tol" => 1e-6,        # 收敛容差
    "max_it" => 30,          # 最大迭代次数
    "verbose" => 1           # 输出详细信息
)

# 使用自定义选项执行潮流计算
results = runpf(case, opt)
```

### 执行直流潮流计算

对于某些应用，直流潮流计算提供了更快的近似解：

```julia
# 执行直流潮流计算
dc_results = rundcpf(case)

# 比较直流和交流潮流计算结果
println("交流潮流计算的系统总损耗: $(sum(results.branch[:, PL])) MW")
println("直流潮流计算的系统总损耗: $(sum(dc_results.branch[:, PL])) MW")
```

## 结果分析与可视化

### 电压分析

分析系统中的电压分布：

```julia
# 提取电压幅值和相角
vm = results.bus[:, VM]
va = results.bus[:, VA]

# 查找最高和最低电压
max_vm_idx = argmax(vm)
min_vm_idx = argmin(vm)

println("最高电压: $(vm[max_vm_idx]) pu，位于母线 $(results.bus[max_vm_idx, BUS_I])")
println("最低电压: $(vm[min_vm_idx]) pu，位于母线 $(results.bus[min_vm_idx, BUS_I])")

# 绘制电压分布图
using Plots
bar(results.bus[:, BUS_I], vm, 
    title="母线电压分布", 
    xlabel="母线编号", 
    ylabel="电压幅值(pu)",
    legend=false)
```

### 支路功率流分析

分析系统中的功率流分布：

```julia
# 提取支路有功功率流
branch_flow = results.branch[:, PF]
branch_from = results.branch[:, F_BUS]
branch_to = results.branch[:, T_BUS]

# 找出功率流最大的支路
max_flow_idx = argmax(abs.(branch_flow))
println("最大功率流: $(branch_flow[max_flow_idx]) MW，位于支路 $(branch_from[max_flow_idx])-$(branch_to[max_flow_idx])")

# 绘制功率流分布图
bar(1:length(branch_flow), abs.(branch_flow),
    title="支路功率流分布",
    xlabel="支路编号",
    ylabel="功率流大小(MW)",
    legend=false)
```

## 系统修改与重新计算

### 修改负荷

模拟负荷变化并重新计算潮流：

```julia
# 保存原始负荷数据
original_load = copy(case.bus[:, PD])

# 增加所有负荷20%
case.bus[:, PD] *= 1.2

# 重新计算潮流
new_results = runpf(case)

# 比较结果
println("原始总负荷: $(sum(original_load)) MW")
println("修改后总负荷: $(sum(case.bus[:, PD])) MW")
println("原始总损耗: $(sum(results.branch[:, PL])) MW")
println("修改后总损耗: $(sum(new_results.branch[:, PL])) MW")

# 恢复原始负荷
case.bus[:, PD] = original_load
```

### 修改发电机出力

调整发电机出力并观察系统响应：

```julia
# 保存原始发电机数据
original_gen = copy(case.gen)

# 修改特定发电机的出力
gen_idx = 1  # 第一台发电机
original_pg = case.gen[gen_idx, PG]
case.gen[gen_idx, PG] *= 0.8  # 降低到原来的80%

# 重新计算潮流
gen_results = runpf(case)

# 比较结果
println("发电机 $(case.gen[gen_idx, GEN_BUS]) 原始出力: $(original_pg) MW")
println("发电机 $(case.gen[gen_idx, GEN_BUS]) 修改后出力: $(case.gen[gen_idx, PG]) MW")
println("平衡节点原始出力: $(results.gen[findall(x -> x == 3, results.bus[:, BUS_TYPE])[1], PG]) MW")
println("平衡节点修改后出力: $(gen_results.gen[findall(x -> x == 3, gen_results.bus[:, BUS_TYPE])[1], PG]) MW")

# 恢复原始发电机数据
case.gen = original_gen
```

## 时间序列潮流计算

### 基本时间序列分析

执行一天的时间序列潮流计算：

```julia
# 准备负荷数据文件
# load_data.xlsx 应包含24小时的负荷数据

# 执行时间序列潮流计算
ts_results = runtdpf(case, "load_data.xlsx")

# 检查计算是否成功
success_count = sum(ts_results.success)
println("成功计算的时间点数量: $success_count / $(length(ts_results.time_points))")

# 绘制总负荷随时间的变化
time_hours = ts_results.time_points
total_load = [sum(result.bus[:, PD]) for result in ts_results.bus_results]

plot(time_hours, total_load,
     title="系统总负荷随时间变化",
     xlabel="时间(小时)",
     ylabel="总负荷(MW)",
     marker=:circle,
     legend=false)
```

### 包含可再生能源的时间序列分析

执行包含光伏发电的时间序列潮流计算：

```julia
# 准备光照强度数据文件
# irradiance_data.xlsx 应包含24小时的光照强度数据

# 执行包含光照数据的时间序列潮流计算
pv_results = runtdpf(case, "load_data.xlsx", "irradiance_data.xlsx")

# 分析光伏发电输出
plot_time_day = 1  # 分析第一天的数据

# 绘制电压分析图
plot_voltage_time_series(pv_results, case, plot_time_day)

# 绘制系统损耗分析图
plot_losses_time_series(pv_results, case, plot_time_day)
```

### 动态经济调度

执行考虑电价的动态经济调度：

```julia
# 准备电价数据文件
# price_data.xlsx 应包含24小时的电价数据

# 设置动态调度选项
dispatch_opt = Dict("mode" => "dynamic_dispatch")

# 执行动态经济调度计算
dispatch_results = runtdpf(case, "load_data.xlsx", "irradiance_data.xlsx", "price_data.xlsx", dispatch_opt)

# 分析调度结果
plot_PD_time_series(dispatch_results, case, 1)

# 分析功率流越限情况
flow_violations = plot_flow_violations(dispatch_results, case, 1, 3.0)
```

## 辅助功能应用

### 内部/外部编号转换

处理系统中的编号转换：

```julia
# 转换为内部编号
internal_case, i2e = ext2int(case)

# 执行潮流计算
internal_results = runpf(internal_case)

# 转换回外部编号
external_results = int2ext(i2e, internal_results)

# 验证转换结果
println("内部编号的母线数量: $(size(internal_results.bus, 1))")
println("外部编号的母线数量: $(size(external_results.bus, 1))")
```

### 孤岛识别与提取

识别系统中的孤岛：

```julia
# 创建包含孤岛的测试系统
test_case = copy(case)

# 断开某些线路以创建孤岛
disconnect_branches = [1, 5, 10]
for br_idx in disconnect_branches
    test_case.branch[br_idx, BR_STATUS] = 0
end

# 识别孤岛
islands = find_islands(test_case)
println("识别到 $(length(islands)) 个孤岛")

# 提取最大的孤岛
main_island = extract_islands(test_case, islands[1])
println("主孤岛包含 $(size(main_island.bus, 1)) 个母线")

# 对主孤岛执行潮流计算
island_results = runpf(main_island)
```

## 混合交直流系统分析

### 重新编号混合系统

处理混合交直流系统：

```julia
# 加载混合交直流系统案例
hybrid_case = load_case("hybrid_case.m")

# 重新编号混合系统
renumbered_case = renumber_hybrid_system(hybrid_case)

# 执行混合系统潮流计算
hybrid_results = runpf(renumbered_case)
```

## 结论

以上示例展示了 PowerFlow 库的基本用法，包括潮流计算、结果分析、时间序列分析和辅助功能。这些示例可以作为开发更复杂应用的起点，用户可以根据自己的需求进行扩展和定制。

要了解更多高级功能，请参考其他示例文档和完整的 API 参考文档。