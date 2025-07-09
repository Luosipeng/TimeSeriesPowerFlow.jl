# TimeSeriesPowerFlow 模块 API文档

TimeSeriesPowerFlow模块提供了电力系统时间序列潮流计算和分析的完整功能集，包括时间序列潮流计算、动态调度、结果可视化以及违例分析等功能。该模块特别适用于含有可再生能源和负荷变化的电力系统分析。

## 模块概述

TimeSeriesPowerFlow模块实现了多种时间序列分析算法，包括线性分布潮流、最大功率点跟踪、动态调度等，并提供了丰富的结果可视化和分析工具。该模块可用于分析电力系统在一天或更长时间内的动态运行特性。

## 主要功能

### 时间序列潮流计算

#### 时间序列潮流计算主函数

```julia
runtdpf(case, load_file, irradiance_file=nothing, price_file=nothing, opt=Dict())
```
执行时间序列潮流计算的主函数。
- `case`: 电力系统案例数据
- `load_file`: 负荷数据文件路径
- `irradiance_file`: 光照强度数据文件路径（可选）
- `price_file`: 电价数据文件路径（可选）
- `opt`: 可选参数设置
- 返回: 包含时间序列潮流计算结果的数据结构

#### 单日运行函数

```julia
run_single_day(old_jpc, opt, day_load_matrix, day_price_line, day_irradiance_line)
```
执行单日时间序列潮流计算。
- `old_jpc`: 电力系统案例数据
- `opt`: 可选参数设置
- `day_load_matrix`: 日负荷数据矩阵
- `day_price_line`: 日电价数据
- `day_irradiance_line`: 日光照强度数据
- 返回: 包含单日潮流计算结果的数据结构

### 优化调度函数

#### 线性分布潮流

```julia
runlindistflow(new_jpc, opt, Cld_ac, Cld_dc, loadAC_PD, loadAC_QD, loadDC_PD, genAC_PG, Cgen_ac, Cconv_ac, Cconv_dc, Pij_limit, Pij_limit_dc, solver)
```
执行线性分布潮流计算。
- `new_jpc`: 电力系统案例数据
- `opt`: 可选参数设置
- `Cld_ac`, `Cld_dc`: 交直流负荷关联矩阵
- `loadAC_PD`, `loadAC_QD`: 交流有功和无功负荷
- `loadDC_PD`: 直流有功负荷
- `genAC_PG`: 交流发电机有功出力
- `Cgen_ac`: 发电机关联矩阵
- `Cconv_ac`, `Cconv_dc`: 换流器关联矩阵
- `Pij_limit`, `Pij_limit_dc`: 交直流线路功率限制
- `solver`: 优化求解器
- 返回: 线性分布潮流计算结果

#### 动态调度

```julia
run_dynamic_dispatch(new_jpc, Cld_ac, Cld_dc, loadAC_PD, loadAC_QD, loadDC_PD, genAC_PG, Cgen_ac, Cconv_ac, Cconv_dc, Pij_limit, Pij_limit_dc, price, irradiance, opt)
```
执行动态经济调度计算。
- `new_jpc`: 电力系统案例数据
- `Cld_ac`, `Cld_dc`: 交直流负荷关联矩阵
- `loadAC_PD`, `loadAC_QD`: 交流有功和无功负荷
- `loadDC_PD`: 直流有功负荷
- `genAC_PG`: 交流发电机有功出力
- `Cgen_ac`: 发电机关联矩阵
- `Cconv_ac`, `Cconv_dc`: 换流器关联矩阵
- `Pij_limit`, `Pij_limit_dc`: 交直流线路功率限制
- `price`: 电价数据
- `irradiance`: 光照强度数据
- `opt`: 可选参数设置
- 返回: 动态调度计算结果

#### 最大功率点跟踪

```julia
runmppt(new_jpc, opt, Cld_ac, Cld_dc, loadAC_PD, loadAC_QD, loadDC_PD, genAC_PG, Cgen_ac, Cconv_ac, Cconv_dc, Pij_limit, Pij_limit_dc, irradiance)
```
执行光伏系统最大功率点跟踪计算。
- `new_jpc`: 电力系统案例数据
- `opt`: 可选参数设置
- `Cld_ac`, `Cld_dc`: 交直流负荷关联矩阵
- `loadAC_PD`, `loadAC_QD`: 交流有功和无功负荷
- `loadDC_PD`: 直流有功负荷
- `genAC_PG`: 交流发电机有功出力
- `Cgen_ac`: 发电机关联矩阵
- `Cconv_ac`, `Cconv_dc`: 换流器关联矩阵
- `Pij_limit`, `Pij_limit_dc`: 交直流线路功率限制
- `irradiance`: 光照强度数据
- 返回: 最大功率点跟踪计算结果

#### 自动潮流优化

```julia
runautopf(jpc, opt, mode="lindistflow")
```
执行自动潮流优化计算。
- `jpc`: 电力系统案例数据
- `opt`: 可选参数设置
- `mode`: 优化模式，默认为"lindistflow"
- 返回: 自动潮流优化计算结果

### 结果可视化函数

#### 功率流越限分析

```julia
plot_flow_violations(results, case, time_day, flow_limit=3.0, plot_type="summary", flow_direction="max", show_plot=true, save_path=nothing)
```
计算并绘制系统潮流越限的时间序列图。
- `results`: 结果数据集
- `case`: 系统案例
- `time_day`: 天数
- `flow_limit`: 线路功率流上限(MW)，默认为3.0
- `plot_type`: 绘图类型，可选"summary"、"worst"或"all"
- `flow_direction`: 功率流方向，可选"max"、"both"、"forward"或"reverse"
- `show_plot`: 是否显示图形，默认为true
- `save_path`: 保存图形的路径，默认为nothing
- 返回: 越限统计信息和图形对象

#### 系统损耗分析

```julia
plot_losses_time_series(results, case, time_day, plot_type="total", loss_type="active", show_plot=true, save_path=nothing)
```
计算并绘制系统损耗的时间序列图。
- `results`: 结果数据集
- `case`: 系统案例
- `time_day`: 天数
- `plot_type`: 绘图类型，可选"total"或"branch"
- `loss_type`: 损耗类型，可选"active"或"reactive"
- `show_plot`: 是否显示图形，默认为true
- `save_path`: 保存图形的路径，默认为nothing
- 返回: 损耗统计信息和图形对象

#### 有功负荷分析

```julia
plot_PD_time_series(results, case, time_day, show_plot=true, save_path=nothing)
```
绘制有功负荷的时间序列图。
- `results`: 结果数据集
- `case`: 系统案例
- `time_day`: 天数
- `show_plot`: 是否显示图形，默认为true
- `save_path`: 保存图形的路径，默认为nothing
- 返回: 图形对象

#### 电压分析

```julia
plot_voltage_time_series(results, case, time_day, plot_type="magnitude", show_plot=true, save_path=nothing)
```
绘制电压幅值和相角的时间序列图。
- `results`: 结果数据集
- `case`: 系统案例
- `time_day`: 天数
- `plot_type`: 绘图类型，可选"magnitude"或"angle"
- `show_plot`: 是否显示图形，默认为true
- `save_path`: 保存图形的路径，默认为nothing
- 返回: 图形对象

```julia
record_voltage_violation(results, bus_name, case, time_day, bus_type="AC")
```
记录并分析电压越限情况。
- `results`: 结果数据集
- `bus_name`: 母线名称
- `case`: 系统案例
- `time_day`: 天数
- `bus_type`: 母线类型，默认为"AC"
- 返回: 电压越限统计信息和图形对象

### 系统处理函数

```julia
renumber_hybrid_system(jpc)
```
重新编号混合交直流系统。
- `jpc`: 电力系统案例数据
- 返回: 重新编号后的案例数据

## 数据结构

### TimeSeriesResults
表示时间序列潮流计算结果的数据结构。
```julia
struct TimeSeriesResults
    time_points::Vector{Float64}      # 时间点
    time_strings::Vector{String}      # 时间字符串
    bus_results::Vector{Any}          # 母线结果
    branch_results::Vector{Any}       # 支路结果
    gen_results::Vector{Any}          # 发电机结果
    converter_results::Vector{Any}    # 换流器结果
    success::Vector{Bool}             # 计算成功标志
end
```

### OptimizationResults
表示优化计算结果的数据结构。
```julia
struct OptimizationResults
    objective::Float64                # 目标函数值
    V_ac::Vector{Float64}             # 交流电压
    V_dc::Vector{Float64}             # 直流电压
    Pij_ac::Vector{Float64}           # 交流支路功率
    Pij_dc::Vector{Float64}           # 直流支路功率
    Pg_ac::Vector{Float64}            # 交流发电机功率
    Pc_ac::Vector{Float64}            # 交流换流器功率
    Pc_dc::Vector{Float64}            # 直流换流器功率
    success::Bool                     # 计算成功标志
end
```

## 使用示例

### 基本时间序列潮流计算
```julia
using TimeSeriesPowerFlow

# 加载案例数据
case = load_case("case14.m")

# 执行时间序列潮流计算
results = runtdpf(case, "load_data.xlsx")

# 分析结果
plot_voltage_time_series(results, case, 1)
```

### 带有光伏的时间序列潮流计算
```julia
using TimeSeriesPowerFlow

# 加载案例数据
case = load_case("case_with_pv.m")

# 执行时间序列潮流计算，包含光照数据
results = runtdpf(case, "load_data.xlsx", "irradiance_data.xlsx")

# 分析结果
plot_flow_violations(results, case, 1, 3.0)
plot_losses_time_series(results, case, 1)
```

### 动态经济调度
```julia
using TimeSeriesPowerFlow

# 加载案例数据
case = load_case("hybrid_case.m")

# 执行时间序列潮流计算，包含电价和光照数据
opt = Dict("mode" => "dynamic_dispatch")
results = runtdpf(case, "load_data.xlsx", "irradiance_data.xlsx", "price_data.xlsx", opt)

# 分析结果
plot_PD_time_series(results, case, 1)
```

### 电压越限分析
```julia
using TimeSeriesPowerFlow

# 加载案例数据
case = load_case("case14.m")

# 执行时间序列潮流计算
results = runtdpf(case, "load_data.xlsx")

# 分析特定母线的电压越限情况
violation_info = record_voltage_violation(results, "Bus 5", case, 1)
println("电压越限次数: ", violation_info.violation_count)
```

## 注意事项

1. 时间序列数据文件（负荷、光照、电价）应采用Excel格式，并按照特定格式组织。
2. 对于混合交直流系统，请确保使用`renumber_hybrid_system`函数进行正确的编号处理。
3. 线性分布潮流算法适用于辐射状网络，对于网状结构可能存在精度问题。
4. 动态调度和MPPT模式需要提供光照数据才能正确模拟光伏系统的行为。
5. 结果可视化函数支持中文显示，但需要系统安装相应的中文字体。

## 参考文献

1. "Distribution System Analysis and the Future Smart Grid", IEEE Transactions on Industry Applications.
2. "Linear Power Flow in Distribution Networks Considering Host Capability and Network Constraints", IEEE Transactions on Smart Grid.
3. "Time-Series Analysis of Large-Scale PV Integration in Distribution Systems", IEEE Transactions on Power Systems.