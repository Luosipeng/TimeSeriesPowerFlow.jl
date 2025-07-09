# DistributionPowerFlow 库 API 概述文档

DistributionPowerFlow 是一个全面的电力系统分析库，提供了从基础组件建模到高级时间序列分析的完整功能集。该库由四个核心模块组成，每个模块专注于电力系统分析的不同方面。

## 模块结构

DistributionPowerFlow 库由以下四个主要模块组成：

1. **ComponentModel** - 电力系统组件建模
2. **PowerFlow** - 潮流计算核心功能
3. **TimeSeriesPowerFlow** - 时间序列潮流分析
4. **Utils** - 辅助工具和功能

## 主要特性

- 支持平衡和不平衡电力系统分析
- 混合交直流系统建模和计算
- 时间序列潮流计算和动态调度
- 可再生能源集成分析
- 结果可视化和违例分析
- CPU和GPU加速计算

## 模块功能概览

### ComponentModel 模块

ComponentModel 模块提供了电力系统组件的数据结构和基本操作。

**主要组件：**
- `Bus` - 交流母线节点
- `BusDC` - 直流母线节点
- `Line` - 交流线路
- `LineDC` - 直流线路
- `Generator` - 发电机
- `Load` - 负荷
- `Converter` - 换流器

**示例用法：**
```julia
using ComponentModel

# 创建母线
bus = Bus(1, "Bus1", 110.0, "PQ", "A1", "Z1", true, 1.0, 0.0)

# 创建线路
line = Line(1, "Line1", 1, 2, 10.0, 0.1, 0.4, 10.0)
```

### PowerFlow 模块

PowerFlow 模块实现了多种潮流计算算法，是库的核心计算引擎。

**主要功能：**
- `runpf()` - 交流潮流计算
- `rundcpf()` - 直流潮流计算
- `runupf()` - 不平衡潮流计算

**示例用法：**
```julia
using PowerFlow

# 加载案例
case = load_case("case14.m")

# 执行交流潮流计算
results = runpf(case)

# 执行直流潮流计算
dc_results = rundcpf(case)
```

### TimeSeriesPowerFlow 模块

TimeSeriesPowerFlow 模块扩展了基本潮流计算，支持时间序列分析和动态调度。

**主要功能：**
- `runtdpf()` - 时间序列潮流计算
- `run_dynamic_dispatch()` - 动态经济调度
- `runlindistflow()` - 线性分布潮流计算
- 结果可视化函数

**示例用法：**
```julia
using TimeSeriesPowerFlow

# 执行时间序列潮流计算
results = runtdpf(case, "load_data.xlsx", "irradiance_data.xlsx")

# 分析结果
plot_voltage_time_series(results, case, 1)
```

### Utils 模块

Utils 模块提供了各种辅助功能，支持其他模块的操作。

**主要功能：**
- `ext2int()/int2ext()` - 内部/外部编号转换
- `find_islands()/extract_islands()` - 孤岛识别与提取
- `merge_virtual_nodes()` - 虚拟节点处理

**示例用法：**
```julia
using Utils

# 转换为内部编号
internal_case, i2e = ext2int(case)

# 识别孤岛
islands_info = find_islands(case)
```

## 典型工作流程

1. **系统建模**：使用 ComponentModel 模块创建或导入电力系统模型
2. **基本潮流分析**：使用 PowerFlow 模块进行静态潮流计算
3. **时间序列分析**：使用 TimeSeriesPowerFlow 模块进行动态分析
4. **结果处理**：使用 Utils 模块和可视化函数分析结果

## 系统要求

- Julia 1.6 或更高版本
- 支持的操作系统: Windows, macOS, Linux
- 可选: CUDA工具包（用于GPU加速）

## 安装方法

通过Julia的包管理器安装:

```julia
using Pkg
Pkg.add("DistributionPowerFlow")
```

或者在Julia REPL的包管理模式中安装:

```
] add DistributionPowerFlow
```

## 快速入门示例

```julia
# 导入所需模块
using PowerFlow
using TimeSeriesPowerFlow
using Utils

# 加载案例数据
case = load_case("case14.m")

# 执行潮流计算
pf_results = runpf(case)

# 执行时间序列潮流计算
ts_results = runtdpf(case, "load_data.xlsx")

# 分析结果
plot_voltage_time_series(ts_results, case, 1)
plot_losses_time_series(ts_results, case, 1)
```

## 文档和资源

- 完整API文档: [链接到文档]
- 示例库: [链接到示例]
- 问题报告: [链接到问题跟踪器]

## 许可证

DistributionPowerFlow 库在 [许可证名称] 下发布。详情请参阅 LICENSE 文件。