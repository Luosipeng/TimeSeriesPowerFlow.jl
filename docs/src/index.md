# DistributionPowerFlow 库文档

欢迎使用 PowerFlow 库，这是一个用于电力系统潮流计算和分析的高性能 Julia 工具包。本库提供了全面的功能，从基础的潮流计算到高级的时间序列分析和可视化，旨在满足研究人员、工程师和教育工作者的需求。

## 功能概述

DistributionPowerFlow 库提供以下核心功能：

- **基础潮流计算**：支持交流潮流、直流潮流和线性分布潮流计算
- **高级算法**：实现了牛顿-拉夫森法、快速解耦法、高斯-赛德尔法等多种求解算法
- **时间序列分析**：支持基于时间序列的潮流计算和动态分析
- **可再生能源集成**：内置对光伏发电、风力发电等可再生能源的建模支持
- **经济调度**：提供动态经济调度和考虑网络约束的最优潮流计算
- **混合交直流系统**：支持混合交直流系统的建模与计算
- **数据可视化**：丰富的可视化工具，用于结果展示和系统分析
- **大规模系统支持**：针对大规模电力系统的高性能计算优化

## 安装指南

### 系统要求

- Julia 1.6 或更高版本
- 支持的操作系统：Windows、macOS、Linux

### 安装步骤

使用 Julia 的包管理器安装 DistributionPowerFlow 库：

```julia
using Pkg
Pkg.add("PowerFlow")
```

要安装完整的工具包（包括时间序列分析和可视化模块）：

```julia
using Pkg
Pkg.add(["PowerFlow", "TimeSeriesPowerFlow", "PowerFlowVisualization"])
```

## 快速入门

以下是一个简单的例子，展示如何加载测试系统并执行潮流计算：

```julia
using PowerFlow

# 加载 IEEE 14 节点测试系统
case = load_case("case14.m")

# 执行交流潮流计算
results = runpf(case)

# 检查结果
if results.success
    println("潮流计算成功收敛！")
    println("系统总发电: $(sum(results.gen[:, PG])) MW")
    println("系统总负荷: $(sum(results.bus[:, PD])) MW")
    println("系统总损耗: $(sum(results.branch[:, PL])) MW")
else
    println("潮流计算未收敛，请检查系统参数")
end
```

## 库结构

PowerFlow 库由以下主要模块组成：

- **PowerFlow**：核心模块，提供基础潮流计算功能
- **TimeSeriesPowerFlow**：时间序列分析模块，用于动态系统分析
- **ComponentModel**：电力系统组件建模模块
- **Utils**：实用工具和辅助函数
- **Visualization**：结果可视化和图形生成工具

## 文档导航

- [基础教程](basic.md)：入门级教程，介绍基本概念和功能
- [PowerFlow 模块](powerflow.md)：核心潮流计算模块的详细文档
- [TimeSeriesPowerFlow 模块](timeseriesflow.md)：时间序列分析模块的详细文档
- [可视化工具](visualization.md)：可视化功能的详细文档
- [API 参考](api.md)：完整的 API 参考文档
- [示例](examples/)：各种应用场景的示例代码

## 示例概览

PowerFlow 库提供了丰富的示例，帮助用户快速上手：

- **基础潮流计算**：展示如何执行和分析基本的潮流计算
- **时间序列分析**：演示如何进行基于时间序列的系统分析
- **可再生能源集成**：说明如何在系统中集成可再生能源
- **经济调度**：展示如何执行经济调度和最优潮流计算
- **混合交直流系统**：演示如何处理混合交直流系统
- **大规模系统**：针对大规模系统的计算优化示例
- **可视化应用**：展示各种可视化功能的应用

## 支持的文件格式

PowerFlow 库支持多种电力系统数据格式：

- **MATPOWER 格式**：兼容 MATPOWER 的 `.m` 文件格式
- **IEEE 通用格式**：支持 IEEE 通用数据格式
- **PSS/E 格式**：支持 PSS/E 的 `.raw` 文件格式
- **自定义 JSON 格式**：支持基于 JSON 的自定义数据格式
- **Excel 格式**：支持从 Excel 文件导入数据

## 性能优化

PowerFlow 库针对性能进行了多方面优化：

- **并行计算**：利用多核处理器加速计算
- **稀疏矩阵**：高效的稀疏矩阵算法
- **GPU 加速**：支持 GPU 加速大规模系统计算
- **自适应算法**：根据系统特性自动选择最优算法

## 贡献指南

我们欢迎社区贡献，无论是功能改进、bug 修复还是文档更新。请参考以下步骤：

1. Fork 项目仓库
2. 创建您的功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交您的更改 (`git commit -m 'Add some amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 打开 Pull Request

## 引用

如果您在研究中使用了 PowerFlow 库，请引用：

```
@software{PowerFlow2023,
  author = {PowerFlow Development Team},
  title = {PowerFlow: A High-Performance Power Flow Analysis Library in Julia},
  year = {2023},
  url = {https://github.com/powerflow/PowerFlow.jl}
}
```

## 许可证

PowerFlow 库采用 MIT 许可证。详情请参见 [LICENSE](LICENSE) 文件。

## 联系方式

- **项目主页**：[https://github.com/powerflow/PowerFlow.jl](https://github.com/powerflow/PowerFlow.jl)
- **文档**：[https://powerflow.github.io/docs](https://powerflow.github.io/docs)
- **问题报告**：[https://github.com/powerflow/PowerFlow.jl/issues](https://github.com/powerflow/PowerFlow.jl/issues)
- **邮件列表**：powerflow-users@googlegroups.com

## 致谢

PowerFlow 库的开发得到了以下机构和项目的支持：

- 国家自然科学基金委员会
- 国家重点研发计划
- 电力系统仿真与优化实验室
- Julia 开源社区

我们也感谢所有为项目做出贡献的开发者和用户。

## 更新日志

### v1.2.0 (2023-06-15)

- 新增混合交直流系统支持
- 改进时间序列分析性能
- 增强可视化功能
- 修复多个 bug

### v1.1.0 (2023-03-10)

- 新增经济调度模块
- 新增可再生能源集成支持
- 改进大规模系统计算性能
- 更新文档和示例

### v1.0.0 (2022-12-01)

- 首次正式发布
- 实现核心潮流计算功能
- 提供基础时间序列分析
- 包含基本可视化工具

查看完整的更新历史，请访问 [CHANGELOG.md](CHANGELOG.md)。