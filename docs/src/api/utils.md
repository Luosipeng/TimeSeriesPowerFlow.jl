# Utils 模块 API文档

Utils模块提供了电力系统分析中的各种辅助功能，包括数据转换、拓扑分析、孤岛提取等功能。这些工具函数支持PowerFlow模块的核心功能，帮助处理电力系统数据的各种特殊情况。

## 模块概述

Utils模块包含了电力系统数据处理的辅助工具，主要涉及以下方面：
- 内部/外部编号转换
- 孤岛识别与提取
- 拓扑分析
- 数据提取与转换
- 虚拟节点处理

## 主要功能

### 内部/外部编号转换

```julia
ext2int(jpc)
```
将外部编号转换为内部连续编号。
- `jpc`: 电力系统案例数据
- 返回: 转换后的案例数据和内部到外部的映射

```julia
int2ext(i2e, results)
```
将内部编号转换回外部编号。
- `i2e`: 内部到外部的映射
- `results`: 使用内部编号的结果数据
- 返回: 使用外部编号的结果数据

### 孤岛识别与提取

```julia
find_islands(jpc)
```
在电力系统中识别电气孤岛。
- `jpc`: 电力系统案例数据
- 返回: 孤岛信息，包括每个孤岛中的节点列表

```julia
extract_islands(jpc)
```
从电力系统中提取各个孤岛的子系统。
- `jpc`: 电力系统案例数据
- 返回: 包含各个孤岛子系统的列表

### 拓扑分析

```julia
create_node_mapping(case)
```
为节点创建编号映射，检测重复节点。
- `case`: JuliaPowerCase格式的电力系统数据
- 返回: 节点名称到ID的映射字典

```julia
resolve_node_mapping(node_id, node_merge_map)
```
递归解析节点映射，确保多层合并正确处理。
- `node_id`: 节点ID
- `node_merge_map`: 节点合并映射
- 返回: 解析后的节点ID

### 数据提取与转换

```julia
extract_data(worksheet, start_row, end_row, start_col, end_col)
```
从Excel工作表中提取数据。
- `worksheet`: Excel工作表对象
- `start_row`, `end_row`: 起始和结束行
- `start_col`, `end_col`: 起始和结束列
- 返回: 提取的数据数组

```julia
JuliaPowerCase2Jpc(case)
```
将JuliaPowerCase格式转换为JPC格式。
- `case`: JuliaPowerCase格式的电力系统数据
- 返回: JPC格式的电力系统数据

```julia
JuliaPowerCase2Jpc_3ph(case)
```
将JuliaPowerCase格式转换为三相JPC格式。
- `case`: JuliaPowerCase格式的电力系统数据
- 返回: 三相JPC格式的电力系统数据

### 虚拟节点处理

```julia
merge_virtual_nodes(case)
```
合并虚拟节点两侧的节点，去除虚拟节点和虚拟连接。
- `case`: JuliaPowerCase格式的电力系统数据
- 返回: 更新后的case对象，不含虚拟节点和虚拟连接

### 索引常量定义

```julia
idx_bus()
```
返回母线数据列的索引常量。

```julia
idx_gen()
```
返回发电机数据列的索引常量。

```julia
idx_branch()
```
返回支路数据列的索引常量。

```julia
idx_brch()
```
返回支路数据列的索引常量（兼容性别名）。

```julia
idx_cost()
```
返回成本数据列的索引常量。

```julia
idx_dcline()
```
返回直流线路数据列的索引常量。

## 数据结构

### JPC
表示标准电力系统案例数据的结构。
```julia
mutable struct JPC
    baseMVA::Float64
    busAC::Matrix{Float64}
    genAC::Matrix{Float64}
    branchAC::Matrix{Float64}
    loadAC::Matrix{Float64}
    # 其他字段...
end
```

### JPC_3ph
表示三相电力系统案例数据的结构。
```julia
mutable struct JPC_3ph
    baseMVA::Float64
    basef::Float64
    busAC::Matrix{Float64}
    genAC::Matrix{Float64}
    branchAC::Matrix{Float64}
    loadAC::Matrix{Float64}
    # 三相特定字段...
end
```

### MicrogridPlanningProblem
表示微电网规划问题的结构。
```julia
mutable struct MicrogridPlanningProblem
    # 常数
    bigM::Float64
    
    # 规划参数
    nPV::Int
    nBESS::Int
    DeltaPV::Float64
    DeltaBESS::Float64
    # 其他字段...
end
```

## 使用示例

### 内部/外部编号转换
```julia
using PowerFlow

# 加载案例数据
case = load_case("case14.m")

# 转换为内部编号
internal_case, i2e = ext2int(case)

# 执行潮流计算
results = runpf(internal_case)

# 转换回外部编号
external_results = int2ext(i2e, results)
```

### 孤岛识别与处理
```julia
using PowerFlow

# 加载案例数据
case = load_case("case_with_islands.m")

# 识别孤岛
islands_info = find_islands(case)
println("系统中包含 $(length(islands_info)) 个孤岛")

# 提取孤岛子系统
island_subsystems = extract_islands(case)

# 对每个孤岛单独进行潮流计算
for (i, island) in enumerate(island_subsystems)
    println("计算孤岛 $i 的潮流...")
    results = runpf(island)
    # 处理结果...
end

# 合并结果
merged_results = merge_results(results_list, islands_info)
```

### 虚拟节点处理
```julia
using PowerFlow

# 加载包含虚拟节点的案例
case = load_case("case_with_virtual_nodes.m")

# 合并虚拟节点
case_without_virtual = merge_virtual_nodes(case)

# 执行潮流计算
results = runpf(case_without_virtual)
```

### 数据格式转换
```julia
using PowerFlow

# 加载JuliaPowerCase格式的数据
case = load_julia_power_case("my_system.json")

# 转换为JPC格式
jpc = JuliaPowerCase2Jpc(case)

# 执行潮流计算
results = runpf(jpc)

# 对于三相系统
jpc_3ph = JuliaPowerCase2Jpc_3ph(case)
results_3ph = runupf(jpc, jpc_3ph, gs_eg, bs_eg, opt)
```

## 注意事项

1. 在处理大型电网时，孤岛识别和提取功能可能需要较长时间，建议对大系统进行优化。
2. 虚拟节点合并可能会改变系统拓扑，请确保合并后的系统仍然符合预期。
3. 内部/外部编号转换在处理结果时非常重要，尤其是将结果与原始数据对比时。
4. 拓扑分析功能可用于检测系统中的潜在问题，如重复节点或断开的连接。

## 参考文献

1. Power System Analysis, J.J. Grainger and W.D. Stevenson Jr.
2. Graph Theory with Applications to Engineering and Computer Science, N. Deo.
3. MATPOWER: A MATLAB Power System Simulation Package, R.D. Zimmerman, C.E. Murillo-Sánchez, and R.J. Thomas.