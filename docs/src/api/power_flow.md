# PowerFlow模块 API文档

PowerFlow模块提供了电力系统潮流计算的完整功能集，包括交流潮流计算、直流潮流计算、不平衡潮流计算以及相关的辅助功能。本模块支持常规CPU计算和GPU加速计算。

## 模块概述

PowerFlow模块实现了多种潮流计算算法，包括牛顿-拉夫森法、直流潮流计算等，并提供了丰富的电网模型构建和结果处理功能。该模块可用于分析平衡和不平衡电力系统的稳态运行特性。

## 主要功能

### 潮流计算主函数

#### 交流潮流计算

```julia
runpf(case, opt=Dict())
```
执行交流潮流计算的主函数。
- `case`: 电力系统案例数据
- `opt`: 可选参数设置
- 返回: 包含潮流计算结果的数据结构

#### 直流潮流计算

```julia
rundcpf(case, opt=Dict())
```
执行直流潮流计算的主函数。
- `case`: 电力系统案例数据
- `opt`: 可选参数设置
- 返回: 包含直流潮流计算结果的数据结构

#### 不平衡潮流计算

```julia
runupf(case, jpc_3ph, gs_eg, bs_eg, opt)
```
执行不平衡潮流计算的主函数。
- `case`: 平衡系统案例数据
- `jpc_3ph`: 三相系统数据
- `gs_eg`, `bs_eg`: 接地导纳参数
- `opt`: 可选参数设置
- 返回: 包含不平衡潮流计算结果的数据结构

#### 快速解耦潮流计算

```julia
runhpf(jpc, opt)
```
执行快速解耦潮流计算的主函数。
- `jpc`: 电力系统案例数据
- `opt`: 可选参数设置
- 返回: 包含潮流计算结果的数据结构

### 求解器函数

#### 牛顿-拉夫森法求解器

```julia
newtonpf(baseMVA, bus, gen, branch, Ybus, Yf, Yt, Vm, Va, ref, pv, pq, mpopt)
```
使用牛顿-拉夫森法求解交流潮流方程。
- `baseMVA`: 基准功率
- `bus`, `gen`, `branch`: 电力系统组件数据
- `Ybus`, `Yf`, `Yt`: 节点导纳矩阵和支路导纳矩阵
- `Vm`, `Va`: 初始电压幅值和相角
- `ref`, `pv`, `pq`: 节点类型索引
- `mpopt`: 求解器选项
- 返回: 求解结果

```julia
newtonpf_gpu(baseMVA, bus, gen, branch, Ybus, Yf, Yt, Vm, Va, ref, pv, pq, mpopt)
```
使用GPU加速的牛顿-拉夫森法求解交流潮流方程。

#### 直流潮流求解器

```julia
newtondcpf(baseMVA, bus, gen, branch, B, Pbus, Va0, ref, pv, pq, mpopt)
```
使用牛顿法求解直流潮流方程。

```julia
newtondcpf_sp(baseMVA, bus, gen, branch, B, Pbus, Va0, ref, pv, pq, mpopt)
```
使用稀疏矩阵技术的牛顿法求解直流潮流方程。

```julia
dcpf(B, Pbus, Va0, ref, pv, pq, mpopt)
```
基本直流潮流求解函数。

### 电网模型构建函数

#### 节点导纳矩阵构建

```julia
makeYbus(baseMVA, bus, branch)
```
构建节点导纳矩阵和支路导纳矩阵。
- `baseMVA`: 基准功率
- `bus`, `branch`: 节点和支路数据
- 返回: `Ybus`, `Yf`, `Yt` 导纳矩阵

#### 直流潮流B矩阵构建

```julia
makeBdc(baseMVA, bus, branch)
```
构建直流潮流计算所需的B矩阵。
- `baseMVA`: 基准功率
- `bus`, `branch`: 节点和支路数据
- 返回: B矩阵和相关数据

#### 功率注入计算

```julia
makeSbus(baseMVA, bus, gen, Vm, load, pvarray; dc=false, Sg=nothing, return_derivative=false)
```
计算节点功率注入。
- `baseMVA`: 基准功率
- `bus`, `gen`: 节点和发电机数据
- `Vm`: 电压幅值
- `load`: 负荷数据
- `pvarray`: 光伏阵列数据
- 返回: 节点功率注入向量

```julia
makeSbus_gpu(baseMVA, bus_gpu, gen_gpu, gen, Vm, load_gpu, load, pvarray_gpu, pvarray; dc=false, Sg=nothing, return_derivative=false)
```
使用GPU计算节点功率注入。

### 电网元素构建函数

```julia
calculate_bus(net, jpc, sequence, slack_bus, opt)
```
构建节点数据结构。

```julia
calculate_line_parameter(net, jpc, sequence, opt)
```
计算线路参数。

```julia
build_gen(net, jpc)
```
构建发电机数据结构。

### 节点类型处理

```julia
bustypes(bus, gen)
```
确定节点类型（参考节点、PV节点、PQ节点）。

```julia
dcbustypes(bus, gen)
```
确定直流潮流计算的节点类型。

### 结果处理函数

```julia
pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt)
```
处理交流潮流计算结果。

```julia
dcpfsoln(baseMVA, bus, gen, branch, Bbus, Bf, Pbusinj, V, ref, pv, pq, mpopt)
```
处理直流潮流计算结果。

```julia
process_result(results, isolated, file_path)
```
处理潮流计算结果并可选保存到文件。

```julia
merge_results(results, isolated)
```
合并孤岛的潮流计算结果。

```julia
write_system_summary(f, mpc, area, isolated)
```
生成系统摘要报告。

### 辅助函数

```julia
dSbus_dV(Ybus, V)
```
计算功率对电压的雅可比矩阵。

```julia
eliminate_element(jpc)
```
消除电网中的特定元素。

```julia
total_load(bus, load)
```
计算系统总负荷。

```julia
process_pv_acsystem(pv_acsystem, pv_acsystem_bus, pv_acsystem_gen, opt)
```
处理光伏交流系统数据。

```julia
compare_voltage_results(results, case, reference_file; tolerance_mag=1e-4, tolerance_ang=1e-3)
```
比较潮流计算结果与参考结果。

### GPU加速计算函数

```julia
gpu_gmres(A, b, x0; tol=1e-6, maxiter=100, restart=20, M=nothing)
```
使用GPU加速的GMRES求解线性方程组。

```julia
makeSdzip(Sbusd, V, bus, opt)
```
计算ZIP负荷模型的功率注入。

```julia
makeSdzip_gpu(Sbusd, V, bus, opt)
```
使用GPU计算ZIP负荷模型的功率注入。

## 设置与预处理

```julia
settings(opt)
```
设置潮流计算的参数。

```julia
runprepf(mpc, opt)
```
执行潮流计算的预处理。

```julia
dc_preprocess(mpc, opt)
```
执行直流潮流计算的预处理。

## 数据结构

### Sd
表示ZIP负荷模型的数据结构。
```julia
mutable struct Sd
    z::Vector{ComplexF64}  # 阻抗负荷
    i::Vector{ComplexF64}  # 电流负荷
    p::Vector{ComplexF64}  # 功率负荷
end
```

### Sd_gpu
表示GPU计算中的ZIP负荷模型数据结构。
```julia
mutable struct Sd_gpu
    z::CuVector{ComplexF64}  # GPU上的阻抗负荷
    i::CuVector{ComplexF64}  # GPU上的电流负荷
    p::CuVector{ComplexF64}  # GPU上的功率负荷
end
```

## 使用示例

### 基本交流潮流计算
```julia
using PowerFlow

# 加载案例数据
case = load_case("case14.m")

# 执行潮流计算
results = runpf(case)

# 分析结果
println("收敛状态: ", results.success)
println("迭代次数: ", results.iterations)
```

### 直流潮流计算
```julia
using PowerFlow

# 加载案例数据
case = load_case("case14.m")

# 执行直流潮流计算
results = rundcpf(case)

# 分析结果
println("收敛状态: ", results.success)
```

### GPU加速计算
```julia
using PowerFlow
using CUDA

# 确保CUDA可用
if CUDA.functional()
    # 加载案例数据
    case = load_case("case14.m")
    
    # 设置GPU选项
    opt = Dict("PF" => Dict("alg" => "NR", "use_gpu" => true))
    
    # 执行GPU加速的潮流计算
    results = runpf(case, opt)
    
    println("使用GPU计算，迭代次数: ", results.iterations)
else
    println("CUDA不可用，无法执行GPU加速计算")
end
```

## 注意事项

1. GPU加速计算需要安装CUDA工具包和兼容的GPU硬件。
2. 大型电网模型计算时，建议使用稀疏矩阵技术以提高计算效率。
3. 对于不平衡系统，请使用`runupf`函数并提供三相系统数据。
4. 结果比较功能可用于验证计算结果与其他软件（如ETAP）的一致性。

## 参考文献

1. Power System Analysis, J.J. Grainger and W.D. Stevenson Jr.
2. Power System Load Flow Analysis, J.D. Glover, M.S. Sarma, and T.J. Overbye.
3. MATPOWER: A MATLAB Power System Simulation Package, R.D. Zimmerman, C.E. Murillo-Sánchez, and R.J. Thomas.