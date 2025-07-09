# 交直流混合潮流迭代算法技术说明

## 1. 概述

本文档描述了一个用于求解交直流混合电力系统潮流的迭代算法。该算法通过交替求解交流和直流潮流，并根据换流器的不同工作模式更新功率分配，最终达到系统收敛。

## 2. 算法架构

### 2.1 主要函数

- `runhpf_iteration(jpc, opt)`: 主迭代函数
- `update_mode_X_converters!()`: 各种换流器模式的更新函数
- `update_dc_load!()`: 直流侧负荷更新函数
- `update_ac_load!()`: 交流侧负荷更新函数

### 2.2 换流器工作模式

算法支持5种换流器工作模式：

| 模式 | 控制变量 | 说明 |
|------|----------|------|
| 1 | constant δs、Us | 恒定相角差和交流电压 |
| 2 | constant Ps、Qs | 恒定有功和无功功率（默认模式） |
| 3 | constant Ps、Us | 恒定有功功率和交流电压 |
| 4 | constant Udc、Qs | 恒定直流电压和无功功率 |
| 5 | constant Udc、Us | 恒定直流电压和交流电压 |

## 3. 算法流程

### 3.1 初始化阶段

```julia
# 设置迭代参数
max_iterations = 20          # 最大迭代次数
convergence_tolerance = 1e-4 # 收敛容差

# 保存初始负荷数据，避免无限叠加
initial_busAC = deepcopy(jpc.busAC)
initial_busDC = deepcopy(jpc.busDC)
initial_loadAC = deepcopy(jpc.loadAC)
initial_loadDC = deepcopy(jpc.loadDC)
```

**关键设计要点：**
- 保存初始负荷状态，确保每次迭代都从原始负荷开始计算
- 避免功率值在迭代过程中无限累积

### 3.2 主迭代循环

每次迭代包含以下步骤：

#### 步骤1：重置负荷状态
```julia
# 重置为初始负荷状态，避免累积叠加
result_jpc.busAC = deepcopy(initial_busAC)
result_jpc.busDC = deepcopy(initial_busDC)
result_jpc.loadAC = deepcopy(initial_loadAC)
result_jpc.loadDC = deepcopy(initial_loadDC)
```

#### 步骤2：准备系统数据
- 创建独立的交流系统数据结构 `jpc1`
- 创建独立的直流系统数据结构 `jpc2`

#### 步骤3：更新换流器功率
根据换流器的工作模式，调用相应的更新函数：
- `update_mode_1_converters!()` - 处理模式1换流器
- `update_mode_3_converters!()` - 处理模式3换流器  
- `update_mode_4_converters!()` - 处理模式4换流器
- `update_mode_5_converters!()` - 处理模式5换流器

#### 步骤4：潮流计算
```julia
# 交流潮流计算
if !isempty(jpc1.busAC)
    jpc1 = PowerFlow.runpf(jpc1, opt)
end

# 直流潮流计算
if !isempty(jpc2.busDC)
    jpc2 = PowerFlow.rundcpf(jpc2, opt)
end
```

#### 步骤5：收敛性检查
- 比较当前迭代与上一次迭代的功率值差异
- 如果最大差值小于收敛容差，则认为收敛
- 更新结果数据结构

## 4. 换流器功率更新机制

### 4.1 功率计算原理

不同模式的换流器采用不同的功率计算方法：

**模式1和模式3（交流侧控制）：**
```julia
P_ac = -jpc1.genAC[gen_row, PG]  # 从交流发电机获取功率
Q_ac = -jpc1.genAC[gen_row, QG]

# 根据效率计算直流侧功率
if P_ac < 0
    P_dc = -P_ac/conv[CONV_EFF]  # 整流模式
else
    P_dc = -P_ac * conv[CONV_EFF]  # 逆变模式
end
```

**模式4和模式5（直流侧控制）：**
```julia
P_dc = -jpc2.genDC[gen_row, PG]  # 从直流发电机获取功率

# 根据效率计算交流侧功率
if P_dc < 0
    P_ac = -P_dc/conv[CONV_EFF]  # 整流模式
else
    P_ac = -P_dc * conv[CONV_EFF]  # 逆变模式
end
```

### 4.2 负荷更新策略

#### 直流侧负荷更新
```julia
function update_dc_load!(result_jpc, jpc2, dc_bus_id, P_dc)
    # 更新母线功率
    result_jpc.busDC[dc_bus_rows, PD] .+= P_dc
    jpc2.busDC[dc_bus_rows_jpc2, PD] .+= P_dc
    
    # 更新或创建负荷
    if !isempty(dc_load_rows)
        result_jpc.loadDC[dc_load_rows, LOAD_PD] .+= P_dc
    else
        # 创建虚拟负荷
        new_load_row = create_virtual_load(dc_bus_id, P_dc)
        result_jpc.loadDC = vcat(result_jpc.loadDC, new_load_row)
    end
end
```

**关键特性：**
- 同时更新 `result_jpc` 和 `jpc2` 的负荷数据
- 如果指定母线没有负荷，自动创建虚拟负荷
- 确保潮流计算使用正确的负荷数据

## 5. 收敛性控制

### 5.1 收敛判据

算法采用功率差值作为收敛判据：

```julia
max_diff = 0.0
for (key, value) in current_power_values
    if haskey(prev_power_values, key)
        diff = abs(value - prev_power_values[key])
        max_diff = max(max_diff, diff)
    end
end

if max_diff < convergence_tolerance
    # 收敛
    break
end
```

### 5.2 收敛参数

- **收敛容差**：`1e-4`
- **最大迭代次数**：`20`
- **监控变量**：换流器的交流和直流功率

## 6. 数据结构设计

### 6.1 换流器数据结构

换流器数据通过 `idx_conv()` 函数获取索引常量：

```julia
CONV_ACBUS, CONV_DCBUS, CONV_INSERVICE, CONV_P_AC, 
CONV_Q_AC, CONV_P_DC, CONV_EFF, CONV_MODE = idx_conv()
```

### 6.2 功率值存储

使用字典存储功率值用于收敛判断：

```julia
current_power_values = Dict()
current_power_values["conv_$(conv_id)_P_ac"] = P_ac
current_power_values["conv_$(conv_id)_Q_ac"] = Q_ac
current_power_values["conv_$(conv_id)_P_dc"] = P_dc
```

## 7. 错误处理

### 7.1 潮流计算失败处理

```julia
if !isempty(jpc2.busDC)
    if !(jpc1.success && jpc2.success)
        @warn "潮流计算在第 $iter 次迭代中失败"
        break
    end
elseif !jpc1.success
    @warn "交流潮流计算在第 $iter 次迭代中失败"
    break
end
```

### 7.2 收敛性监控

- 输出收敛信息：`@info "迭代收敛于第 $iter 次迭代"`
- 警告未收敛情况：`@warn "达到最大迭代次数但未收敛"`

## 8. 算法优势

### 8.1 数值稳定性
- 每次迭代重置负荷状态，避免数值累积误差
- 独立的交直流系统数据结构，避免相互干扰

### 8.2 灵活性
- 支持多种换流器工作模式
- 自动处理虚拟负荷创建
- 可配置的收敛参数

### 8.3 鲁棒性
- 完善的错误处理机制
- 潮流计算失败时的优雅退出
- 详细的日志输出

## 9. 使用示例

```julia
# 基本用法
result = runhpf_iteration(jpc, opt)

# 检查结果
if result.success
    println("混合潮流计算成功")
    println("交流迭代次数: ", result.iterationsAC)
    println("直流迭代次数: ", result.iterationsDC)
else
    println("混合潮流计算失败")
end
```

## 10. 注意事项

1. **函数依赖性**：算法依赖 `idx_conv()` 函数提供换流器索引常量，该函数不可修改
2. **数据一致性**：确保输入数据中换流器、母线、发电机数据的一致性
3. **收敛性调优**：根据具体系统特性，可能需要调整收敛容差和最大迭代次数
4. **内存使用**：算法使用大量 `deepcopy` 操作，对于大型系统需要注意内存使用

## 11. 扩展性

算法设计具有良好的扩展性：

- 可以轻松添加新的换流器工作模式
- 支持自定义收敛判据
- 可以集成其他类型的电力电子设备模型