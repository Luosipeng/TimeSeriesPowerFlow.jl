# ComponentModel 模块 API 参考文档

## 安装指南

### 系统要求
- Julia 1.6 或更高版本
- 支持的操作系统: Windows, macOS, Linux

### 安装方法

通过Julia的包管理器安装:

```julia
using Pkg
Pkg.add("ComponentModel")
```

或者在Julia REPL的包管理模式中安装:

```
] add ComponentModel
```

### 导入模块

```julia
using ComponentModel
```

## 基础组件

### Bus
```julia
mutable struct Bus <: AbstractComponent
    index::Int
    name::String
    vn_kv::Float64
    type::String
    area::String
    zone::String
    in_service::Bool
    vm_pu::Float64
    va_degree::Float64
    # 其他属性...
end
```
表示电力系统中的交流母线节点。

属性:
- `index`: 唯一标识符
- `name`: 母线名称
- `vn_kv`: 额定电压(kV)
- `type`: 母线类型(PQ, PV, REF)
- `area`: 所属区域
- `zone`: 所属分区
- `in_service`: 运行状态
- `vm_pu`: 电压幅值(标幺值)
- `va_degree`: 电压相角(度)

### BusDC
```julia
mutable struct BusDC <: AbstractComponent
    index::Int
    name::String
    vn_kv::Float64
    type::String
    area::String
    zone::String
    in_service::Bool
    vm_pu::Float64
    # 其他属性...
end
```
表示电力系统中的直流母线节点。

属性:
- `index`: 唯一标识符
- `name`: 母线名称
- `vn_kv`: 额定电压(kV)
- `type`: 母线类型(P, DC_REF)
- `area`: 所属区域
- `zone`: 所属分区
- `in_service`: 运行状态
- `vm_pu`: 电压幅值(标幺值)

### Line
```julia
mutable struct Line <: AbstractComponent
    index::Int
    name::String
    from_bus::Int
    to_bus::Int
    length_km::Float64
    r_ohm_per_km::Float64
    x_ohm_per_km::Float64
    c_nf_per_km::Float64
    r0_ohm_per_km::Float64
    x0_ohm_per_km::Float64
    c0_nf_per_km::Float64
    g_us_per_km::Float64
    max_i_ka::Float64
    type::String
    max_loading_percent::Float64
    parallel::Int
    df::Float64
    in_service::Bool
    mtbf_hours::Float64
    mttr_hours::Float64
    sw_hours::Float64
    rp_hours::Float64
end
```
表示交流输电线路或配电线路。

属性:
- `index`: 唯一标识符
- `name`: 线路名称
- `from_bus`: 起始母线ID
- `to_bus`: 终止母线ID
- `length_km`: 线路长度(km)
- `r_ohm_per_km`: 正序电阻(Ω/km)
- `x_ohm_per_km`: 正序电抗(Ω/km)
- `c_nf_per_km`: 正序电容(nF/km)
- `r0_ohm_per_km`: 零序电阻(Ω/km)
- `x0_ohm_per_km`: 零序电抗(Ω/km)
- `c0_nf_per_km`: 零序电容(nF/km)
- `g_us_per_km`: 电导(μS/km)
- `max_i_ka`: 最大电流容量(kA)
- `type`: 线路类型(cs-电缆, ol-架空线)
- `max_loading_percent`: 最大负载百分比
- `parallel`: 并联线路数量
- `df`: 分布因子
- `in_service`: 运行状态
- `mtbf_hours`: 平均故障间隔时间(小时)
- `mttr_hours`: 平均修复时间(小时)
- `sw_hours`: 开关操作时间(小时)
- `rp_hours`: 修复准备时间(小时)

### LineDC
```julia
mutable struct LineDC <: AbstractComponent
    index::Int
    name::String
    from_bus::Int
    to_bus::Int
    length_km::Float64
    r_ohm_per_km::Float64
    max_i_ka::Float64
    max_loading_percent::Float64
    parallel::Int
    in_service::Bool
    # 其他属性...
end
```
表示直流输电线路或配电线路。

属性:
- `index`: 唯一标识符
- `name`: 线路名称
- `from_bus`: 起始母线ID
- `to_bus`: 终止母线ID
- `length_km`: 线路长度(km)
- `r_ohm_per_km`: 电阻(Ω/km)
- `max_i_ka`: 最大电流容量(kA)
- `max_loading_percent`: 最大负载百分比
- `parallel`: 并联线路数量
- `in_service`: 运行状态

## 发电组件

### StaticGenerator
```julia
mutable struct StaticGenerator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    k::Float64
    rx::Float64
    in_service::Bool
    type::String
    controllable::Bool
end
```
表示静态发电机，如风电、光伏等。

属性:
- `index`: 唯一标识符
- `name`: 发电机名称
- `bus`: 连接母线ID
- `p_mw`: 有功功率输出(MW)
- `q_mvar`: 无功功率输出(MVar)
- `scaling`: 功率缩放因子
- `max_p_mw`: 最大有功功率(MW)
- `min_p_mw`: 最小有功功率(MW)
- `max_q_mvar`: 最大无功功率(MVar)
- `min_q_mvar`: 最小无功功率(MVar)
- `k`: 无功功率系数
- `rx`: 电抗/电阻比
- `in_service`: 运行状态
- `type`: 发电机类型(WP-风电, PV-光伏, CHP-热电联产等)
- `controllable`: 是否可控

### Generator
```julia
mutable struct Generator <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    vm_pu::Float64
    sn_mva::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    vn_kv::Float64
    xdss_pu::Float64
    rdss_pu::Float64
    cos_phi::Float64
    controllable::Bool
    in_service::Bool
    type::String
    generator_type::String
    fuel_type::String
    startup_time_cold_h::Float64
    startup_time_warm_h::Float64
    startup_time_hot_h::Float64
    min_up_time_h::Float64
    min_down_time_h::Float64
    ramp_up_rate_mw_per_min::Float64
    ramp_down_rate_mw_per_min::Float64
    # 其他属性...
end
```
表示常规发电机组。

属性:
- `index`: 唯一标识符
- `name`: 发电机名称
- `bus`: 连接母线ID
- `p_mw`: 有功功率输出(MW)
- `vm_pu`: 电压控制设定值(标幺值)
- `sn_mva`: 额定视在功率(MVA)
- `scaling`: 功率缩放因子
- `max_p_mw`: 最大有功功率(MW)
- `min_p_mw`: 最小有功功率(MW)
- `max_q_mvar`: 最大无功功率(MVar)
- `min_q_mvar`: 最小无功功率(MVar)
- `vn_kv`: 额定电压(kV)
- `xdss_pu`: 次暂态电抗(标幺值)
- `rdss_pu`: 次暂态电阻(标幺值)
- `cos_phi`: 功率因数
- `controllable`: 是否可控
- `in_service`: 运行状态
- `type`: 发电机类型
- `generator_type`: 发电机技术类型
- `fuel_type`: 燃料类型
- `startup_time_cold_h`: 冷启动时间(小时)
- `startup_time_warm_h`: 温启动时间(小时)
- `startup_time_hot_h`: 热启动时间(小时)
- `min_up_time_h`: 最小开机时间(小时)
- `min_down_time_h`: 最小关机时间(小时)
- `ramp_up_rate_mw_per_min`: 爬坡上升速率(MW/分钟)
- `ramp_down_rate_mw_per_min`: 爬坡下降速率(MW/分钟)

### StaticGeneratorDC
```julia
mutable struct StaticGeneratorDC <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    scaling::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    in_service::Bool
    type::String
    controllable::Bool
    # 其他属性...
end
```
表示直流系统中的静态发电机。

属性:
- `index`: 唯一标识符
- `name`: 发电机名称
- `bus`: 连接母线ID
- `p_mw`: 功率输出(MW)
- `scaling`: 功率缩放因子
- `max_p_mw`: 最大功率(MW)
- `min_p_mw`: 最小功率(MW)
- `in_service`: 运行状态
- `type`: 发电机类型
- `controllable`: 是否可控

## 负荷组件

### Load
```julia
mutable struct Load <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    q_mvar::Float64
    const_z_percent::Float64
    const_i_percent::Float64
    const_p_percent::Float64
    scaling::Float64
    in_service::Bool
    # 其他属性...
end
```
表示交流系统中的对称负荷。

属性:
- `index`: 唯一标识符
- `name`: 负荷名称
- `bus`: 连接母线ID
- `p_mw`: 有功功率(MW)
- `q_mvar`: 无功功率(MVar)
- `const_z_percent`: 恒阻抗负荷百分比
- `const_i_percent`: 恒电流负荷百分比
- `const_p_percent`: 恒功率负荷百分比
- `scaling`: 负荷缩放因子
- `in_service`: 运行状态

### LoadDC
```julia
mutable struct LoadDC <: AbstractComponent
    index::Int
    name::String
    bus::Int
    p_mw::Float64
    const_z_percent::Float64
    const_i_percent::Float64
    const_p_percent::Float64
    scaling::Float64
    in_service::Bool
    # 其他属性...
end
```
表示直流系统中的负荷。

属性:
- `index`: 唯一标识符
- `name`: 负荷名称
- `bus`: 连接母线ID
- `p_mw`: 功率(MW)
- `const_z_percent`: 恒阻抗负荷百分比
- `const_i_percent`: 恒电流负荷百分比
- `const_p_percent`: 恒功率负荷百分比
- `scaling`: 负荷缩放因子
- `in_service`: 运行状态

## 储能组件

### Storage
```julia
mutable struct Storage <: AbstractComponent
    index::Int
    name::String
    bus::Int
    power_capacity_mw::Float64
    energy_capacity_mwh::Float64
    soc_init::Float64
    min_soc::Float64
    max_soc::Float64
    efficiency::Float64
    in_service::Bool
    type::String
    controllable::Bool
    # 其他属性...
end
```
表示储能系统。

属性:
- `index`: 唯一标识符
- `name`: 储能系统名称
- `bus`: 连接母线ID
- `power_capacity_mw`: 功率容量(MW)
- `energy_capacity_mwh`: 能量容量(MWh)
- `soc_init`: 初始荷电状态
- `min_soc`: 最小荷电状态
- `max_soc`: 最大荷电状态
- `efficiency`: 充放电效率
- `in_service`: 运行状态
- `type`: 储能类型
- `controllable`: 是否可控

## 转换器组件

### Converter
```julia
mutable struct Converter <: AbstractComponent
    index::Int
    name::String
    bus_ac::Int
    bus_dc::Int
    p_mw::Float64
    q_mvar::Float64
    vm_ac_pu::Float64
    vm_dc_pu::Float64
    loss_percent::Float64
    max_p_mw::Float64
    min_p_mw::Float64
    max_q_mvar::Float64
    min_q_mvar::Float64
    in_service::Bool
    type::String
    mode::String
    droop_kp::Float64
    # 其他属性...
end
```
表示AC/DC转换器。

属性:
- `index`: 唯一标识符
- `name`: 转换器名称
- `bus_ac`: 交流侧母线ID
- `bus_dc`: 直流侧母线ID
- `p_mw`: 有功功率(MW)
- `q_mvar`: 无功功率(MVar)
- `vm_ac_pu`: 交流侧电压(标幺值)
- `vm_dc_pu`: 直流侧电压(标幺值)
- `loss_percent`: 损耗百分比
- `max_p_mw`: 最大有功功率(MW)
- `min_p_mw`: 最小有功功率(MW)
- `max_q_mvar`: 最大无功功率(MVar)
- `min_q_mvar`: 最小无功功率(MVar)
- `in_service`: 运行状态
- `type`: 转换器类型
- `mode`: 控制模式
- `droop_kp`: 下垂控制系数

## 电动汽车组件

### ChargingStation
```julia
mutable struct ChargingStation <: AbstractComponent
    index::Int
    name::String
    bus::Int
    location::String
    operator::String
    num_chargers::Int
    max_power_kw::Float64
    in_service::Bool
end
```
表示电动汽车充电站。

属性:
- `index`: 唯一标识符
- `name`: 充电站名称
- `bus`: 连接母线ID
- `location`: 地理位置
- `operator`: 运营商
- `num_chargers`: 充电桩数量
- `max_power_kw`: 最大功率(kW)
- `in_service`: 运行状态

### Charger
```julia
mutable struct Charger <: AbstractComponent
    index::Int
    name::String
    station_id::Int
    type::String
    power_kw::Float64
    voltage_v::Float64
    current_a::Float64
    efficiency::Float64
    in_service::Bool
    # 其他属性...
end
```
表示充电桩。

属性:
- `index`: 唯一标识符
- `name`: 充电桩名称
- `station_id`: 所属充电站ID
- `type`: 充电桩类型
- `power_kw`: 功率(kW)
- `voltage_v`: 电压(V)
- `current_a`: 电流(A)
- `efficiency`: 效率
- `in_service`: 运行状态

## 微电网组件

### Microgrid
```julia
mutable struct Microgrid <: AbstractComponent
    index::Int
    name::String
    description::String
    control_area::String
    capacity_mw::Float64
    energy_mwh::Float64
    response_time_s::Float64
    ramp_rate_mw_per_min::Float64
    availability_percent::Float64
    operator::String
    in_service::Bool
    # 其他属性...
end
```
表示微电网系统。

属性:
- `index`: 唯一标识符
- `name`: 微电网名称
- `description`: 描述
- `control_area`: 控制区域
- `capacity_mw`: 容量(MW)
- `energy_mwh`: 能量(MWh)
- `response_time_s`: 响应时间(秒)
- `ramp_rate_mw_per_min`: 爬坡速率(MW/分钟)
- `availability_percent`: 可用性百分比
- `operator`: 运营商
- `in_service`: 运行状态

## 虚拟电厂组件

### VirtualPowerPlant
```julia
mutable struct VirtualPowerPlant <: AbstractComponent
    index::Int
    name::String
    description::String
    control_area::String
    capacity_mw::Float64
    energy_mwh::Float64
    response_time_s::Float64
    ramp_rate_mw_per_min::Float64
    availability_percent::Float64
    operator::String
    in_service::Bool
    resource_type::String
    resource_id::Int
    capacity_share_percent::Float64
    control_priority::Int
    resource_response_time_s::Float64
    max_duration_h::Float64
    timestamp::DateTime
    p_mw::Float64
    q_mvar::Float64
    flexibility_up_mw::Float64
    flexibility_down_mw::Float64
    flexibility_duration_h::Float64
    # 其他属性...
end
```
表示虚拟电厂。

属性:
- `index`: 唯一标识符
- `name`: 虚拟电厂名称
- `description`: 描述
- `control_area`: 控制区域
- `capacity_mw`: 容量(MW)
- `energy_mwh`: 能量(MWh)
- `response_time_s`: 响应时间(秒)
- `ramp_rate_mw_per_min`: 爬坡速率(MW/分钟)
- `availability_percent`: 可用性百分比
- `operator`: 运营商
- `in_service`: 运行状态
- `resource_type`: 资源类型
- `resource_id`: 资源ID
- `capacity_share_percent`: 容量分担百分比
- `control_priority`: 控制优先级
- `resource_response_time_s`: 资源响应时间(秒)
- `max_duration_h`: 最大持续时间(小时)
- `timestamp`: 时间戳
- `p_mw`: 有功功率(MW)
- `q_mvar`: 无功功率(MVar)
- `flexibility_up_mw`: 上调灵活性(MW)
- `flexibility_down_mw`: 下调灵活性(MW)
- `flexibility_duration_h`: 灵活性持续时间(小时)

## 碳排放组件

### CarbonTimeSeries
```julia
mutable struct CarbonTimeSeries <: AbstractComponent
    index::Int
    timestamp::DateTime
    grid_carbon_intensity_kgCO2e_per_MWh::Float64
    renewable_generation_carbon_intensity_kgCO2e_per_MWh::Float64
    storage_carbon_intensity_kgCO2e_per_MWh::Float64
    # 其他属性...
end
```
表示碳排放时间序列。

属性:
- `index`: 唯一标识符
- `timestamp`: 时间戳
- `grid_carbon_intensity_kgCO2e_per_MWh`: 电网碳强度(kgCO2e/MWh)
- `renewable_generation_carbon_intensity_kgCO2e_per_MWh`: 可再生能源碳强度(kgCO2e/MWh)
- `storage_carbon_intensity_kgCO2e_per_MWh`: 储能碳强度(kgCO2e/MWh)

### EquipmentCarbon
```julia
mutable struct EquipmentCarbon <: AbstractComponent
    index::Int
    element_type::String
    element_id::Int
    carbon_intensity_kgCO2e_per_MWh::Float64
    lifecycle_emissions_kgCO2e::Float64
    operational_emissions_kgCO2e_per_year::Float64
    # 其他属性...
end
```
表示设备碳排放。

属性:
- `index`: 唯一标识符
- `element_type`: 设备类型
- `element_id`: 设备ID
- `carbon_intensity_kgCO2e_per_MWh`: 碳强度(kgCO2e/MWh)
- `lifecycle_emissions_kgCO2e`: 生命周期排放(kgCO2e)
- `operational_emissions_kgCO2e_per_year`: 运行排放(kgCO2e/年)

## 使用示例

### 创建电力系统模型

```julia
using ComponentModel

# 创建母线
bus1 = Bus(1, "Bus 1", 110.0, "PQ", "Area1", "Zone1", true, 1.0, 0.0)
bus2 = Bus(2, "Bus 2", 110.0, "PV", "Area1", "Zone1", true, 1.02, 0.0)
bus3 = Bus(3, "Bus 3", 110.0, "REF", "Area1", "Zone1", true, 1.05, 0.0)

# 创建线路
line1 = Line(1, "Line 1-2", 1, 2, 10.0, 0.1, 0.4, 10.0, 0.3, 1.2, 8.0, 0.0, 0.5, "ol", 80.0, 1, 1.0, true, 8760.0, 4.0, 1.0, 2.0)
line2 = Line(2, "Line 2-3", 2, 3, 15.0, 0.1, 0.4, 10.0, 0.3, 1.2, 8.0, 0.0, 0.5, "ol", 80.0, 1, 1.0, true, 8760.0, 4.0, 1.0, 2.0)

# 创建发电机
gen1 = Generator(1, "Generator 1", 3, 100.0, 1.05, 120.0, 1.0, 150.0, 50.0, 100.0, -100.0, 110.0, 0.2, 0.01, 0.85, true, true, "Thermal", "Steam", "Coal", 4.0, 2.0, 1.0, 4.0, 2.0, 5.0, 5.0)

# 创建负荷
load1 = Load(1, "Load 1", 1, 80.0, 30.0, 30.0, 30.0, 40.0, 1.0, true)
```

### 分析和计算

```julia
# 假设有一个计算潮流的函数
function run_power_flow(buses, lines, generators, loads)
    # 潮流计算实现
    # ...
    return results
end

# 使用组件进行潮流计算
results = run_power_flow([bus1, bus2, bus3], [line1, line2], [gen1], [load1])
```

### 数据导入导出

```julia
# 从CSV文件导入母线数据
function import_buses_from_csv(filename)
    buses = Bus[]
    # 读取CSV文件并创建Bus对象
    # ...
    return buses
end

# 导出结果到CSV文件
function export_results_to_csv(results, filename)
    # 将结果写入CSV文件
    # ...
end

# 使用示例
buses = import_buses_from_csv("buses.csv")
export_results_to_csv(results, "results.csv")
```

## 依赖项

- `Dates`: 用于处理时间戳
- `DataFrames`: 用于数据处理和分析
- `CSV`: 用于文件导入导出
- `Plots`: 用于结果可视化

## 版本历史

- v0.1.0: 初始版本，基本组件定义
- v0.2.0: 添加直流系统组件
- v0.3.0: 添加微电网和虚拟电厂组件
- v0.4.0: 添加碳排放相关组件
- v0.5.0: 性能优化和bug修复

## 贡献指南

欢迎通过以下方式为ComponentModel做出贡献:

1. 提交问题报告和功能请求
2. 提交代码修复和新功能
3. 改进文档和示例

请确保所有提交都包含适当的测试和文档。

## 许可证

ComponentModel使用MIT许可证。详情请参见LICENSE文件。