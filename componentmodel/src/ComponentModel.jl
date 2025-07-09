"""
微电网系统数据结构定义
"""
module ComponentModel

using SparseArrays
using DataFrames
using Dates

abstract type AbstractComponent end

# Include all the necessary modules
include("../AC/Generators.jl")
include("../AC/Lines.jl")
include("../AC/Loads.jl")
include("../AC/Nodes.jl")
include("../AC/Transformers.jl")
include("../DC/DCLines.jl")
include("../DC/DCStorages.jl")
include("../DC/DCNodes.jl")
include("../DC/DCGenerators.jl")
include("../DC/DCLoads.jl")
include("../Hybrid/Converters.jl")
include("../Hybrid/EVCharging.jl")
include("../Hybrid/VirtualPlant.jl")
include("../Hybrid/Microgrids.jl")
include("../Carbon/CarbonScenarios.jl")
include("../Carbon/EmissionModels.jl")

# 导出所有公共接口
export AbstractComponent

# 基础组件
export Bus, Line, LineDC, Switch, ThreePhaseBranch, BusDC, HighVoltageCircuitBreaker

# 变压器组件
export Transformer2W, Transformer3W,Transformer2Wetap

# 发电组件
export StaticGenerator, Generator, StaticGeneratorDC

# 负荷组件
export Load, AsymmetricLoad, LoadDC

# 储能组件
export Storage, MobileStorage, Storageetap, StorageAC

# 电网组件
export ExternalGrid, Converter

# 碳排放组件
export EquipmentCarbon, CarbonTimeSeries, CarbonScenario

# 虚拟电厂组件
export VirtualPowerPlant, FlexLoad

# 电动汽车组件
export ChargingStation, Charger, EVAggregator, V2GService

# 微电网组件
export Microgrid

# 光伏组件
export PVArray, ACPVSystem

end # module MicrogridSystem
