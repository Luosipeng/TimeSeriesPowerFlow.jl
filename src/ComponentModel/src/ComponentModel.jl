"""
ComponentModel

This module defines the component model structure for microgrid systems,
including AC/DC components, hybrid elements, and carbon emission models.
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

# Export all public interfaces
export AbstractComponent

# Basic components
export Bus, Line, LineDC, Switch, ThreePhaseBranch, BusDC, HighVoltageCircuitBreaker

# Transformer components
export Transformer2W, Transformer3W, Transformer2Wetap

# Generation components
export StaticGenerator, Generator, StaticGeneratorDC

# Load components
export Load, AsymmetricLoad, LoadDC

# Storage components
export Storage, MobileStorage, Storageetap, StorageAC

# Grid components
export ExternalGrid, Converter

# Carbon emission components
export EquipmentCarbon, CarbonTimeSeries, CarbonScenario

# Virtual power plant components
export VirtualPowerPlant, FlexLoad

# Electric vehicle components
export ChargingStation, Charger, EVAggregator, V2GService

# Microgrid components
export Microgrid

# Photovoltaic components
export PVArray, ACPVSystem

end # module ComponentModel
