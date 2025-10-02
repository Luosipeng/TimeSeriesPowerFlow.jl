# API Reference

```@meta
CurrentModule = HyDistFlow
```

## AC Components
```@autodocs
Modules = [ComponentModel]
Pages   = ["Nodes.jl","Generators.jl","Loads.jl","Lines.jl","Transformers.jl"]
Order   = [:type, :function, :macro, :constant]
```

## DC Components
```@autodocs
Modules = [ComponentModel]
Pages   = ["DCGenerators.jl","DCLines.jl","DCLoads.jl","DCNodes.jl","DCStorages.jl"]
Order   = [:type, :function, :macro, :constant]
Filter  = t -> !(t in [StaticGeneratorDC, LineDC, LoadDC, PVArray, BusDC])
```

## Hybrid Components
```@autodocs
Modules = [ComponentModel]
Pages   = ["Converters.jl","EVCharging.jl","Microgrids.jl","VirtualPlant.jl"]
Order   = [:type, :function, :macro, :constant]
```

## Carbon Components
```@autodocs
Modules = [ComponentModel]
Pages   = ["CarbonScenarios.jl","EmissionModels.jl"]
Order   = [:type, :function, :macro, :constant]
```