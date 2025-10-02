# API Reference

```@meta
CurrentModule = HyDistFlow
```

## Power Flow Algorithm
```@autodocs
Modules = [PowerFlow]
Pages   = ["adaptive_damped_newton.jl","currentinjectionpf.jl","newtonpf.jl","newtonpf_gpu.jl"]
Order   = [:type, :function, :macro, :constant]
```

## DC Power Flow Algorithm
```@autodocs
Modules = [PowerFlow]
Pages   = ["newtondcpf.jl","newtondcpf_sp.jl"]
Order   = [:type, :function, :macro, :constant]
```

## Linear Problem Solution Algorithm
```@autodocs
Modules = [PowerFlow]
Pages   = ["julinsolve.jl"]
Order   = [:type, :function, :macro, :constant]
```

## GPU Acceleration
```@autodocs
Modules = [PowerFlow]
Pages   = ["makeSbus_gpu.jl","makeSdzip_gpu.jl"]
Order   = [:type, :function, :macro, :constant]
```

## Other Functions
```@autodocs
Modules = [PowerFlow]
Pages   = ["makeSbus.jl","makeSdzip.jl","makeYbus.jl","merge_results.jl","pf_summary.jl","pfsoln.jl","process_result.jl","result_compare_etap.jl","rundcpf.jl","runhpf.jl","runpf.jl","runupf.jl","total_load.jl"]
Order   = [:type, :function, :macro, :constant]
Filter = f -> !(f in [PowerFlow.newtondcpf])
```
