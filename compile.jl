using PackageCompiler

create_app("demo",  # 你的应用目录
                 "PowerFlowApp",         # 输出目录
           precompile_execution_file="C:/Users/DELL/Desktop/DistributionPowerFlow/demo/precompile_app.jl",
           filter_stdlibs=false,
           include_lazy_artifacts=true,  # 添加这个参数
           force=true
       )