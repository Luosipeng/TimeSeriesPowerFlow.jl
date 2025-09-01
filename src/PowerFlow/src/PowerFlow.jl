"""
    Define the power flow module with different functions
"""

module PowerFlow
    using Printf
    using SparseArrays
    using LinearAlgebra
    using PrettyTables
    using AMD
    using SuiteSparse
    using IterativeSolvers
    using IncompleteLU
    using KrylovKit
    using Krylov
    using KrylovPreconditioners
    using LinearOperators
    using CUDA
    using CUDA.CUSPARSE
    # using AlgebraicMultigrid
    using XLSX
    using DataFrames
    using StatsBase
    using Graphs
    using DataStructures  
    using Dates
    using Base.Threads
    using Printf 
    # using CUDA, CUDA.CUSPARSE
    # using CUDSS
    using Plots
    using Test
    # using KrylovPreconditioners
    # using different packages based on the operating system
    # ...other files in component_models directory...
    
    # ...other files in Utils directory...
    include("../../Utils/src/Utils.jl")
    using .Utils
    using .Utils: JuliaPowerCase2Jpc,JuliaPowerCase2Jpc_3ph
    const ComponentModel = Utils.ComponentModel

    # ... other files directly under src ...

    include(joinpath(@__DIR__,"bustypes.jl"))
    # include(joinpath(@__DIR__,"ext2int.jl"))
    include(joinpath(@__DIR__,"makeYbus.jl"))
    include(joinpath(@__DIR__,"newtonpf.jl"))
    include(joinpath(@__DIR__,"makeSbus.jl"))
    include(joinpath(@__DIR__,"makeSdzip.jl"))
    include(joinpath(@__DIR__,"julinsolve.jl"))
    include(joinpath(@__DIR__,"pfsoln.jl"))
    include(joinpath(@__DIR__,"dSbus_dV.jl"))
    include(joinpath(@__DIR__,"runpf.jl"))
    include(joinpath(@__DIR__,"settings.jl"))
    include(joinpath(@__DIR__,"rundcpf.jl"))
    include(joinpath(@__DIR__,"makeBdc.jl"))
    include(joinpath(@__DIR__,"dcpf.jl"))
    include(joinpath(@__DIR__,"total_load.jl"))
    include(joinpath(@__DIR__,"runprepf.jl"))
    include(joinpath(@__DIR__,"dc_preprocess.jl"))
    include(joinpath(@__DIR__,"build_branch.jl"))
    include(joinpath(@__DIR__,"build_bus.jl"))
    include(joinpath(@__DIR__,"build_gen.jl"))
    include(joinpath(@__DIR__,"build_load.jl"))
    include(joinpath(@__DIR__,"process_pv_acsystem.jl"))
    include(joinpath(@__DIR__,"adaptive_damped_newton.jl"))
    include(joinpath(@__DIR__,"currentinjectionpf.jl"))
    include(joinpath(@__DIR__,"newtondcpf_sp.jl"))
    include(joinpath(@__DIR__,"eliminate_element.jl"))
    # include(joinpath(@__DIR__,"run_single_day.jl"))
    # include(joinpath(@__DIR__,"run_dynamic_dispatch.jl"))

    # include(joinpath(@__DIR__,"extract_data.jl"))

    include(joinpath(@__DIR__,"runhpf.jl"))
    include(joinpath(@__DIR__,"dcbustypes.jl"))
    include(joinpath(@__DIR__,"newtondcpf.jl"))
    include(joinpath(@__DIR__,"dcpfsoln.jl"))

    # include(joinpath(@__DIR__,"model_examine.jl"))
    include(joinpath(@__DIR__,"pf_summary.jl"))
    include(joinpath(@__DIR__,"merge_results.jl"))
    include(joinpath(@__DIR__,"process_result.jl"))
    include(joinpath(@__DIR__,"result_compare_etap.jl"))
    include(joinpath(@__DIR__,"makeSbus_gpu.jl"))
    include(joinpath(@__DIR__,"newtonpf_gpu.jl"))
    include(joinpath(@__DIR__,"makeSdzip_gpu.jl"))
    include(joinpath(@__DIR__,"newtondcpf_gpu.jl"))
    include(joinpath(@__DIR__,"runupf.jl"))


    # ... other files in models directory ...

    # ... other files in test directory ...
    # include(joinpath(dirname(@__DIR__), "test","lindistflow_inverter_evaluate.jl"))
    # include(joinpath(dirname(@__DIR__), "test","ac_element_validate.jl"))
    # include(joinpath(dirname(@__DIR__), "test","dc_element_validate.jl"))
    # include(joinpath(dirname(@__DIR__), "test","acdc_power_flow_compared.jl"))
    # include(joinpath(dirname(@__DIR__), "test","ac_power_flow_compared.jl"))
    # include(joinpath(dirname(@__DIR__), "test","modelingexamine.jl"))
    # include(joinpath(dirname(@__DIR__), "test","resultcompare.jl"))
    # include(joinpath(dirname(@__DIR__), "test", "loadflow_result_ETAP.jl"))
    # include(joinpath(dirname(@__DIR__), "test","acdccompare.jl"))
    # # ... other files in ios directory ...
    # include(joinpath(dirname(@__DIR__), "ios","matlab2julia.jl"))
    # include(joinpath(dirname(@__DIR__), "ios","ETAPImporter.jl"))


    
    
    export idx_bus, idx_brch, idx_gen, bustypes, makeYbus, newtonpf, makeSbus, makeSdzip, mplinsolve, total_load, pfsoln, dSbus_dV, MPC
    export runpf, rundcpf, makeBdc, dcpf, int2ext, runprepf, dc_preprocess, build_branch, build_bus
    export build_gen, options, runhpf_iteration,runpf,runhpf, runprepf, runupf,rundcpf
    export runhpf, dcbustypes, newtondcpf, dcpfsoln
    export  pf_summary, merge_results, process_result
end
