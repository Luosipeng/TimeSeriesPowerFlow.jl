push!(LOAD_PATH, pwd()*"/src/")
using PowerFlow
using BenchmarkTools
using Statistics
using Printf
using MATLAB
using Base.Threads

# 设置MATLAB路径，指向MATPOWER数据目录
mat"addpath('C:/Users/13733/Desktop/matpower-8.0/data')"
mat"addpath('C:/Users/13733/Desktop/matpower-8.0')"

# 定义测试用的案例列表
test_cases = [
    "case14",
    "case30",
    "case300",
    "case1888rte",
    "case6468rte",
    "case13659pegase"  # 已添加case13659pegase
]

# 定义要测试的求解器
solvers = [
    "LU3", 
    "gpuLU"
]

# 创建结果存储结构
results = Dict()
for case_name in test_cases
    results[case_name] = Dict()
    for solver in solvers
        results[case_name][solver] = Dict(
            "time" => Float64[], 
            "success" => Bool[],
            "Vm" => [],
            "Va" => []
        )
    end
end

# 设置基准测试参数
benchmark_samples = 10

# 运行基准测试
for case_name in test_cases
    println("Testing case: $case_name")
    
    for solver in solvers
        println("  Solver: $solver")
        
        # 设置选项
        opt = PowerFlow.options()
        opt["PF"]["ENFORCE_Q_LIMS"] = 0
        
        if solver == "gpuLU"
            opt["PF"]["NR_ALG"] = "gpuLU"
            opt["PF"]["GPU_ACCELERATION"] = 1
            opt["PF"]["DC_PREPROCESS"] = 1
        else
            opt["PF"]["NR_ALG"] = solver
            opt["PF"]["GPU_ACCELERATION"] = 0
            opt["PF"]["DC_PREPROCESS"] = 0
        end
        
        # 运行基准测试
        for i = 1:benchmark_samples
            try
                # 关键修改：每次测试都重新加载MATPOWER案例
                # 使用您确认有效的方式加载案例
                mpc = nothing  # 先清空之前的对象
                if case_name == "case14"
                    mpc = mat"case14"
                elseif case_name == "case30"
                    mpc = mat"case30"
                elseif case_name == "case300"
                    mpc = mat"case300"
                elseif case_name == "case1888rte"
                    mpc = mat"case1888rte"
                elseif case_name == "case6468rte"
                    mpc = mat"case6468rte"
                elseif case_name == "case13659pegase"
                    mpc = mat"case13659pegase"  # 添加case13659pegase的加载
                end
                
                # 测量运行时间
                t = @elapsed begin
                    result = PowerFlow.runpf(mpc, opt)
                    success = result["success"]
                end
                
                push!(results[case_name][solver]["time"], t)
                push!(results[case_name][solver]["success"], success)
                
                # 存储电压幅值和相角结果
                if i == 1 && success
                    results[case_name][solver]["Vm"] = result["bus"][:,8]
                    results[case_name][solver]["Va"] = result["bus"][:,9]
                end
                
                print(".")
            catch e
                push!(results[case_name][solver]["time"], NaN)
                push!(results[case_name][solver]["success"], false)
                println("\nError with $solver on $case_name: $e")
                println("Error details: ", sprint(showerror, e, catch_backtrace()))
            end
        end
        println()
    end
end

# 计算并显示性能结果
println("\n===== PERFORMANCE BENCHMARK RESULTS =====")
println("Case\tSolver\tAvg Time (ms)\tSuccess Rate\tSpeedup vs LU3")

reference_times = Dict()

for case_name in test_cases
    # 先计算LU3的时间作为参考
    if haskey(results[case_name], "LU3")
        lu3_times = filter(!isnan, results[case_name]["LU3"]["time"])
        if length(lu3_times) > 0
            reference_times[case_name] = mean(lu3_times) * 1000  # 转换为毫秒
        end
    end
    
    for solver in solvers
        times = filter(!isnan, results[case_name][solver]["time"])
        success_rate = mean(results[case_name][solver]["success"]) * 100
        
        if length(times) > 0
            avg_time = mean(times) * 1000  # 转换为毫秒
            
            # 计算相对于LU3的加速比
            speedup = haskey(reference_times, case_name) ? reference_times[case_name] / avg_time : NaN
            
            @printf("%s\t%s\t%.2f\t\t%.1f%%\t\t%.2fx\n", 
                    case_name, solver, avg_time, success_rate, speedup)
        else
            println("$case_name\t$solver\tFailed\t\t0.0%\t\tN/A")
        end
    end
    println()
end

# 添加结果比较部分 - 比较LU3和gpuLU的解的差异
println("\n===== SOLUTION COMPARISON =====")
println("Case\tVm Diff (max)\tVa Diff (max)\tVm Diff (avg)\tVa Diff (avg)")

for case_name in test_cases
    if haskey(results[case_name]["LU3"], "Vm") && haskey(results[case_name]["gpuLU"], "Vm")
        # 确保两个求解器都有有效结果
        if !isempty(results[case_name]["LU3"]["Vm"]) && !isempty(results[case_name]["gpuLU"]["Vm"])
            # 计算电压幅值和相角的差异
            vm_diff = abs.(results[case_name]["LU3"]["Vm"] - results[case_name]["gpuLU"]["Vm"])
            va_diff = abs.(results[case_name]["LU3"]["Va"] - results[case_name]["gpuLU"]["Va"])
            
            max_vm_diff = maximum(vm_diff)
            max_va_diff = maximum(va_diff)
            avg_vm_diff = mean(vm_diff)
            avg_va_diff = mean(va_diff)
            
            @printf("%s\t%.2e\t%.2e\t%.2e\t%.2e\n", 
                    case_name, max_vm_diff, max_va_diff, avg_vm_diff, avg_va_diff)
        else
            println("$case_name\tNo valid solution comparison available")
        end
    else
        println("$case_name\tNo valid solution comparison available")
    end
end
