using CUDA, CUDA.CUSPARSE
using CUDSS
using SparseArrays, LinearAlgebra
using BenchmarkTools

function compare_solvers(n, density=0.005)
    # 创建相同的稀疏矩阵和向量（只创建一次）
    T = Float64
    A_cpu = sprand(T, n, n, density) + I
    b_cpu = rand(T, n)
    
    println("矩阵大小: $n x $n, 非零元素数: $(nnz(A_cpu)), 密度: $(round(density*100, digits=2))%")
    
    # ===== CPU求解（明确使用LU分解）=====
    cpu_time = @elapsed begin
        # 明确使用LU分解而不是默认的\运算符
        F = lu(A_cpu)
        x_cpu = F \ b_cpu
    end
    
    # 计算CPU残差
    r_cpu = b_cpu - A_cpu * x_cpu
    cpu_residual = norm(r_cpu)
    
    # ===== GPU求解（使用LU分解）=====
    # 转移到GPU
    A_gpu = CuSparseMatrixCSR(A_cpu)
    x_gpu = CuVector(zeros(T, n))
    b_gpu = CuVector(b_cpu)
    
    gpu_time = @elapsed begin
        # 创建求解器，使用'F'参数指定LU分解
        solver = CudssSolver(A_gpu, "G", 'F')  # 'F'表示LU分解
        
        # 求解过程
        cudss("analysis", solver, x_gpu, b_gpu)
        cudss("factorization", solver, x_gpu, b_gpu)
        cudss("solve", solver, x_gpu, b_gpu)
    end
    
    # 计算GPU残差
    r_gpu = b_gpu - A_gpu * x_gpu
    gpu_residual = norm(r_gpu)
    
    # 计算解的差异
    x_gpu_cpu = Array(x_gpu)  # 将GPU解转回CPU
    solution_diff = norm(x_cpu - x_gpu_cpu) / norm(x_cpu)
    
    # 输出结果
    println("  CPU时间(LU分解): $(round(cpu_time, digits=6))秒")
    println("  GPU时间(LU分解): $(round(gpu_time, digits=6))秒")
    println("  加速比: $(round(cpu_time/gpu_time, digits=2))倍")
    println("  CPU残差: $(cpu_residual)")
    println("  GPU残差: $(gpu_residual)")
    println("  解的相对差异: $(solution_diff)")
    println()
    
    # 添加多次求解的测试
    multi_result = nothing
    if n <= 5000  # 对较大矩阵限制多次求解测试，避免内存问题
        multi_result = test_multiple_solves(A_cpu, n)
    end
    
    return (n=n, cpu_time=cpu_time, gpu_time=gpu_time, 
            speedup=cpu_time/gpu_time, cpu_residual=cpu_residual, 
            gpu_residual=gpu_residual, solution_diff=solution_diff,
            multi_result=multi_result)
end

# 测试多次求解的情况（重用分解）
function test_multiple_solves(A_cpu, n, num_solves=10)
    T = Float64
    A_gpu = CuSparseMatrixCSR(A_cpu)
    
    # GPU多次求解（只进行一次分析和分解）
    solver = CudssSolver(A_gpu, "G", 'F')
    x_gpu = CuVector(zeros(T, n))
    b_gpu = CuVector(rand(T, n))
    
    # 初始化
    cudss("analysis", solver, x_gpu, b_gpu)
    cudss("factorization", solver, x_gpu, b_gpu)
    
    # 多次求解不同的b
    gpu_multi_time = @elapsed begin
        for i in 1:num_solves
            b_gpu = CuVector(rand(T, n))
            cudss("solve", solver, x_gpu, b_gpu)
        end
    end
    
    # CPU多次求解（重用LU分解）
    F = lu(A_cpu)  # 只分解一次
    
    cpu_multi_time = @elapsed begin
        for i in 1:num_solves
            b_cpu = rand(T, n)
            x_cpu = F \ b_cpu  # 重用LU分解
        end
    end
    
    println("  多次求解($(num_solves)次):")  # 修正这里的错误
    println("    CPU总时间(重用LU): $(round(cpu_multi_time, digits=6))秒")
    println("    GPU总时间(重用LU): $(round(gpu_multi_time, digits=6))秒")
    println("    平均每次加速比: $(round(cpu_multi_time/gpu_multi_time, digits=2))倍")
    
    return (cpu_multi_time=cpu_multi_time, gpu_multi_time=gpu_multi_time, 
            multi_speedup=cpu_multi_time/gpu_multi_time)
end

# 测试不同大小的矩阵
function run_comparison()
    matrix_sizes = [100, 500, 1000, 2000, 5000, 10000]
    results = []
    
    for n in matrix_sizes
        try
            result = compare_solvers(n)
            push!(results, result)
        catch e
            println("矩阵大小 $n 测试失败: $e")
        end
    end
    
    # 打印汇总表格
    println("\n性能比较汇总:")
    println("| 矩阵大小 | CPU时间(秒) | GPU时间(秒) | 加速比 | CPU残差 | GPU残差 | 解的相对差异 |")
    println("|----------|------------|------------|--------|---------|---------|--------------|")
    
    for r in results
        println("| $(r.n) | $(round(r.cpu_time, digits=6)) | $(round(r.gpu_time, digits=6)) | $(round(r.speedup, digits=2)) | $(r.cpu_residual) | $(r.gpu_residual) | $(r.solution_diff) |")
    end
    
    return results
end

# 运行比较
results = run_comparison()
