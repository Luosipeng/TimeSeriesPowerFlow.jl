# todo list:
# 1. test the gpu based algorithms
# 2. test the cgs solvers

function julinsolve(A, b, solver = "", opt = nothing)
    info = nothing

    if solver in ["", "\\"]
        x = A \ b
    elseif solver == "LU3"
        q = amd(A)
        if issparse(A)
            # 使用与示例代码相同的方式进行LU分解
            F = lu(A[q,q])
            x = zeros(size(A, 1))
            x[q] = F \ b[q]  # 直接使用分解结果求解
        else
            L, U, p = lu(A[q,q])
            x = zeros(size(A, 1))
            x[q] = U \ (L \ b[q[p]])
        end
    
    elseif solver == "LU3a"
        q = amd(A)
        L, U, p = lu(A[q,q])
        x = zeros(size(A, 1))
        x[q] = U \ (L \ b[q[p]])
    elseif solver == "LU4"
        L, U, p, q = lu(A)
        x = zeros(size(A, 1))
        x[q] = U \ (L \ b[p])
    elseif solver == "LU5"
        L, U, p, q, R = lu(A)
        x = zeros(size(A, 1))
        x[q] = U \ (L \ (R[:, p] \ b))
    elseif solver == "cholesky"
        factor = cholesky(A)
        x = factor \ b
    elseif solver == "gmres"
        ilu_fact = ilu(A)
        x = IterativeSolvers.gmres(A, b, Pl=ilu_fact, reltol=1e-8, maxiter = 1000)
    elseif solver == "bicgstab"
        n = size(A,1)
        F = ilu(A, τ = 0.05)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> forward_substitution!(y, F, v))
        opN = LinearOperator(Float64, n, n, false, false, (y, v) -> backward_substitution!(y, F, v))
        x , stats = bicgstab(A, b, history=false, M=opM, N=opN)
    elseif solver == "cgs"
        x = cgs(A, b, rtol=1e-8, itmax=1000)
    elseif solver == "gpuLU"
        T = Float64
        n= size(A,1)
        solver = CudssSolver(A, "G", 'F')
        x = CUDA.zeros(T, n)
        cudss("analysis", solver, x, b)
        cudss("factorization", solver, x, b)
        cudss("solve", solver, x, b)
    else
        error("mplinsolve: '$solver' is not a valid value for SOLVER, using default.")
        x = A \ b
    end

    return x, info
end
