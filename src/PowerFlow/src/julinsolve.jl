"""
    julinsolve(A, b, solver = "", opt = nothing)

Solve a linear system Ax = b using various solver methods.

# Arguments
- `A`: Coefficient matrix
- `b`: Right-hand side vector
- `solver`: String specifying the solution method (default: `""` uses Julia's backslash operator)
- `opt`: Optional parameters for solvers (not currently used)

# Returns
- `x`: Solution vector
- `info`: Additional information about the solution process (currently always `nothing`)

# Available Solvers
- `""` or `"\\"`: Julia's default backslash operator
- `"LU3"`: AMD ordering with sparse or dense LU factorization
- `"LU3a"`: AMD ordering with LU factorization (L, U, p)
- `"LU4"`: LU factorization with row and column permutations (L, U, p, q)
- `"LU5"`: LU factorization with row and column permutations and row scaling (L, U, p, q, R)
- `"cholesky"`: Cholesky factorization for symmetric positive definite matrices
- `"gmres"`: Generalized Minimal Residual method with ILU preconditioning
- `"bicgstab"`: Biconjugate Gradient Stabilized method with ILU preconditioning
- `"cgs"`: Conjugate Gradient Squared method
- `"gpuLU"`: LU factorization using CUDA for GPU acceleration

# Description
This function provides a unified interface to various linear system solvers in Julia.
It supports direct methods (like LU and Cholesky factorizations) and iterative methods
(like GMRES, BiCGSTAB, and CGS). For sparse systems, appropriate ordering and 
preconditioning techniques are applied to improve performance.

The GPU-based solver requires CUDA support and appropriate hardware.

# Examples
```julia
A = rand(100, 100)
b = rand(100)
x, _ = julinsolve(A, b)                  # Default solver
x, _ = julinsolve(A, b, "cholesky")      # Cholesky factorization
x, _ = julinsolve(A, b, "gmres")         # GMRES iterative method
"""

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
