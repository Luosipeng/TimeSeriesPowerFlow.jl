using CUDA, CUDA.CUSPARSE
using CUDSS
using SparseArrays, LinearAlgebra

T = Float64
n = 100
A_cpu = sprand(T, n, n, 0.05) + I
x_cpu = zeros(T, n)
b_cpu = rand(T, n)

A_gpu = CuSparseMatrixCSR(A_cpu)
x_gpu = CuVector(x_cpu)
b_gpu = CuVector(b_cpu)

solver = CudssSolver(A_gpu, "G", 'F')

cudss("analysis", solver, x_gpu, b_gpu)
cudss("factorization", solver, x_gpu, b_gpu)
cudss("solve", solver, x_gpu, b_gpu)

r_gpu = b_gpu - A_gpu * x_gpu
norm(r_gpu)

# In-place LU
d_gpu = rand(T, n) |> CuVector
A_gpu = A_gpu + Diagonal(d_gpu)
cudss_set(solver, A_gpu)

c_cpu = rand(T, n)
c_gpu = CuVector(c_cpu)

cudss("refactorization", solver, x_gpu, c_gpu)
cudss("solve", solver, x_gpu, c_gpu)

r_gpu = c_gpu - A_gpu * x_gpu
norm(r_gpu)