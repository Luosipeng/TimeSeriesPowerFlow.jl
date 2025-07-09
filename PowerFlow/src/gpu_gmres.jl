using CUDA
using IncompleteLU
using LinearAlgebra
using SparseArrays
using Krylov
using Plots

#ILU preconditioner
function gpu_gmres(A,x0,tol,b)
    # Compute the ILU preconditioner
    lu = ilu(A, τ = 0.1)  # τ is the drop tolerance
    j = 0
    r = lu \ (b - A * x0)  # Apply the preconditioner here
    beta = norm(r)
    V = zeros(length(r), length(r))
    V[:, 1] = r / norm(r)
    Eps = beta * CUDA.ones(Float32, size(x0, 1), 1)
    H = CUDA.zeros(size(A))  # Initialize H here
    c = CUDA.zeros(size(A, 1))  # Initialize c here
    s = CUDA.zeros(size(A, 1))  # Initialize s here
    m = 0  # Initialize m here
    residual_error = 0  # Initialize residual_error
    while true
        j += 1
        w = zeros(length(r), j)
        w[:, j] = lu \ (A * V[:, j])  # Apply the preconditioner here
        w=CuArray(w)
        V=CuArray(V)
        for i = 1:j
            # Use broadcasting to assign the result to H[i, j]
            H[i:i, j:j] = CUDA.dot(w[:, j], V[:, i])
        
            # Perform the subtraction using broadcasting
            w[:, j:j] .= w[:, j:j] .-  H[i:i, j:j] .* V[:, i:i]
        end

        H[j+1:j+1, j:j] = norm(w[:, j])

        if H[j+1:j+1, j:j] == 0
            m = j
            break
        else
            V[:, j+1] = w[:, j] ./ H[j+1:j+1, j:j]
        end

        if all(abs.(H[j:j, j:j]) .> abs.(H[j+1:j+1, j:j]))
            tau = H[j+1:j+1, j:j] / H[j:j, j:j]
            c[j:j] = 1 ./ sqrt(1 .+ tau.^2)
            s[j:j] = c[j:j] * tau
        else
            tau = H[j:j, j:j] / H[j+1:j+1, j:j]
            s[j:j] = 1 ./ sqrt(1 .+ tau.^2)
            c[j:j] = s[j:j] * tau
        end

        H[j:j, j:j] = c[j:j] * H[j:j, j:j] + s[j:j] * H[j+1:j+1, j:j]
        H[j+1:j+1, j:j] = 0

        Eps[j:j+1] = [c[j:j] s[j:j]; -s[j:j] c[j:j]] * [Eps[j:j]; 0]

        residual_error = abs.(Eps[j+1:j+1]) * beta
        if  all(abs.(Eps[j+1:j+1]) * beta .< tol)
            m = j
            break
        end
    end

    y = H[1:m, 1:m] \ Eps[1:m]
    x = CuArray(x0) + V[:, 1:length(y)] * y  # Use length(y) instead of size(y)

    return Array(x), Array(residual_error), m
end

#jacobi preconditioner
function gpu_gmres_jacobi(A,x0,tol,b)
     # Compute the Jacobi preconditioner
     M_inv = 1 ./ diag(A)
     j = 0
     r = M_inv .* (b - A * x0)  # Apply the preconditioner here
     beta = norm(r)
     V = zeros(length(r), length(r))
     V[:, 1] = r / norm(r)
     Eps = beta * CUDA.ones(Float32, size(x0, 1), 1)
     H = CUDA.zeros(size(A))  # Initialize H here
     c = CUDA.zeros(size(A, 1))  # Initialize c here
     s = CUDA.zeros(size(A, 1))  # Initialize s here
     m = 0  # Initialize m here
     residual_error = 0  # Initialize residual_error
     while true
         j += 1
         w = zeros(length(r), j)
         w[:, j] = M_inv .* Array(A * Array(V[:, j]))  # Apply the preconditioner here
         w=CuArray(w)
         V=CuArray(V)
         for i = 1:j
             # Use broadcasting to assign the result to H[i, j]
             H[i:i, j:j] = CUDA.dot(w[:, j], V[:, i])
         
             # Perform the subtraction using broadcasting
             w[:, j:j] .= w[:, j:j] .-  H[i:i, j:j] .* V[:, i:i]
         end
 
         H[j+1:j+1, j:j] = norm(w[:, j])
 
         if H[j+1:j+1, j:j] == 0
             m = j
             break
         else
             V[:, j+1] = w[:, j] ./ H[j+1:j+1, j:j]
         end
 
         if all(abs.(H[j:j, j:j]) .> abs.(H[j+1:j+1, j:j]))
             tau = H[j+1:j+1, j:j] / H[j:j, j:j]
             c[j:j] = 1 ./ sqrt(1 .+ tau.^2)
             s[j:j] = c[j:j] * tau
         else
             tau = H[j:j, j:j] / H[j+1:j+1, j:j]
             s[j:j] = 1 ./ sqrt(1 .+ tau.^2)
             c[j:j] = s[j:j] * tau
         end
 
         H[j:j, j:j] = c[j:j] * H[j:j, j:j] + s[j:j] * H[j+1:j+1, j:j]
         H[j+1:j+1, j:j] = 0
 
         Eps[j:j+1] = [c[j:j] s[j:j]; -s[j:j] c[j:j]] * [Eps[j:j]; 0]
 
         residual_error = abs.(Eps[j+1:j+1]) * beta
         if  all(abs.(Eps[j+1:j+1]) * beta .< tol)
             m = j
             break
         end
     end
 
     y = H[1:m, 1:m] \ Eps[1:m]
     x = CuArray(x0) + V[:, 1:length(y)] * y  # Use length(y) instead of size(y)
 
     return Array(x), Array(residual_error), m
 end
# n = 1000
# A = sprand(n, n, 1.0) 
# #A = A + I*n
# b = sparse(rand(n))
# x0=zeros(length(b))
# maxit=1000
# tol=1e-6
# @time x1,res,m=gpu_gmres_jacobi(A,x0,tol,b)
# @time x2 , stats = Krylov.gmres(A, b)
# mis=x1-x2
# plot(mis)