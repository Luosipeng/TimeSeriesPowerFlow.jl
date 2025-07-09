using LinearAlgebra
using SparseArrays
using IncompleteLU

function GMRES_economic(A, x_0, e, b)
    j = 0
    r = b - A * x_0
    beta = norm(r)
    V = zeros(length(r), length(r))
    V[:, 1] = r / norm(r)
    Eps = beta * Matrix{Float64}(I, size(x_0, 1), 1)
    H = zeros(size(A))  # Initialize H here
    c = zeros(size(A, 1))  # Initialize c here
    s = zeros(size(A, 1))  # Initialize s here
    m = 0  # Initialize m here
    residual_error = 0  # Initialize residual_error
    while true
        j += 1
        w = zeros(length(r), j)
        w[:, j] = A * V[:, j]

        for i = 1:j
            H[i, j] = dot(w[:, j], V[:, i])
            w[:, j] = w[:, j] - H[i, j] * V[:, i]
        end

        H[j+1, j] = norm(w[:, j])

        if H[j+1, j] == 0
            m = j
            break
        else
            V[:, j+1] = w[:, j] / H[j+1, j]
        end

        if abs(H[j, j]) > abs(H[j+1, j])
            tau = H[j+1, j] / H[j, j]
            c[j] = 1 / sqrt(1 + tau^2)
            s[j] = c[j] * tau
        else
            tau = H[j, j] / H[j+1, j]
            s[j] = 1 / sqrt(1 + tau^2)
            c[j] = s[j] * tau
        end

        H[j, j] = c[j] * H[j, j] + s[j] * H[j+1, j]
        H[j+1, j] = 0

        Eps[j:j+1] = [c[j] s[j]; -s[j] c[j]] * [Eps[j]; 0]

        residual_error = abs(Eps[j+1]) * beta
        if abs(Eps[j+1]) * beta < e
            m = j
            break
        end
        end

        y = H[1:m, 1:m] \ Eps[1:m]
        x = x_0 + V[:, 1:length(y)] * y  # Use length(y) instead of size(y)

    return x, residual_error, m
end
function GMRES_ILU(A, x_0, e, b)
    # Compute the ILU preconditioner
    lu = ilu(A, τ = 0.1)  # τ is the drop tolerance

    j = 0
    r = lu \ (b - A * x_0)  # Apply the preconditioner here
    beta = norm(r)
    V = zeros(length(r), length(r))
    V[:, 1] = r / norm(r)
    Eps = beta * Matrix{Float64}(I, size(x_0, 1), 1)
    H = zeros(size(A))  # Initialize H here
    c = zeros(size(A, 1))  # Initialize c here
    s = zeros(size(A, 1))  # Initialize s here
    m = 0  # Initialize m here
    residual_error = 0  # Initialize residual_error
    while true
        j += 1
        w = zeros(length(r), j)
        w[:, j] = lu \ (A * V[:, j])  # Apply the preconditioner here

        for i = 1:j
            H[i, j] = dot(w[:, j], V[:, i])
            w[:, j] = w[:, j] - H[i, j] * V[:, i]
        end

        H[j+1, j] = norm(w[:, j])

        if H[j+1, j] == 0
            m = j
            break
        else
            V[:, j+1] = w[:, j] / H[j+1, j]
        end

        if abs(H[j, j]) > abs(H[j+1, j])
            tau = H[j+1, j] / H[j, j]
            c[j] = 1 / sqrt(1 + tau^2)
            s[j] = c[j] * tau
        else
            tau = H[j, j] / H[j+1, j]
            s[j] = 1 / sqrt(1 + tau^2)
            c[j] = s[j] * tau
        end

        H[j, j] = c[j] * H[j, j] + s[j] * H[j+1, j]
        H[j+1, j] = 0

        Eps[j:j+1] = [c[j] s[j]; -s[j] c[j]] * [Eps[j]; 0]

        residual_error = abs(Eps[j+1]) * beta
        if abs(Eps[j+1]) * beta < e
            m = j
            break
        end
    end

    y = H[1:m, 1:m] \ Eps[1:m]
    x = x_0 + V[:, 1:length(y)] * y  # Use length(y) instead of size(y)

    return x, residual_error, m
end
# n = 1000
# A = sprand(n, n, 1.0) 
# A = A + I*n
# b = sparse(rand(n))
# x_0=zeros(length(b))
# maxit=1000
# e=1e-6
# @time x, residual_error, m = GMRES_economic(A, x_0, e, b)
# @time x, residual_error, m ==GMRES_ILU(A, x_0, e, b)
# #Reasons for non-convergence
