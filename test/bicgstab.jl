using LinearAlgebra
using IncompleteLU
using SparseArrays
function bicgstab(A, b, maxiter=1000, tol=1e-6)
    n = size(A, 1)
    x = zeros(n)
    r = b - A * x
    rhat = r
    rho = alpha = omega = 1.0
    v = p = zeros(n)
    
    # Preconditioner M
    # M = ilu(A, τ=0.1)  # τ is the drop tolerance
    D = Diagonal(A)
    M_inv = inv(D)

    for iter = 1:maxiter
        rho_1 = rho
        rho = dot(rhat, r)
        beta = (rho/rho_1) * (alpha/omega)
        p = r + beta * (p - omega * v)
        
        # Apply preconditioner
        # phat = similar(p)
        # ldiv!(phat, M, p)
        phat = M_inv * p
        v = A * phat
        
        alpha = rho / dot(rhat, v)
        s = r - alpha * v
        
        # Apply preconditioner
        # shat = similar(s)
        # ldiv!(shat, M, s)
        shat = M_inv * s
        t = A * shat
        
        omega = dot(t, s) / dot(t, t)
        x = x + alpha * phat + omega * shat
        r = s - omega * t
        
        if norm(r) < tol
            return x
        end
    end
    
    error("BiCGSTAB did not converge within the maximum number of iterations")
end