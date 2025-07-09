function newtondcpf_sp(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    
    # 初始化
    converged = false
    i = 0
    V = V0
    
    # 设置索引
    np = length(p)
    
    # 计算初始功率不平衡
    # 直流系统中，功率 P = V * I，其中 I = G * V (G是电导矩阵)
    Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
    # 直流系统只关心有功功率
    P_calc = real.(V .* (Ybus * V))
    P_spec = real.(Sbus)
    
    # 计算功率不平衡
    F = P_spec[p] - P_calc[p]
    
    # 检查收敛性
    normF = norm(F, Inf)
    if normF < tol
        converged = true
    end
    
    # 牛顿迭代
    while (!converged && i < max_it)
        # 更新迭代计数
        i += 1
        
        # 构建雅可比矩阵
        # 对于直流系统，J = ∂P/∂V
        J = zeros(np, np)
        for j in 1:np
            for k in 1:np
                if j == k
                    # 对角元素
                    J[j,j] = 2 * real(Ybus[p[j],p[j]]) * V[p[j]] + 
                             sum(real(Ybus[p[j],m]) * V[m] for m in 1:length(V) if m != p[j])
                else
                    # 非对角元素
                    J[j,k] = real(Ybus[p[j],p[k]]) * V[p[j]]
                end
            end
        end
        
        # 计算更新步长
        dV = J \ F
        
        # 更新电压
        V[p] += dV
        
        # 重新计算功率不平衡
        Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        P_calc = real.(V .* (Ybus * V))
        P_spec = real.(Sbus)
        F = P_spec[p] - P_calc[p]
        
        # 检查收敛性
        normF = norm(F, Inf)
        println("Iteration $i: normF = $normF")
        if normF < tol
            converged = true
        end
    end
    
    return V, converged, i
end
