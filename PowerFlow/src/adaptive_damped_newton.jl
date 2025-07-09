function adaptive_damped_newton(baseMVA, bus, gen, load, pvarray, Ybus, V0, ref, p, tol0, max_it0, alg="")
    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    
    # Initialize
    converged = false
    i = 0
    V = V0
    
    # Set up indexing for updating V
    np = length(p)
    j1 = 1; j2 = np; # j1:j2 - V angle of pv buses
    
    # Evaluate F(x0)
    mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
    F = real(mis[p])
    
    # Check tolerance
    normF = norm(F, Inf)
    if normF < tol
        converged = true
    end
    
    # Do Newton iterations
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1

        # Evaluate Jacobian
        dSbus_dVa, dSbus_dVm = PowerFlow.dSbus_dV(Ybus, V)
        # neg_dSd_dVm = PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray, return_derivative=true)
        # dSbus_dVm .-= neg_dSd_dVm

        J = real(dSbus_dVm[p,p])

        # Compute update step
        dx, info = PowerFlow.julinsolve(J, -F, alg)
        # dx = J \ -F
        
        # 创建电压更新的副本，用于测试不同阻尼因子
        V_temp = copy(V)
        
        # 初始阻尼因子
        damping = 1.0
        
        # 应用初始更新并计算新的不平衡量
        if np > 0
            V_temp[p] .+= damping * dx[j1:j2]
        end
        
        mis_temp = V_temp .* conj.(Ybus * V_temp) - PowerFlow.makeSbus(baseMVA, bus, gen, V_temp, load, pvarray)
        F_temp = real(mis_temp[p])
        new_normF = norm(F_temp, Inf)
        
        # 自适应调整阻尼因子
        # 如果新的不平衡量更大，则减小阻尼因子
        min_damping = 0.1  # 最小阻尼因子限制
        max_attempts = 5   # 最大尝试次数
        attempt = 0
        
        while new_normF > normF && damping > min_damping && attempt < max_attempts
            # 减小阻尼因子
            damping *= 0.5
            attempt += 1
            
            # 使用新的阻尼因子重新计算
            V_temp = copy(V)
            if np > 0
                V_temp[p] .+= damping * dx[j1:j2]
            end
            
            mis_temp = V_temp .* conj.(Ybus * V_temp) - PowerFlow.makeSbus(baseMVA, bus, gen, V_temp, load, pvarray)
            F_temp = real(mis_temp[p])
            new_normF = norm(F_temp, Inf)
        end
        
        # 应用最终确定的阻尼因子更新电压
        if np > 0
            V[p] .+= damping * dx[j1:j2]
        end
        
        # 评估新状态下的不平衡量
        mis = V .* conj.(Ybus * V) - PowerFlow.makeSbus(baseMVA, bus, gen, V, load, pvarray)
        F = real(mis[p])

        # 检查收敛性
        normF = norm(F, Inf)
        
        # 可选：输出迭代信息，包括阻尼因子
        println("Iteration $i: normF = $normF, damping = $damping")
        λ_min, v_min, info = KrylovKit.eigsolve(J, 1, :SR, tol=1e-10)
        println("最小特征值: ", λ_min[1])
        
        if normF < tol
            converged = true
        end
    end

    return V, converged, i
end
