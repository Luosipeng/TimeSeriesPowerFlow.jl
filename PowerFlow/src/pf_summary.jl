"""
写入系统摘要部分
"""
function write_system_summary(f::IOStream, mpc::JPC, area, isolated)
    # 提取必要的数据
    baseMVA = mpc.baseMVA
    
    # 计算基本统计数据
    if hasproperty(mpc, :busAC)
        n_buses = size(mpc.busAC, 1)
    else
        # 从其他数据推断母线数量
        bus_ids = Set()
        if hasproperty(mpc, :branchAC)
            for i in 1:size(mpc.branchAC, 1)
                push!(bus_ids, mpc.branchAC[i, 1])
                push!(bus_ids, mpc.branchAC[i, 2])
            end
        end
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                push!(bus_ids, mpc.genAC[i, 1])
            end
        end
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                push!(bus_ids, mpc.loadAC[i, 1])
            end
        end
        n_buses = length(bus_ids)
    end
    n_isolated = length(isolated)
    n_buses = n_buses + n_isolated
    
    # 获取发电机数量
    n_gens = hasproperty(mpc, :genAC) ? size(mpc.genAC, 1) : 0
    
    # 获取负荷数量
    n_loads = hasproperty(mpc, :loadAC) ? size(mpc.loadAC, 1) : 0
    
    # 获取支路数量
    n_branches = hasproperty(mpc, :branchAC) ? size(mpc.branchAC, 1) : 0
    
    # 计算变压器数量 (假设branch矩阵中有tap列)
    n_transformers = 0
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 9
        for i in 1:size(mpc.branchAC, 1)
            if abs(mpc.branchAC[i, 9]) > 1e-6
                n_transformers += 1
            end
        end
    end
    
    # 计算发电总量和负荷总量
    total_gen_p = 0.0
    total_gen_q = 0.0
    if hasproperty(mpc, :genAC)
        # 假设gen矩阵的第2列是有功功率，第3列是无功功率
        total_gen_p = sum(mpc.genAC[:, 2]) 
        total_gen_q = sum(mpc.genAC[:, 3]) 
    end
    
    total_load_p = 0.0
    total_load_q = 0.0
    if hasproperty(mpc, :busAC)
        # 假设busAC矩阵的第3列是有功负荷，第4列是无功负荷
        if size(mpc.busAC, 2) >= 4
            total_load_p = sum(mpc.busAC[:, 3]) 
            total_load_q = sum(mpc.busAC[:, 4]) 
        end
    elseif hasproperty(mpc, :loadAC)
        # 从loadAC中获取负荷数据
        if size(mpc.loadAC, 2) >= 4
            total_load_p = sum(mpc.loadAC[:, 3])
            total_load_q = sum(mpc.loadAC[:, 4])
        end
    end
    
    # 计算损耗
    total_p_loss = total_gen_p - total_load_p
    total_q_loss = 0.0
    charging_q = 0.0
    
    if hasproperty(mpc, :branchAC)
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
                
            # 获取线路参数
            r = mpc.branchAC[i, 3] 
            x = mpc.branchAC[i, 4]
            b = mpc.branchAC[i, 5]  # 线路半充电电纳
                
            # 获取实际的母线电压值和相角
            v_from = 1.0  # 默认值，如果找不到实际电压
            v_to = 1.0    # 默认值，如果找不到实际电压
            ang_from = 0.0
            ang_to = 0.0
                
            # 从母线数据中查找实际电压值和相角
            if hasproperty(mpc, :busAC)
                for j in 1:size(mpc.busAC, 1)
                    if Int(mpc.busAC[j, 1]) == from_bus
                        v_from = mpc.busAC[j, 8]  # 使用实际电压幅值
                        ang_from = mpc.busAC[j, 9] * pi/180  # 转换为弧度
                    elseif Int(mpc.busAC[j, 1]) == to_bus
                        v_to = mpc.busAC[j, 8]    # 使用实际电压幅值
                        ang_to = mpc.busAC[j, 9] * pi/180  # 转换为弧度
                    end
                end
            end
                
            # 相角差
            angle_diff = ang_from - ang_to
                
            # 计算线路导纳
            y = 1 / complex(r, x)
            y_abs = abs(y)
                
            # 直接计算线路电流幅值平方
            i_mag_squared = (v_from^2 + v_to^2 - 2*v_from*v_to*cos(angle_diff)) * y_abs^2
                
            # 计算无功损耗 - 使用电抗和电流平方
            q_loss = x * i_mag_squared * baseMVA

            charging_from = 0.5 * b * v_from^2 * baseMVA
            charging_to = 0.5 * b * v_to^2 * baseMVA
            charging_q += charging_from + charging_to
            total_q_loss += q_loss
        end
    end
    
    # 计算发电机容量
    total_gen_pmax = 0.0
    total_gen_qmin = 0.0
    total_gen_qmax = 0.0
    if hasproperty(mpc, :genAC) && size(mpc.genAC, 2) >= 9
        total_gen_pmax = sum(mpc.genAC[:, 9]) 
        total_gen_qmin = sum(mpc.genAC[:, 5])
        total_gen_qmax = sum(mpc.genAC[:, 4])
    end
    
    # 写入系统摘要
    write(f, "================================================================================\n")
    write(f, "|     System Summary                                                           |\n")
    write(f, "================================================================================\n\n")
    
    @printf(f, "How many?                How much?              P (MW)            Q (MVAr)\n")
    @printf(f, "---------------------    -------------------  -------------  -----------------\n")
    @printf(f, "Buses            %3d     Total Gen Capacity    %7.1f       %7.1f to %7.1f\n", 
            n_buses, total_gen_pmax, total_gen_qmin, total_gen_qmax)
    @printf(f, "Generators       %3d     On-line Capacity      %7.1f       %7.1f to %7.1f\n", 
            n_gens, total_gen_pmax, total_gen_qmin, total_gen_qmax)
    @printf(f, "Committed Gens   %3d     Generation (actual)   %7.1f             %7.1f\n", 
            n_gens, total_gen_p, total_gen_q)
    @printf(f, "Loads            %3d     Load                  %7.1f            %7.1f\n", 
            n_loads, total_load_p, total_load_q)
    @printf(f, "  Fixed          %3d       Fixed               %7.1f            %7.1f\n", 
            n_loads, total_load_p, total_load_q)
    @printf(f, "  Dispatchable    %2d       Dispatchable          %4.1f of %4.1f      %5.1f\n", 
            0, 0.0, 0.0, 0.0)
    @printf(f, "Shunts           %3d     Shunt (inj)             %5.1f              %5.1f\n", 
            0, 0.0, 0.0)  # 假设没有分路元件
    @printf(f, "Branches         %3d     Losses (I^2 * Z)       %6.2f            %6.2f\n", 
            n_branches, total_p_loss, total_q_loss)  # 无功损耗需要计算
    @printf(f, "Transformers     %3d     Branch Charging (inj)     -             %6.1f\n", 
            n_transformers, charging_q)  # 需要计算实际值
    @printf(f, "Inter-ties        %2d     Total Inter-tie Flow     %4.1f               %4.1f\n", 
            0, 0.0, 0.0)
    @printf(f, "Areas             %2d\n\n", area)
    
    # 电压和相角的最大最小值
    # 这部分需要从bus数据中提取，如果没有完整的bus数据，可以省略或使用估计值
    if hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 9
        min_vm = Inf
        max_vm = -Inf
        min_vm_bus = 0
        max_vm_bus = 0
        min_va = Inf
        max_va = -Inf
        min_va_bus = 0
        max_va_bus = 0
        
        for i in 1:size(mpc.busAC, 1)
            vm = mpc.busAC[i, 8]
            va = mpc.busAC[i, 9]
            
            if vm < min_vm
                min_vm = vm
                min_vm_bus = Int(mpc.busAC[i, 1])
            end
            if vm > max_vm
                max_vm = vm
                max_vm_bus = Int(mpc.busAC[i, 1])
            end
            
            if va < min_va
                min_va = va
                min_va_bus = Int(mpc.busAC[i, 1])
            end
            if va > max_va
                max_va = va
                max_va_bus = Int(mpc.busAC[i, 1])
            end
        end
        
        @printf(f, "                          Minimum                      Maximum\n")
        @printf(f, "                 -------------------------  --------------------------------\n")
        @printf(f, "Voltage Magnitude   %5.3f p.u. @ bus %3d         %5.3f p.u. @ bus %3d  \n", 
                min_vm, min_vm_bus, max_vm, max_vm_bus)
        @printf(f, "Voltage Angle      %6.2f deg   @ bus %3d        %6.2f deg   @ bus %3d  \n", 
                min_va, min_va_bus, max_va, max_va_bus)
    end
    
    # 线路损耗信息
    # 如果有详细的支路损耗数据，可以添加这部分
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 18
        # 处理支路损耗数据
        max_p_loss= -Inf
        max_q_loss = -Inf
        max_p_line = 0
        max_q_line = 0
        p_loss = mpc.branchAC[:, 15] + mpc.branchAC[:, 17]
        q_loss = mpc.branchAC[:, 16] + mpc.branchAC[:, 18]
        for i in 1:size(mpc.branchAC, 1)
            ploss= p_loss[i]
            qloss= q_loss[i]
            if ploss > max_p_loss
                max_p_loss = ploss
                max_p_line = i
            end
            if qloss > max_q_loss
                max_q_loss = qloss
                max_q_line = i
            end
        end
        # 使用估计值或省略
        @printf(f, "P Losses (I^2*R)             -                  %5.2f MW    @ line %s\n", 
        max_p_loss, string(Int(mpc.branchAC[max_p_line, 1]), "-", Int(mpc.branchAC[max_p_line, 2])))
        @printf(f, "Q Losses (I^2*X)             -                 %5.2f MVAr  @ line %s\n", 
        max_q_loss, string(Int(mpc.branchAC[max_q_line, 1]), "-", Int(mpc.branchAC[max_q_line, 2])))
    else
        # 使用估计值或省略
        @printf(f, "P Losses (I^2*R)             -                  %5.2f MW    @ line %s\n", 
                0.0, "X-X")
        @printf(f, "Q Losses (I^2*X)             -                 %5.2f MVAr  @ line %s\n", 
                0.0, "X-X")
    end
    
    write(f, "\n")
end

"""
写入母线数据部分
"""
function write_bus_data(f::IOStream, mpc::JPC, isolated)
    baseMVA = mpc.baseMVA
    
    write(f, "================================================================================\n")
    write(f, "|     Bus Data                                                                 |\n")
    write(f, "================================================================================\n")
    write(f, " Bus      Voltage          Generation             Load        \n")
    write(f, "  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n")
    write(f, "----- ------- --------  --------  --------  --------  --------\n")
    
    # 假设我们有完整的bus数据
    if hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 9
        # 创建发电机和负荷查找表
        gen_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                bus_id = Int(mpc.genAC[i, 1])
                pg = mpc.genAC[i, 2] 
                qg = mpc.genAC[i, 3] 
                gen_lookup[bus_id] = (pg, qg)
            end
        end
        
        load_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                bus_id = Int(mpc.loadAC[i, 1])
                pd = mpc.loadAC[i, 3] 
                qd = mpc.loadAC[i, 4] 
                load_lookup[bus_id] = (pd, qd)
            end
        elseif hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 4
            for i in 1:size(mpc.busAC, 1)
                bus_id = Int(mpc.busAC[i, 1])
                pd = mpc.busAC[i, 3] 
                qd = mpc.busAC[i, 4] 
                if pd != 0.0 || qd != 0.0
                    load_lookup[bus_id] = (pd, qd)
                end
            end
        end
        
        # 添加孤岛
        bus_data = copy(mpc.busAC)
        for i in eachindex(isolated)
            isolated_bus = zeros(1, size(bus_data, 2))
            isolated_bus[1] = isolated[i]
            isolated_bus[2] = 1.0  # 类型（PQ节点）
            isolated_bus[8] = 0.0  # 默认电压幅值
            isolated_bus[9] = 0.0  # 默认相角
            bus_data = vcat(bus_data, isolated_bus)
        end
        
        # 根据母线编号对bus_data进行排序
        bus_ids = bus_data[:, 1]
        sorted_indices = sortperm(bus_ids)
        bus_data = bus_data[sorted_indices, :]

        # 总计
        total_pg = 0.0
        total_qg = 0.0
        total_pd = 0.0
        total_qd = 0.0
        
        # 遍历所有母线
        for i in 1:size(bus_data, 1)
            bus_id = Int(bus_data[i, 1])
            vm = bus_data[i, 8]
            va = bus_data[i, 9]
            
            # 获取发电数据
            pg_str = "-"
            qg_str = "-"
            if haskey(gen_lookup, bus_id)
                pg, qg = gen_lookup[bus_id]
                pg_str = @sprintf("%.2f", pg)
                qg_str = @sprintf("%.2f", qg)
                total_pg += pg
                total_qg += qg
            end
            
            # 获取负荷数据
            pd = 0.0
            qd = 0.0
            if haskey(load_lookup, bus_id)
                pd, qd = load_lookup[bus_id]
                total_pd += pd
                total_qd += qd
            end
            
            # 打印母线数据
            @printf(f, "%5d  %5.3f   %6.3f   %8s   %8s   %7.2f   %7.2f \n", 
                    bus_id, vm, va, pg_str, qg_str, pd, qd)
        end
        
        # 打印总计
        @printf(f, "                        --------  --------  --------  --------\n")
        @printf(f, "               Total:   %7.2f   %7.2f   %7.2f   %7.2f\n", 
                total_pg, total_qg, total_pd, total_qd)
    else
        # 如果没有完整的bus数据，尝试从其他数据构建
        bus_ids = Set{Int}()
        
        # 从gen、load和branch数据中收集所有母线ID
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                push!(bus_ids, Int(mpc.genAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                push!(bus_ids, Int(mpc.loadAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :branchAC)
            for i in 1:size(mpc.branchAC, 1)
                push!(bus_ids, Int(mpc.branchAC[i, 1]))
                push!(bus_ids, Int(mpc.branchAC[i, 2]))
            end
        end
        
        # 创建发电机和负荷查找表
        gen_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                bus_id = Int(mpc.genAC[i, 1])
                pg = mpc.genAC[i, 2] * baseMVA
                qg = mpc.genAC[i, 3] * baseMVA
                gen_lookup[bus_id] = (pg, qg)
            end
        end
        
        load_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                bus_id = Int(mpc.loadAC[i, 1])
                pd = mpc.loadAC[i, 3] * baseMVA
                qd = mpc.loadAC[i, 4] * baseMVA
                load_lookup[bus_id] = (pd, qd)
            end
        end
        
        # 总计
        total_pg = 0.0
        total_qg = 0.0
        total_pd = 0.0
        total_qd = 0.0
        
        # 遍历所有收集到的母线ID
        for bus_id in sort(collect(bus_ids))
            # 假设电压数据
            vm = 1.0  # 默认值
            va = 0.0  # 默认值
            
            # 获取发电数据
            pg_str = "-"
            qg_str = "-"
            if haskey(gen_lookup, bus_id)
                pg, qg = gen_lookup[bus_id]
                pg_str = @sprintf("%.2f", pg)
                qg_str = @sprintf("%.2f", qg)
                total_pg += pg
                total_qg += qg
            end
            
            # 获取负荷数据
            pd = 0.0
            qd = 0.0
            if haskey(load_lookup, bus_id)
                pd, qd = load_lookup[bus_id]
                total_pd += pd
                total_qd += qd
            end
            
            # 打印母线数据
            @printf(f, "%5d  %5.3f   %6.3f   %8s   %8s   %7.2f   %7.2f \n", 
                    bus_id, vm, va, pg_str, qg_str, pd, qd)
        end
        
        # 打印总计
        @printf(f, "                        --------  --------  --------  --------\n")
        @printf(f, "               Total:   %7.2f   %7.2f   %7.2f   %7.2f\n", 
                total_pg, total_qg, total_pd, total_qd)
    end
    
    write(f, "\n")
end

"""
写入支路数据部分
"""
function write_branch_data(f::IOStream, mpc::JPC)
    baseMVA = mpc.baseMVA
    
    write(f, "================================================================================\n")
    write(f, "|     Branch Data                                                              |\n")
    write(f, "================================================================================\n")
    write(f, "Brnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  \n")
    write(f, "  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n")
    write(f, "-----  -----  -----  --------  --------  --------  --------  --------  --------\n")
    
    # 检查是否有支路数据
    if !hasproperty(mpc, :branchAC) || size(mpc.branchAC, 1) == 0
        @printf(f, "No branch data available.\n")
        return
    end
    
    # 检查是否有支路功率流数据
    has_flow_data = false
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 18
        has_flow_data = true
    end
    
    total_p_loss = 0.0
    total_q_loss = 0.0
    
    if has_flow_data
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
            
            # 获取功率流数据
            pf = mpc.branchAC[i, 15] 
            qf = mpc.branchAC[i, 16] 
            pt = mpc.branchAC[i, 17] 
            qt = mpc.branchAC[i, 18] 
            
            # 计算损耗
            p_loss = pf + pt
            
            # 获取线路参数
            r = mpc.branchAC[i, 3] 
            x = mpc.branchAC[i, 4]
            
            # 获取实际的母线电压值和相角
            v_from = 1.0  # 默认值，如果找不到实际电压
            v_to = 1.0    # 默认值，如果找不到实际电压
            ang_from = 0.0
            ang_to = 0.0
            
            # 从母线数据中查找实际电压值和相角
            if hasproperty(mpc, :busAC)
                for j in 1:size(mpc.busAC, 1)
                    if Int(mpc.busAC[j, 1]) == from_bus
                        v_from = mpc.busAC[j, 8]  # 使用实际电压幅值
                        ang_from = mpc.busAC[j, 9] * pi/180  # 转换为弧度
                    elseif Int(mpc.busAC[j, 1]) == to_bus
                        v_to = mpc.busAC[j, 8]    # 使用实际电压幅值
                        ang_to = mpc.busAC[j, 9] * pi/180  # 转换为弧度
                    end
                end
            end
            
            # 相角差
            angle_diff = ang_from - ang_to
            
            # 计算线路导纳
            y = 1 / complex(r, x)
            y_abs = abs(y)
            
            # 直接计算线路电流幅值平方
            i_mag_squared = (v_from^2 + v_to^2 - 2*v_from*v_to*cos(angle_diff)) * y_abs^2
            
            # 计算无功损耗 - 使用电抗和电流平方
            q_loss = x * i_mag_squared * baseMVA
            
            total_p_loss += p_loss
            total_q_loss += q_loss
            
            # 打印支路数据
            @printf(f, "%5d  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n", 
                    branch_id, from_bus, to_bus, pf, qf, pt, qt, p_loss, q_loss)
        end
    else
        # 如果没有功率流数据，只打印支路拓扑信息
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
            
            # 打印支路数据（无功率流信息）
                        # 打印支路数据（无功率流信息）
                        @printf(f, "%5d  %5d  %5d  %8s  %8s  %8s  %8s  %8s  %8s\n", 
                        branch_id, from_bus, to_bus, "-", "-", "-", "-", "-", "-")
            end
        end
        
        # 打印总计
        @printf(f, "                                                             --------  --------\n")
        @printf(f, "                                                    Total:   %8.2f  %8.2f\n", 
                total_p_loss, total_q_loss)
        
        write(f, "\n")
    end
    
    """
    从您的数据结构中提取母线数据
    """
    function extract_bus_data(mpc::JPC)
        # 如果已经有bus数据，直接返回
        if hasproperty(mpc, :busAC)
            return mpc.busAC
        end
        
        # 否则，尝试从gen、load和branch数据构建基本的bus数据
        bus_ids = Set{Int}()
        
        # 从gen、load和branch数据中收集所有母线ID
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                push!(bus_ids, Int(mpc.genAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                push!(bus_ids, Int(mpc.loadAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :branchAC)
            for i in 1:size(mpc.branchAC, 1)
                push!(bus_ids, Int(mpc.branchAC[i, 1]))
                push!(bus_ids, Int(mpc.branchAC[i, 2]))
            end
        end
        
        # 创建基本的bus数据矩阵
        # 列：[bus_id, Vm, Va]
        bus_data = zeros(length(bus_ids), 3)
        
        for (i, bus_id) in enumerate(sort(collect(bus_ids)))
            bus_data[i, 1] = bus_id
            bus_data[i, 2] = 1.0  # 默认电压幅值
            bus_data[i, 3] = 0.0  # 默认相角
        end
        
        return bus_data
    end
    
    """
    将潮流计算结果格式化为MATPOWER风格的报告并保存为文本文件
    """
    function generate_matpower_report(mpc::JPC, area, execution_time, isolated, output_file::String="PowerFlow_report.txt")
        # 打开文件用于写入
        open(output_file, "w") do f
            # 写入报告头部
            write(f, "JUPOWER Version 0.01, $(Dates.format(now(), "dd-u-yyyy"))\n")
            write(f, "Power Flow -- AC-polar-power formulation\n\n")
            
            # 写入收敛信息
            write(f, "Newton's method converged in $(mpc.iterationsAC) iterations.\n")
            if mpc.success
                write(f, "PF successful\n\n")
            else
                write(f, "PF NOT successful\n\n")
            end
            
            # 假设计算时间
            write(f, "Converged in $(execution_time) seconds\n")
            
            # 系统摘要
            write_system_summary(f, mpc, area, isolated)
            
            # 母线数据
            write_bus_data(f, mpc, isolated)
            
            # 支路数据
            write_branch_data(f, mpc)
        end
        
        println("报告已保存至 $output_file")
    end
    
