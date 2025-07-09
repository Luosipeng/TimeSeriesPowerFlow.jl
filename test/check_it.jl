function check_it(net)
    """检查结果是否符合预期"""
    # 获取母线电压数据并过滤掉NaN值
    mask_bus = .!isnan.(net["res_bus_3ph"].vm_a_pu)
    bus_pp = hcat(abs.(net["res_bus_3ph"].vm_a_pu)[mask_bus],
                   abs.(net["res_bus_3ph"].vm_b_pu)[mask_bus],
                   abs.(net["res_bus_3ph"].vm_c_pu)[mask_bus])
    #, "vm_b_pu", "vm_c_pu"
    
    bus_pf = abs.([0.96742893 1.01302766 1.019784;
                   0.74957533 1.09137945 1.05124282])
    
    diff_bus = maximum(abs.(bus_pp .- bus_pf))
    # println("母线电压最大误差: $diff_bus")

    # 获取线路数据并过滤掉NaN值
    # mask_line = .!isnan.(net["res_line_3ph"]["i_a_from_ka"])
    # line_cols = [:i_a_from_ka, :i_a_to_ka, :i_b_from_ka, :i_b_to_ka,
    #              :i_c_from_ka, :i_c_to_ka,
    #              :p_a_from_mw, :p_a_to_mw, :q_a_from_mvar, :q_a_to_mvar,
    #              :p_b_from_mw, :p_b_to_mw, :q_b_from_mvar, :q_b_to_mvar,
    #              :p_c_from_mw, :p_c_to_mw, :q_c_from_mvar, :q_c_to_mvar,
    #              :loading_a_percent, :loading_b_percent, :loading_c_percent,
    #              :loading_percent]
    
    # line_pp = abs.(Matrix(net.res_line_3ph[mask_line, line_cols]))
    
    # line_pf = abs.([1.34212045 1.48537916 0.13715552 0.26009611 
    #                 0.22838401 0.1674634 
    #                 55.70772301 (-49.999992954) 60.797262682 (-49.999959283) 
    #                 8.7799379802 (-9.9999996625) (-0.88093549983) (-15.000000238) 
    #                 9.3739293122 (-10.000000161) (-11.441663679) (-4.9999997418)
    #                 154.2452 27.00894 23.71589 
    #                 154.2452])
    
    # diff_line = maximum(abs.(line_pp .- line_pf))
    # println("线路参数最大误差: $diff_line")

    # if diff_bus < 1.1e-6
    #     println("✓ 检查通过！结果符合预期。")
    # else
    #     println("✗ 检查失败！结果与预期不符。")
    # end

    # # 打印一些关键结果
    # println("\n母线电压结果:")
    # display(hcat(net["res_bus_3ph"].vm_a_pu,
    # net["res_bus_3ph"].vm_b_pu,
    # net["res_bus_3ph"].vm_c_pu))

    # println("\n线路结果:")
    # display(net.res_line_3ph[:, [:i_a_from_ka, :i_b_from_ka, :i_c_from_ka,
    #                             :p_a_from_mw, :q_a_from_mvar,
    #                             :loading_a_percent, :loading_percent]])
    
    return diff_bus
end
