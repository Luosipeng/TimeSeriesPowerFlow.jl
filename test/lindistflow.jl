"""
This file is used to test the linear distribution flow model with fixed generation.
Using infinity norm of slack variables as objective function.
Slack variables are unconstrained (can be positive or negative).
"""

using JuMP
using Ipopt
using LinearAlgebra
using Plots

# 定义配电网络参数 - 使用无穷范数目标函数，松弛变量无约束
function create_fixed_gen_distflow_model(
    N::Int,                  # 节点数量
    branches::Vector{Tuple{Int,Int}},  # 支路连接关系 (from_node, to_node)
    r::Vector{Float64},      # 支路电阻
    x::Vector{Float64},      # 支路电抗
    p_load::Vector{Float64}, # 节点有功负荷
    q_load::Vector{Float64}, # 节点无功负荷
    v_0::Float64,            # 根节点电压
    p_gen_fixed::Vector{Float64}, # 节点固定有功发电量
    q_gen_fixed::Vector{Float64}, # 节点固定无功发电量
    v_min::Float64,          # 最小电压限制
    v_max::Float64           # 最大电压限制
)
    # 创建优化模型
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    
    # 构建节点的子节点映射
    children = Dict(i => Int[] for i in 0:N-1)
    for (from, to) in branches
        push!(children[from], to)
    end
    
    # 构建节点的父节点映射
    parent = Dict{Int, Int}()
    for (from, to) in branches
        parent[to] = from
    end
    
    # 支路索引映射
    branch_idx = Dict((from, to) => i for (i, (from, to)) in enumerate(branches))
    
    # 决策变量
    @variable(model, P[1:length(branches)]) # 支路有功功率
    @variable(model, Q[1:length(branches)]) # 支路无功功率
    @variable(model, v_min <= v[1:N-1] <= v_max) # 节点电压幅值的平方
    @variable(model, slack_p[1:N-1])  # 有功功率失配松弛变量，无约束
    @variable(model, slack_q[1:N-1])  # 无功功率失配松弛变量，无约束
    
    # 无穷范数目标函数变量
    @variable(model, max_slack >= 0)
    
    # 目标函数：最小化最大失配量的绝对值（无穷范数）
    @objective(model, Min, max_slack)
    
    # 无穷范数约束（考虑正负值）
    for i in 1:N-1
        @constraint(model, slack_p[i] <= max_slack)
        @constraint(model, -slack_p[i] <= max_slack)
        @constraint(model, slack_q[i] <= max_slack)
        @constraint(model, -slack_q[i] <= max_slack)
    end
    
    # LinDistFlow模型约束
    # 1. 功率平衡约束 (2a) - 使用固定的发电量
    for n in 1:N-1
        # 如果节点有子节点
        if haskey(children, n) && !isempty(children[n])
            child_sum = sum(P[branch_idx[(n, k)]] for k in children[n])
            if haskey(parent, n)
                parent_branch = branch_idx[(parent[n], n)]
                @constraint(model, child_sum == p_gen_fixed[n] - p_load[n] + P[parent_branch] - slack_p[n])
            else
                @constraint(model, child_sum == p_gen_fixed[n] - p_load[n] - slack_p[n])
            end
        # 如果节点是叶节点
        else
            if haskey(parent, n)
                parent_branch = branch_idx[(parent[n], n)]
                @constraint(model, 0 == p_gen_fixed[n] - p_load[n] + P[parent_branch] - slack_p[n])
            else
                @constraint(model, 0 == p_gen_fixed[n] - p_load[n] - slack_p[n])
            end
        end
    end
    
    # 2. 无功功率平衡约束 (2b) - 使用固定的发电量
    for n in 1:N-1
        # 如果节点有子节点
        if haskey(children, n) && !isempty(children[n])
            child_sum = sum(Q[branch_idx[(n, k)]] for k in children[n])
            if haskey(parent, n)
                parent_branch = branch_idx[(parent[n], n)]
                @constraint(model, child_sum == q_gen_fixed[n] - q_load[n] + Q[parent_branch] - slack_q[n])
            else
                @constraint(model, child_sum == q_gen_fixed[n] - q_load[n] - slack_q[n])
            end
        # 如果节点是叶节点
        else
            if haskey(parent, n)
                parent_branch = branch_idx[(parent[n], n)]
                @constraint(model, 0 == q_gen_fixed[n] - q_load[n] + Q[parent_branch] - slack_q[n])
            else
                @constraint(model, 0 == q_gen_fixed[n] - q_load[n] - slack_q[n])
            end
        end
    end
    
    # 3. 电压关系约束 (2c)
    for (i, (from, to)) in enumerate(branches)
        from_voltage = from == 0 ? v_0 : v[from]
        @constraint(model, v[to] == from_voltage - 2*(r[i]*P[i] + x[i]*Q[i]))
    end
    
    return model, P, Q, v, slack_p, slack_q, max_slack
end

# 测试案例：5节点径向配电网，固定发电机出力
function run_fixed_gen_test_case()
    # 网络参数
    N = 5  # 节点数量 (包括根节点0)
    branches = [(0,1), (1,2), (2,3), (3,4)]  # 支路连接
    r = [0.01, 0.01, 0.01, 0.01]  # 支路电阻
    x = [0.02, 0.02, 0.02, 0.02]  # 支路电抗
    p_load = [0.1, 0.2, 0.15, 0.3]  # 节点有功负荷
    q_load = [0.05, 0.1, 0.07, 0.15]  # 节点无功负荷
    v_0 = 1.0  # 根节点电压
    v_min = 0.95  # 最小电压限制
    v_max = 1.05  # 最大电压限制
    
    # 固定的发电机出力 - 这里可以设置您想要测试的值
    p_gen_fixed = [0.1, 0.2, 0.15, 0.1]  # 固定有功发电量
    q_gen_fixed = [0.05, 0.1, 0.07, 0.05]  # 固定无功发电量
    
    # 创建并求解模型
    model, P, Q, v, slack_p, slack_q, max_slack = create_fixed_gen_distflow_model(
        N, branches, r, x, p_load, q_load, v_0, p_gen_fixed, q_gen_fixed, v_min, v_max
    )
    
    optimize!(model)
    
    # 打印结果
    println("优化状态: ", termination_status(model))
    println("\n最优目标值 (最大失配量绝对值): ", objective_value(model))
    
    # 判断潮流是否有解
    if objective_value(model) < 1e-6  # 允许一个小的数值误差
        println("\n结论: 在给定的固定发电机出力下，潮流有解")
    else
        println("\n结论: 在给定的固定发电机出力下，潮流无解或不满足约束")
        println("最大失配量绝对值为 ", objective_value(model), " 单位")
        
        # 找出最大失配量绝对值所在的节点
        abs_p_slack = [abs(value(slack_p[i])) for i in 1:N-1]
        abs_q_slack = [abs(value(slack_q[i])) for i in 1:N-1]
        max_p_slack_idx = argmax(abs_p_slack)
        max_q_slack_idx = argmax(abs_q_slack)
        
        if abs_p_slack[max_p_slack_idx] >= abs_q_slack[max_q_slack_idx]
            println("最大失配量出现在节点 $max_p_slack_idx 的有功功率，值为 $(value(slack_p[max_p_slack_idx]))")
        else
            println("最大失配量出现在节点 $max_q_slack_idx 的无功功率，值为 $(value(slack_q[max_q_slack_idx]))")
        end
    end
    
    println("\n节点固定发电量:")
    for i in 1:N-1
        println("节点 $i: P_gen = $(p_gen_fixed[i]) MW, Q_gen = $(q_gen_fixed[i]) MVar")
    end
    
    println("\n支路功率流:")
    for (i, (from, to)) in enumerate(branches)
        println("支路 $from->$to: P = $(value(P[i])) MW, Q = $(value(Q[i])) MVar")
    end
    
    println("\n节点电压:")
    println("节点 0: V = $v_0 p.u.")
    for i in 1:N-1
        println("节点 $i: V = $(sqrt(value(v[i]))) p.u.")
    end
    
    println("\n功率失配量:")
    for i in 1:N-1
        println("节点 $i: P_slack = $(value(slack_p[i])) MW, Q_slack = $(value(slack_q[i])) MVar")
    end
    
    # 绘制电压分布图
    voltage_values = [v_0; [sqrt(value(v[i])) for i in 1:N-1]]
    nodes = 0:N-1
    
    plt = plot(nodes, voltage_values, 
        marker=:circle, 
        linewidth=2, 
        xlabel="Node Number", 
        ylabel="Voltage (p.u.)",
        title="Distribution Network Voltage Profile with Fixed Generation",
        legend=false,
        ylims=(0.94, 1.06),
        grid=true)
    
    # 添加电压限制线
    hline!([v_min, v_max], linestyle=:dash, color=:red, label="Voltage Limits")
    
    # 绘制失配量分布图（使用绝对值来显示大小）
    p_slack_values = [value(slack_p[i]) for i in 1:N-1]
    q_slack_values = [value(slack_q[i]) for i in 1:N-1]
    
    slack_plt = plot(1:N-1, [p_slack_values q_slack_values], 
        marker=:circle, 
        linewidth=2,
        xlabel="Node Number", 
        ylabel="Slack Value (p.u.)",
        title="Power Mismatch Distribution",
        label=["Active Power Slack" "Reactive Power Slack"],
        grid=true)
    
    # 添加零线
    hline!([0], linestyle=:dash, color=:black, label="Zero Line")
    
    display(plt)
    display(slack_plt)
    
    return model, P, Q, v, slack_p, slack_q, max_slack
end

# 运行测试案例
run_fixed_gen_test_case()
