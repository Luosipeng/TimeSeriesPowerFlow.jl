using JuMP
using Gurobi

const baseMVA = 10.0

# 创建模型，使用 Gurobi 求解器
model = Model(Gurobi.Optimizer)

# AC 变量
@variable(model, Vm_ac[1:2])
@variable(model, Va_ac[1:2])
@variable(model, P_ac[1:2])
@variable(model, Q_ac[1:2])

# DC 变量
@variable(model, Vm_dc[1:2])
@variable(model, P_dc[1:2])

# 转换器变量
@variable(model, Pac_vsc)
@variable(model, Qac_vsc)
@variable(model, Pdc_vsc)

# 添加参数
G_ac = [1.0 -1.0; -1.0 1.0]  # 使用空格或逗号分隔元素
B_ac = [-1.0 1.0; 1.0 -1.0]   # 同样的修改
G_dc = [0.16975 -0.16975 0.0; -0.16975 2.95245 -2.7827; 0.0 -2.7827 2.7827]
eta = 0.9

# 添加约束
# 添加约束
@constraint(model, Va_ac[1] == 0.0)
@constraint(model, Vm_ac[1] == 1.0)

@constraint(model, P_ac[2] == 0.17/baseMVA)
@constraint(model, Q_ac[2] == 0.1054/baseMVA)

@constraint(model, P_dc[1] == 0.2/baseMVA)
@constraint(model, Vm_dc[2] == 1.0)

# 使用新的非线性表达式接口
@constraint(model, Vm_ac[2]*Vm_ac[1]*(G_ac[2,1]*cos(Va_ac[2]-Va_ac[1]) + B_ac[2,1]*sin(Va_ac[2]-Va_ac[1])) + Vm_ac[2]^2*G_ac[2,2] + P_ac[2] - Pac_vsc == 0)
@constraint(model, Vm_ac[2]*Vm_ac[1]*(G_ac[2,1]*sin(Va_ac[2]-Va_ac[1]) - B_ac[2,1]*cos(Va_ac[2]-Va_ac[1])) + Vm_ac[2]^2*(-B_ac[2,2]) + Q_ac[2] - Qac_vsc == 0)

@constraint(model, Vm_dc[1]*Vm_dc[2]*G_dc[1,2] + Vm_dc[1]^2*G_dc[1,1] + Pdc_vsc + P_dc[1]  == 0)


@constraint(model, Pac_vsc -2*Qac_vsc == 0)
# 添加条件约束：当 Pac_vsc > 0 时，Pdc_vsc = Pac_vsc/0.9；当 Pac_vsc < 0 时，Pdc_vsc = Pac_vsc*0.9
# 引入二进制变量来表示 Pac_vsc 的符号
@variable(model, z, Bin)  # z=1 当 Pac_vsc >= 0, z=0 当 Pac_vsc < 0

# 引入足够大的常数 M (big-M 方法)
M = 100.0  # 选择一个足够大的值，确保能覆盖 Pac_vsc 的可能范围

# 约束 z 的值基于 Pac_vsc 的符号
@constraint(model, Pac_vsc <= M * z)           # 当 Pac_vsc > 0 时，z 必须为 1
@constraint(model, Pac_vsc >= -M * (1 - z))    # 当 Pac_vsc < 0 时，z 必须为 0

# 根据 z 的值确定 Pdc_vsc 的计算方式
@constraint(model, Pdc_vsc <= Pac_vsc/0.9 + M*(1-z))  # 当 z=1 (Pac_vsc >= 0) 时生效
@constraint(model, Pdc_vsc >= Pac_vsc/0.9 - M*(1-z))  # 当 z=1 (Pac_vsc >= 0) 时生效

@constraint(model, Pdc_vsc <= Pac_vsc*0.9 + M*z)      # 当 z=0 (Pac_vsc < 0) 时生效
@constraint(model, Pdc_vsc >= Pac_vsc*0.9 - M*z)      # 当 z=0 (Pac_vsc < 0) 时生效

@constraint(model, Vm_ac .>= 0.9)  # 电压幅值下限
@constraint(model, Vm_ac .<= 1.1)  # 电压幅值上限
@constraint(model, -pi .<= Va_ac .<= pi)  # 相角范围

@constraint(model, Vm_dc .>= 0.9)  # 电压幅值下限
@constraint(model, Vm_dc .<= 1.1)  # 电压幅值上限

@constraint(model, Pac_vsc .<= 0.0135)  # Pac_vsc 上限
@constraint(model, Pac_vsc .>= -0.0135)  # Pac_vsc 下限
@constraint(model, Qac_vsc .<= 0.0135)  # Pac_vsc 上限
@constraint(model, Qac_vsc .>= -0.0135)  # Pac_vsc 下限
# 添加目标函数（这里假设最小化某个目标，您可以根据需要修改）
# 例如，最小化网络损耗或最大化经济效益
# 引入辅助变量
@variable(model, max_value)
# 添加约束来定义 max_value
@constraint(model, max_value >= Vm_ac[2]*Vm_ac[1]*(G_ac[2,1]*cos(Va_ac[2]-Va_ac[1]) + B_ac[2,1]*sin(Va_ac[2]-Va_ac[1])) + Vm_ac[2]^2*G_ac[2,2] + P_ac[2] - Pac_vsc)
@constraint(model, max_value >= Vm_ac[2]*Vm_ac[1]*(G_ac[2,1]*sin(Va_ac[2]-Va_ac[1]) - B_ac[2,1]*cos(Va_ac[2]-Va_ac[1])) + Vm_ac[2]^2*(-B_ac[2,2]) + Q_ac[2] - Qac_vsc)

# 修改目标函数
@objective(model, Min, max_value)
# 求解模型
optimize!(model)

# 输出结果
if termination_status(model) == MOI.OPTIMAL
    println("最优解:")
    println("Vm_ac = ", value.(Vm_ac))
    println("Va_ac = ", value.(Va_ac))
    println("P_ac = ", value.(P_ac))
    println("Q_ac = ", value.(Q_ac))
    println("Vm_dc = ", value.(Vm_dc))
    println("P_dc = ", value.(P_dc))
    println("Pac_vsc = ", value(Pac_vsc))
    println("Qac_vsc = ", value(Qac_vsc))
    println("Pdc_vsc = ", value(Pdc_vsc))
    println("z = ", value(z), " (", value(z) == 1 ? "Pac_vsc >= 0" : "Pac_vsc < 0", ")")
else
    println("求解失败: ", termination_status(model))
end
