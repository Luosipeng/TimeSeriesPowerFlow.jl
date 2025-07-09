using Test
using PowerFlow

"""
测试不平衡电力潮流分析在非对称负载节点上的功能。
测试包括以下步骤：
1. 加载测试数据
2. 设置电力潮流分析的选项
3. 运行不平衡电力潮流分析
4. 验证结果的正确性
"""
function test_unbalanced_power_flow()
    @testset verbose=true "不平衡电力潮流分析测试" begin
        
        @testset verbose=true "Dyn变压器测试" begin
            # 加载测试数据
            xlsx_file = "data/test_trafo_Dyn.xlsx"
            net = PowerFlow.import_distribution_system_data(xlsx_file)
            
            # 设置电力潮流分析的选项
            opt = PowerFlow.options()
            pf_config = PowerFlow.Powerflow_config()
            opt["PF"]["NR_ALG"] = "bicgstab"
            opt["PF"]["ENFORCE_Q_LIMS"] = 0
            opt["PF"]["baseMVA"] = 100.0
            opt["PF"]["net_hz"] = 50.0
            
            net["baseMVA"] = pf_config["baseMVA"]
            net["mode"] = pf_config["Mode"]
            
            # 运行不平衡电力潮流分析
            result_net = PowerFlow.runupf(net, opt)
            
            # 验证结果
            @test haskey(result_net, "success")
            @test result_net["success"] == true
            
            # 检查是否有电压结果
            @test haskey(result_net, "res_bus_3ph")
            @test "vm_a_pu" in names(result_net["res_bus_3ph"])
            @test "vm_b_pu" in names(result_net["res_bus_3ph"])
            @test "vm_c_pu" in names(result_net["res_bus_3ph"])
            @test "va_a_degree" in names(result_net["res_bus_3ph"])
            @test "va_b_degree" in names(result_net["res_bus_3ph"])
            @test "va_c_degree" in names(result_net["res_bus_3ph"])
            
            # 检查电压幅值是否在合理范围内 (通常在0.9-1.1标幺值范围)
            @test all(0.8 .<= result_net["res_bus_3ph"].vm_a_pu .<= 1.2)
            @test all(0.8 .<= result_net["res_bus_3ph"].vm_b_pu .<= 1.2)
            @test all(0.8 .<= result_net["res_bus_3ph"].vm_c_pu .<= 1.2)
            
            # 检查收敛性
            @test haskey(result_net, "iterations")
            @test result_net["iterations"] < 50  # 确保在合理迭代次数内收敛
        end
        
        @testset verbose=true "YNyn变压器测试" begin
            # 加载测试数据
            xlsx_file = "data/test_trafo_YNyn.xlsx"
            net = PowerFlow.import_distribution_system_data(xlsx_file)
            
            # 设置电力潮流分析的选项
            opt = PowerFlow.options()
            pf_config = PowerFlow.Powerflow_config()
            opt["PF"]["NR_ALG"] = "bicgstab"
            opt["PF"]["ENFORCE_Q_LIMS"] = 0
            opt["PF"]["baseMVA"] = 100.0
            opt["PF"]["net_hz"] = 50.0
            
            net["baseMVA"] = pf_config["baseMVA"]
            net["mode"] = pf_config["Mode"]
            
            # 运行不平衡电力潮流分析
            result_net = PowerFlow.runupf(net, opt)
            
             # 验证结果
             @test haskey(result_net, "success")
             @test result_net["success"] == true
             
             # 检查是否有电压结果
             @test haskey(result_net, "res_bus_3ph")
             @test "vm_a_pu" in names(result_net["res_bus_3ph"])
             @test "vm_b_pu" in names(result_net["res_bus_3ph"])
             @test "vm_c_pu" in names(result_net["res_bus_3ph"])
             @test "va_a_degree" in names(result_net["res_bus_3ph"])
             @test "va_b_degree" in names(result_net["res_bus_3ph"])
             @test "va_c_degree" in names(result_net["res_bus_3ph"])
             
             # 检查电压幅值是否在合理范围内 (通常在0.9-1.1标幺值范围)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_a_pu .<= 1.2)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_b_pu .<= 1.2)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_c_pu .<= 1.2)
             
             # 检查收敛性
             @test haskey(result_net, "iterations")
             @test result_net["iterations"] < 50  # 确保在合理迭代次数内收敛
        end
        
        @testset verbose=true "Yzn变压器测试" begin
            # 加载测试数据
            xlsx_file = "data/test_trafo_Yzn.xlsx"
            net = PowerFlow.import_distribution_system_data(xlsx_file)
            
            # 设置电力潮流分析的选项
            opt = PowerFlow.options()
            pf_config = PowerFlow.Powerflow_config()
            opt["PF"]["NR_ALG"] = "bicgstab"
            opt["PF"]["ENFORCE_Q_LIMS"] = 0
            opt["PF"]["baseMVA"] = 100.0
            opt["PF"]["net_hz"] = 50.0
            
            net["baseMVA"] = pf_config["baseMVA"]
            net["mode"] = pf_config["Mode"]
            
            # 运行不平衡电力潮流分析
            result_net = PowerFlow.runupf(net, opt)
            
             # 验证结果
             @test haskey(result_net, "success")
             @test result_net["success"] == true
             
             # 检查是否有电压结果
             @test haskey(result_net, "res_bus_3ph")
             @test "vm_a_pu" in names(result_net["res_bus_3ph"])
             @test "vm_b_pu" in names(result_net["res_bus_3ph"])
             @test "vm_c_pu" in names(result_net["res_bus_3ph"])
             @test "va_a_degree" in names(result_net["res_bus_3ph"])
             @test "va_b_degree" in names(result_net["res_bus_3ph"])
             @test "va_c_degree" in names(result_net["res_bus_3ph"])
             
             # 检查电压幅值是否在合理范围内 (通常在0.9-1.1标幺值范围)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_a_pu .<= 1.2)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_b_pu .<= 1.2)
             @test all(0.8 .<= result_net["res_bus_3ph"].vm_c_pu .<= 1.2)
             
             # 检查收敛性
             @test haskey(result_net, "iterations")
             @test result_net["iterations"] < 50  # 确保在合理迭代次数内收敛
        end
        
        @testset verbose=true "两母线测试系统" begin
            # 加载测试数据
            xlsx_file = "data/bus_2_test.xlsx"
            net = PowerFlow.import_distribution_system_data(xlsx_file)
            
            # 设置电力潮流分析的选项
            opt = PowerFlow.options()
            pf_config = PowerFlow.Powerflow_config()
            opt["PF"]["NR_ALG"] = "\\"  
            opt["PF"]["ENFORCE_Q_LIMS"] = 0
            opt["PF"]["baseMVA"] = 100.0
            opt["PF"]["net_hz"] = 50.0
            
            net["baseMVA"] = pf_config["baseMVA"]
            net["mode"] = pf_config["Mode"]
            
            # 运行不平衡电力潮流分析
            result_net = PowerFlow.runupf(net, opt)
            
            # 验证结果
            @test haskey(result_net, "success")
            @test result_net["success"] == true
            
            # 检查是否有电压结果
            @test haskey(result_net, "res_bus_3ph")
            @test "vm_a_pu" in names(result_net["res_bus_3ph"])
            @test "vm_b_pu" in names(result_net["res_bus_3ph"])
            @test "vm_c_pu" in names(result_net["res_bus_3ph"])
            
            # 检查电压误差是否小于1e-5
            mis = PowerFlow.check_it(result_net) 
            @test mis < 1e-5

            # 检查收敛性
            @test haskey(result_net, "iterations")
            @test result_net["iterations"] < 50  # 确保在合理迭代次数内收敛
        end

        @testset verbose=true "三相单相平衡潮流测试" begin
            opt = PowerFlow.options()
            pf_config = PowerFlow.Powerflow_config()

            ##============ one phase power flow==================
            xlsx_file = "data/test_trafo_YNyn_1ph_pf.xlsx"
            net_1ph = PowerFlow.import_distribution_system_data(xlsx_file)

            opt["PF"]["NR_ALG"] = "bicgstab";
            opt["PF"]["ENFORCE_Q_LIMS"] = 0;
            opt["PF"]["baseMVA"] = 100.0;
            opt["PF"]["net_hz"] = 50.0;
            pf_config["Mode"] = "1_ph_pf"

            net_1ph["baseMVA"] = pf_config["baseMVA"]
            net_1ph["mode"] = pf_config["Mode"]
            # Run a single phase power flow analysis
            net_1ph = PowerFlow.runnetpf(net_1ph, opt)

            @test haskey(net_1ph, "success")
            @test net_1ph["success"] == true

            @test haskey(net_1ph, "iterations")
            @test net_1ph["iterations"] < 50  # 确保在合理迭代次数内收敛
            ##============ Three phase power flow================
            xlsx_file = "data/test_trafo_YNyn_3ph_pf.xlsx"
            net_3ph = PowerFlow.import_distribution_system_data(xlsx_file)

            pf_config["Mode"] = "3_ph_pf"
            net_3ph["baseMVA"] = pf_config["baseMVA"]
            net_3ph["mode"] = pf_config["Mode"]
            # # # Run the unbalanced power flow analysis
            net_3ph = PowerFlow.runupf(net_3ph, opt)

            @test haskey(net_3ph, "success")
            @test net_3ph["success"] == true

            @test haskey(net_3ph, "iterations")
            @test net_3ph["iterations"] < 50  # 确保在合理迭代次数内收敛

            # ============ Check the results ===================    
            Vm_1ph = net_1ph["res_bus"].vm_pu
            Vm_3ph = net_3ph["res_bus_3ph"].vm_a_pu

            Va_1ph = net_1ph["res_bus"].va_degree
            Va_3ph = net_3ph["res_bus_3ph"].va_a_degree

            mis_Vm = abs.((Vm_1ph - Vm_3ph))
            mis_Va = abs.((Va_1ph - Va_3ph))

            mis_Vm_max = maximum(mis_Vm)
            mis_Va_max = maximum(mis_Va)

            @test mis_Vm_max < 1e-3
            @test mis_Va_max < 1e-3
        end
    end
end

# 运行测试
test_unbalanced_power_flow()
