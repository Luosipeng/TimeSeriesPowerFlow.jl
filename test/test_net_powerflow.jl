"""
This file is aimed to design a net-based data structure for power flow analysis.

"""

using  PowerFlow
using  Plots
using  Statistics

# # Set the options for the power flow analysis
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
@time net_1ph = PowerFlow.runnetpf(net_1ph, opt)

##============ Three phase power flow================
xlsx_file = "data/test_trafo_YNyn_3ph_pf.xlsx"
net_3ph = PowerFlow.import_distribution_system_data(xlsx_file)

pf_config["Mode"] = "3_ph_pf"
net_3ph["baseMVA"] = pf_config["baseMVA"]
net_3ph["mode"] = pf_config["Mode"]
# # # Run the unbalanced power flow analysis
@time net_3ph = PowerFlow.runupf(net_3ph, opt)

# ============ Check the results ===================    
Vm_1ph = net_1ph["res_bus"].vm_pu
Vm_3ph = net_3ph["res_bus_3ph"].vm_a_pu

Va_1ph = net_1ph["res_bus"].va_degree
Va_3ph = net_3ph["res_bus_3ph"].va_a_degree

mis_Vm = abs.((Vm_1ph - Vm_3ph))
mis_Va = abs.((Va_1ph - Va_3ph))

@time PowerFlow.error_plot(mis_Vm, mis_Va)