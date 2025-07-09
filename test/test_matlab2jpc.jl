using PowerFlow
input_path = "C:/Users/13733/Desktop/matpower-8.0/data/case9.m"
output_path ="data/case9.jl"
mpc = PowerFlow.convert_matpower_case(input_path, output_path)