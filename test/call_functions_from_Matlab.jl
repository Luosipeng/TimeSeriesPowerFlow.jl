# Test the functions test/call_functions_from_Matlab.jl
using MATLAB
mat"print(pathdef())";

mat"addpath('C:\\Codes\\matpower')"
mat"addpath('C:\\Codes\\matpower\\lib')"
mat"addpath('C:\\Codes\\matpower\\data')"
mat"addpath('C:\\Codes\\matpower\\mips')"
mat"addpath('C:\\Codes\\matpower\\extras')"
mat"addpath('C:\\Codes\\matpower\\mptest')"
mat"addpath('C:\\Codes\\matpower\\mptest\\lib')"
mat"addpath('C:\\Codes\\matpower\\mp-opt-model\\lib\')"
mat"addpath('C:\\Codes\\matpower\\mips\\lib')"
mpc = mat"case14";
mat"result = runpf(case14);";
results = @mget result;

print(results)
