"""
    Sd_gpu

A mutable structure to store ZIP load model components on GPU.

# Fields
- `z`: Complex vector on GPU representing constant impedance component of load
- `i`: Complex vector on GPU representing constant current component of load
- `p`: Complex vector on GPU representing constant power component of load

This structure is used to store the components of the ZIP load model in GPU memory
for efficient power flow calculations.
"""
mutable struct Sd_gpu
    z::CuVector{ComplexF64}
    i::CuVector{ComplexF64}
    p::CuVector{ComplexF64}
end

"""
    makeSdzip_gpu(baseMVA, bus_gpu, pw_1, pw_2, pw_3)

Create a ZIP load model structure on GPU from bus data and ZIP percentages.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus_gpu`: Bus data matrix on GPU with columns representing bus parameters
- `pw_1`: Vector of constant power percentages for active power
- `pw_2`: Vector of constant current percentages for active power
- `pw_3`: Vector of constant impedance percentages for active power

# Returns
- `sd`: An Sd_gpu structure containing the ZIP load model components

# Description
This function creates a ZIP (constant impedance, constant current, constant power) load model
on the GPU for power flow calculations. It converts the load data from the bus matrix to per-unit
values based on the system MVA base and the specified percentages for each component.

The function assumes that the same percentage distribution applies to both active and reactive power,
so it uses the same percentages (pw_1, pw_2, pw_3) for both P and Q components.

# Notes
- All power values are converted to per-unit on system MVA base
- PD and QD columns in bus_gpu represent the total active and reactive power demand
- pw_1, pw_2, pw_3 represent the percentages of constant power, constant current, and constant impedance components
"""
function makeSdzip_gpu(baseMVA, bus_gpu, pw_1, pw_2, pw_3)
    qw_1 = pw_1
    qw_2 = pw_2
    qw_3 = pw_3
    z = (bus_gpu[:, PD] .* pw_3  + 1im * bus_gpu[:, QD] .* qw_3) / baseMVA
    i = (bus_gpu[:, PD] .* pw_2  + 1im * bus_gpu[:, QD] .* qw_2) / baseMVA
    p = (bus_gpu[:, PD] .* pw_1  + 1im * bus_gpu[:, QD] .* qw_1) / baseMVA
    sd = Sd_gpu(z[:,1], i[:,1], p[:,1])
    return sd
end
