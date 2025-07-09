"""
    Sd

A mutable structure to store ZIP load model components.

# Fields
- `z`: Complex vector representing constant impedance component of load
- `i`: Complex vector representing constant current component of load
- `p`: Complex vector representing constant power component of load

This structure is used to store the components of the ZIP load model
for efficient power flow calculations.
"""
mutable struct Sd
    z::Vector{ComplexF64}
    i::Vector{ComplexF64}
    p::Vector{ComplexF64}
end

"""
    makeSdzip(baseMVA, bus, pw_1, pw_2, pw_3)

Create a ZIP load model structure from bus data and specified ZIP percentages.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `pw_1`: Vector of constant power percentages for active power
- `pw_2`: Vector of constant current percentages for active power
- `pw_3`: Vector of constant impedance percentages for active power

# Returns
- `sd`: An Sd structure containing the ZIP load model components

# Description
This function creates a ZIP (constant impedance, constant current, constant power) load model
for power flow calculations. It converts the load data from the bus matrix to per-unit
values based on the system MVA base and the specified percentages for each component.

The function assumes that the same percentage distribution applies to both active and reactive power,
so it uses the same percentages (pw_1, pw_2, pw_3) for both P and Q components.

# Notes
- All power values are converted to per-unit on system MVA base
- PD and QD columns in bus represent the total active and reactive power demand
- pw_1, pw_2, pw_3 represent the percentages of constant power, constant current, and constant impedance components
"""
function makeSdzip(baseMVA, bus, pw_1, pw_2, pw_3)

    qw_1 = pw_1
    qw_2 = pw_2
    qw_3 = pw_3
    z = (bus[:, PD] .* pw_3  + 1im * bus[:, QD] .* qw_3) / baseMVA
    i = (bus[:, PD] .* pw_2  + 1im * bus[:, QD] .* qw_2) / baseMVA
    p = (bus[:, PD] .* pw_1  + 1im * bus[:, QD] .* qw_1) / baseMVA
    sd = Sd(z[:,1], i[:,1], p[:,1])
    return sd
end

"""
    makeSdzip(baseMVA, bus)

Create a ZIP load model structure from bus data with default constant power model.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters

# Returns
- `sd`: An Sd structure containing the ZIP load model components

# Description
This is a simplified version of the makeSdzip function that assumes a constant power model
(100% constant power, 0% constant current, 0% constant impedance).

# Notes
- All power values are converted to per-unit on system MVA base
- PD and QD columns in bus represent the total active and reactive power demand
- By default, all load is modeled as constant power (P-Q) load
"""
function makeSdzip(baseMVA, bus)

        pw = [1 0 0]
        qw = pw
    z = (bus[:, PD] * pw[3]  + 1im * bus[:, QD] * qw[3]) / baseMVA
    i = (bus[:, PD] * pw[2]  + 1im * bus[:, QD] * qw[2]) / baseMVA
    p = (bus[:, PD] * pw[1]  + 1im * bus[:, QD] * qw[1]) / baseMVA
    sd = Sd(z, i, p)
    return sd
end
