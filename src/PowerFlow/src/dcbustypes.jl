"""
    dcbustypes(bus::Matrix{Float64}, gen::Matrix{Float64})

Determine the types of buses in a DC power flow model.

# Arguments
- `bus`: Matrix containing bus data
- `gen`: Matrix containing generator data

# Returns
- `ref`: Index of the reference (slack) bus
- `p`: Indices of P buses (buses with specified power injection)

# Description
This function categorizes buses in a DC power flow model into two types:
1. Reference/slack bus - controls the voltage angle reference
2. P buses - buses with specified active power injection

Unlike AC power flow which has three bus types (PQ, PV, slack), DC power flow
only distinguishes between the reference bus and P buses since voltage magnitudes
are assumed to be 1.0 per unit and reactive power is not modeled.
"""
function dcbustypes(bus::Matrix{Float64}, gen::Matrix{Float64})
    # constants

    # get generator status
    nb = size(bus, 1)
    ng = size(gen, 1)
    Cg = sparse(gen[:, GEN_BUS], 1:ng, gen[:, GEN_STATUS] .> 0, nb, ng)  # gen connection matrix
    bus_gen_status = Cg * ones(ng, 1)  # number of generators at each bus that are ON
    
    # form index lists for slack, PV, and PQ buses
    bus_gen_status = vec(bus_gen_status)
    map!(x -> x != 0 ? true : x, bus_gen_status, bus_gen_status)
    ref = findall(bus[:, BUS_TYPE] .== DC_REF .* bus_gen_status )  # reference bus index
    p  = findall(bus[:, BUS_TYPE] .== P  )  # P bus indices

    return ref, p
end
