"""
    bustypes(bus::Matrix{Float64}, gen::Matrix{Float64})

Determine the types of buses in the power system network.

# Arguments
- `bus`: Matrix containing bus data
- `gen`: Matrix containing generator data

# Returns
- `ref`: Index of the reference (slack) bus
- `pv`: Indices of PV buses (voltage controlled)
- `pq`: Indices of PQ buses (load buses)

# Description
This function categorizes buses into three types:
1. Reference/slack bus - controls both voltage magnitude and angle
2. PV buses - voltage magnitude controlled by generators
3. PQ buses - both active and reactive power are specified
"""
function bustypes(bus::Matrix{Float64}, gen::Matrix{Float64})

    # get generator status
    nb = size(bus, 1)
    ng = size(gen, 1)
    Cg = sparse(gen[:, GEN_BUS], 1:ng, gen[:, GEN_STATUS] .> 0, nb, ng)  # gen connection matrix
    bus_gen_status = Cg * ones(ng, 1)  # number of generators at each bus that are ON
    
    # form index lists for slack, PV, and PQ buses
    bus_gen_status = vec(bus_gen_status)
    map!(x -> x != 0 ? true : x, bus_gen_status, bus_gen_status)
    ref = findall(bus[:, BUS_TYPE] .== REF .* bus_gen_status )  # reference bus index
    pv  = findall(bus[:, BUS_TYPE] .== PV  .* bus_gen_status )  # PV bus indices
    m = bus[:, BUS_TYPE] .== PQ;
    n = bus_gen_status;
    x = ones(size(n));
    y = x-n;
    s = xor(m,y)
    s = Bool.(s)
    s = vec(s)
    pq  = findall(s)  # PQ bus indices
    # pick a new reference bus if for some reason there is none (may have been shut down)
    if isempty(ref)
        ref = pv[1]  # use the first PV bus
        pv = pv[2:end]  # delete it from PV list
    end

    return ref, pv, pq
end

"""
    xor(m, n)

Perform element-wise logical XOR operation on two arrays.

# Arguments
- `m`: First array of boolean or numeric values
- `n`: Second array of boolean or numeric values

# Returns
- Array with 1 where exactly one of the inputs is non-zero, 0 otherwise

# Description
This function implements a custom XOR operation that returns:
- 0 if both inputs at position i are 0
- 1 if exactly one input at position i is non-zero
"""
function xor(m, n)
    l = length(m);
    s = ones(1, l);
    for i = 1:l
        if(m[i] == 0 && n[i] == 0)
            s[i] = 0;
        else
            s[i] = 1;
        end
    end
    return s
end
