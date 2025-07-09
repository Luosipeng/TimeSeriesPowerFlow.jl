"""
    build_gen(net, jpc)

Build the generator matrix in the MATPOWER-style case structure.

# Arguments
- `net`: Network data structure containing generator information
- `jpc`: MATPOWER-style power flow case structure

# Returns
- `jpc`: Updated jpc structure with generator data
"""
function build_gen(net, jpc)
    # Initialize
    gen_order = Dict()
    f = 1
    # Process generator data
    for element = ["gen", "ext_grid"]
        f = add_gen_order(gen_order, element, net, f)
    end
    # Initialize the generator matrix   
    init_gen(net, jpc, f-1)
    # Process generator data
    for (element, (f, t)) in pairs(gen_order)
        add_element_to_gen(net, jpc, element, f, t)
    end

    return jpc
end

"""
    add_gen_order(gen_order, element, net, f)

Add generator elements to the ordering dictionary and update the counter.

# Arguments
- `gen_order`: Dictionary tracking the position of each element type in the generator matrix
- `element`: Element type (e.g., "gen", "ext_grid")
- `net`: Network data structure
- `f`: Current position counter

# Returns
- Updated position counter
"""
function add_gen_order(gen_order, element, net, f)
    # Process generator data
    if haskey(net, element)
        i = size(net[element], 1)
        gen_order[element] = (f, f + i - 1)
        f += i
    end
    return f
end

"""
    init_gen(net, jpc, f)

Initialize the generator matrix with default values.

# Arguments
- `net`: Network data structure
- `jpc`: MATPOWER-style power flow case structure
- `f`: Total number of generator elements

# Effects
- Sets up the generator matrix in jpc with default limit values
"""
function init_gen(net, jpc, f)
    jpc["gen"] = zeros(f, 21)

    jpc["gen"][:,PMAX] .= 1000000000.0
    jpc["gen"][:,PMIN] .= -1000000000.0
    jpc["gen"][:,QMAX] .= 1000000000.0
    jpc["gen"][:,QMIN] .= -1000000000.0
end

"""
    add_element_to_gen(net, jpc, element, f, t)

Add specific generator element data to the generator matrix.

# Arguments
- `net`: Network data structure
- `jpc`: MATPOWER-style power flow case structure
- `element`: Element type (e.g., "gen", "ext_grid")
- `f`: Start index in the generator matrix
- `t`: End index in the generator matrix

# Throws
- Error if element type is unknown
"""
function add_element_to_gen(net, jpc, element, f, t)
    # Process generator data
    if element == "ext_grid"
        _build_pp_ext_grid(net, jpc, f, t)
    # elseif element == "gen"
    #     _build_pp_gen(net, jpc, f, t)
    # elseif element == "sgen_controllable"
    #     _build_pp_pq_element(net, jpc, "sgen", f, t)
    # elseif element == "load_controllable"
    #     _build_pp_pq_element(net, jpc, "load", f, t, inverted=true)
    # elseif element == "storage_controllable"
    #     _build_pp_pq_element(net, jpc, "storage", f, t, inverted=true)
    # elseif element == "xward"
    #     _build_pp_xward(net, jpc, f, t)
    else
        error("Unknown element $element")
    end
end

"""
    _build_pp_ext_grid(net, jpc, f, t)

Build external grid data in the generator matrix.

# Arguments
- `net`: Network data structure containing external grid information
- `jpc`: MATPOWER-style power flow case structure
- `f`: Start index in the generator matrix
- `t`: End index in the generator matrix

# Effects
- Populates external grid parameters in the generator matrix
"""
function _build_pp_ext_grid(net, jpc, f, t)
    jpc["gen"][f:t, GEN_BUS] = net["ext_grid"].bus
    jpc["gen"][f:t, VG] = net["ext_grid"].vm_pu
    jpc["gen"][f:t, GEN_STATUS] = net["ext_grid"].in_service .== true
    jpc["gen"][f:t, MBASE] .= 100.0

    # jpc["bus"][net["ext_grid"].bus, VM] = net["ext_grid"].vm_pu
    # jpc["bus"][net["ext_grid"].bus, VA] = net["ext_grid"].va_degree

    jpc["gen"][f:t, QMAX] .= 0.0
    jpc["gen"][f:t, QMIN] .= 0.0
end

# """
#     _build_pp_pq_element(net, jpc, element_type, f, t, inverted=false)
#
# Build power injection elements (like static generators) in the generator matrix.
#
# # Arguments
# - `net`: Network data structure
# - `jpc`: MATPOWER-style power flow case structure
# - `element_type`: Type of power injection element (e.g., "sgen")
# - `f`: Start index in the generator matrix
# - `t`: End index in the generator matrix
# - `inverted`: Whether to invert the sign of power values
#
# # Effects
# - Populates power injection parameters in the generator matrix
# """
# function _build_pp_pq_element(net, jpc, "sgen", f, t)
#     (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN,
#      QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, PW_LINEAR, POLYNOMIAL,
#       MODEL, STARTUP, SHUTDOWN, NCOST, COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, GEN_AREA)= PowerFlow.idx_gen()
#
#     jpc["gen"][f:t, GEN_BUS] = net["sgen"].bus
#     jpc["gen"][f:t, PG] = net["sgen"].p_mw
#     jpc["gen"][f:t, QG] = net["sgen"].q_mvar
#     jpc["gen"][f:t, GEN_STATUS] = net["sgen"].in_service .== true
#     jpc["gen"][f:t, MBASE] .= 100.0
# end
