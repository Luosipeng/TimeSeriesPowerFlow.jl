"""
    ext2int(bus::Matrix{Float64}, gen::Matrix{Float64}, branch::Matrix{Float64}, load::Matrix{Float64}, pvarray::Matrix{Float64})

Convert power system data from external to internal format.

This function performs several important preprocessing steps:
1. Filters out out-of-service components (buses, generators, branches, loads, PV arrays)
2. Removes zero-value loads and generation
3. Creates a mapping between external bus numbers and consecutive internal indices
4. Renumbers all buses consecutively and updates all component references accordingly

# Arguments
- `bus::Matrix{Float64}`: Bus data matrix
- `gen::Matrix{Float64}`: Generator data matrix
- `branch::Matrix{Float64}`: Branch data matrix
- `load::Matrix{Float64}`: Load data matrix
- `pvarray::Matrix{Float64}`: PV array data matrix (optional)

# Returns
- Tuple containing the processed matrices (bus, gen, branch, load, pvarray) and the mapping from internal to external bus numbers (i2e)
"""
function ext2int(bus::Matrix{Float64}, gen::Matrix{Float64}, branch::Matrix{Float64}, load::Matrix{Float64}, pvarray::Matrix{Float64})
    
    # Find the in_service buses, generators, branches and loads  
    gen  = gen[gen[:, GEN_STATUS] .!= 0, :]
    bus  = bus[bus[:, BUS_TYPE] .!= 0, :]
    load = load[load[:, LOAD_STATUS] .!= 0, :]
    branch = branch[branch[:, BR_STATUS] .!= 0, :]
    if size(pvarray, 1) > 0
        pvarray = pvarray[pvarray[:, PV_IN_SERVICE] .!= 0, :]
    end

    # remove zero load and generation
    gen = gen[gen[:, 8] .!= 0, :]
    branch = branch[branch[:, 11] .!= 0, :]
    # create map of external bus numbers to bus indices
    i2e = Int.(bus[:, BUS_I])  # ensure i2e is integer type
    e2i = sparsevec(zeros(Int, Int(maximum(i2e))))
    e2i[Int.(i2e)] = 1:size(bus, 1)
    # renumber buses consecutively
    bus[:, BUS_I] = e2i[bus[:, BUS_I]]
    gen[:, GEN_BUS] = e2i[gen[:, GEN_BUS]]
    branch[:, F_BUS] = e2i[branch[:, F_BUS]]
    branch[:, T_BUS] = e2i[branch[:, T_BUS]]
    load[:, LOAD_CND] = e2i[load[:, LOAD_CND]]
    if size(pvarray, 1) > 0
        pvarray[:, PV_BUS] = e2i[pvarray[:, PV_BUS]]
    end
    return bus, gen, branch, load, pvarray, i2e
end

"""
    ext2int(jpc::JPC)

Convert power system data from external to internal format using a JPC structure.

This function creates a deep copy of the input JPC structure and processes its components:
1. Filters out out-of-service components (buses, generators, branches, loads)
2. Removes zero-value loads and generation
3. Creates a mapping between external bus numbers and consecutive internal indices
4. Renumbers all buses consecutively and updates all component references accordingly
5. Updates the JPC structure with the processed data

# Arguments
- `jpc::JPC`: A JPC structure containing power system data

# Returns
- A tuple containing the processed JPC structure and the mapping from internal to external bus numbers (i2e)
"""
function ext2int(jpc::JPC)

    # Create a copy of JPC to avoid modifying the original data
    new_jpc = deepcopy(jpc)
    
    # Get data from JPC
    bus = new_jpc.busAC
    gen = new_jpc.genAC
    branch = new_jpc.branchAC
    load = new_jpc.loadAC
    
    # Find the in_service buses, generators, branches and loads  
    gen = gen[gen[:, GEN_STATUS] .!= 0, :]
    bus = bus[bus[:, BUS_TYPE] .!= 0, :]
    load = load[load[:, LOAD_STATUS] .!= 0, :]
    branch = branch[branch[:, BR_STATUS] .!= 0, :]

    # remove zero load and generation
    gen = gen[gen[:, GEN_STATUS] .!= 0, :]
    branch = branch[branch[:, BR_STATUS] .!= 0, :]
    
    # create map of external bus numbers to bus indices
    i2e = Int.(bus[:, BUS_I])  # ensure i2e is integer type
    max_bus_num = Int(maximum(i2e))
    e2i = sparsevec(zeros(Int, max_bus_num))
    e2i[Int.(i2e)] = 1:size(bus, 1)
    
    # renumber buses consecutively
    bus[:, BUS_I] = e2i[Int.(bus[:, BUS_I])]
    gen[:, GEN_BUS] = e2i[Int.(gen[:, GEN_BUS])]
    branch[:, F_BUS] = e2i[Int.(branch[:, F_BUS])]
    branch[:, T_BUS] = e2i[Int.(branch[:, T_BUS])]
    load[:, LOAD_CND] = e2i[Int.(load[:, LOAD_CND])]
    
    # Update JPC structure
    new_jpc.busAC = bus
    new_jpc.genAC = gen
    new_jpc.branchAC = branch
    new_jpc.loadAC = load
    
    return new_jpc, i2e
end
