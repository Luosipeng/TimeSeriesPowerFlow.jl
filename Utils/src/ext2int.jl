"""
    Convert the data from external format to internal format
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
    i2e = Int.(bus[:, BUS_I])  # 确保i2e是整数类型
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
    Convert the data from external format to internal format
"""
function ext2int(jpc::JPC)

    # 创建JPC的副本，以免修改原始数据
    new_jpc = deepcopy(jpc)
    
    # 获取JPC中的数据
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
    i2e = Int.(bus[:, BUS_I])  # 确保i2e是整数类型
    max_bus_num = Int(maximum(i2e))
    e2i = sparsevec(zeros(Int, max_bus_num))
    e2i[Int.(i2e)] = 1:size(bus, 1)
    
    # renumber buses consecutively
    bus[:, BUS_I] = e2i[Int.(bus[:, BUS_I])]
    gen[:, GEN_BUS] = e2i[Int.(gen[:, GEN_BUS])]
    branch[:, F_BUS] = e2i[Int.(branch[:, F_BUS])]
    branch[:, T_BUS] = e2i[Int.(branch[:, T_BUS])]
    load[:, LOAD_CND] = e2i[Int.(load[:, LOAD_CND])]
    
    # 更新JPC结构
    new_jpc.busAC = bus
    new_jpc.genAC = gen
    new_jpc.branchAC = branch
    new_jpc.loadAC = load
    
    return new_jpc, i2e
end
