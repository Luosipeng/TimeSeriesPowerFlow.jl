"""
    extract_islands(jpc::JPC)

Extract electrically isolated islands from a power system case.

This function identifies separate electrical islands in the power system and creates
individual JPC objects for each valid island. It also identifies isolated buses that
don't belong to any energized island.

# Arguments
- `jpc::JPC`: Power system case data structure

# Returns
- `Vector{JPC}`: List of JPC objects, each representing an energized island
- `Vector{Int}`: List of isolated buses (buses not connected to any energized island)

# Notes
- An island is considered valid/energized if it contains at least one generator or reference bus
- Islands with only PQ buses and no generation are considered isolated
"""
function extract_islands(jpc::JPC)
    
    # Set up connection matrices
    nb = size(jpc.busAC, 1)     # Number of buses
    nl = size(jpc.branchAC, 1)  # Number of branches
    ng = size(jpc.genAC, 1)     # Number of generators
    nld = size(jpc.loadAC, 1)   # Number of loads

    # Create mapping from external bus numbers to internal indices
    i2e = jpc.busAC[:, BUS_I]
    e2i = Dict{Int, Int}()
    for i in 1:nb
        e2i[Int(i2e[i])] = i
    end
    
    # Find all islands
    groups, isolated = find_islands(jpc)
    
    # Extract each island
    jpc_list = JPC[]
    
    # Create a set to record all bus numbers in valid islands
    all_valid_island_buses = Set{Int}()
    
    # Process each island
    for i in eachindex(groups)
        # Get external bus numbers in island i
        b_external = groups[i]
        
        # Check if the island has generators or reference nodes (i.e., if it is energized)
        has_power_source = false
        for bus_id in b_external
            # Find the corresponding bus row in jpc.busAC
            bus_row = findfirst(x -> Int(x) == bus_id, jpc.busAC[:, BUS_I])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, BUS_TYPE])
                if bus_type == PV || bus_type == REF
                    has_power_source = true
                    break
                end
            end
            
            # Check if there are generators connected to this bus
            for j in 1:ng
                gen_bus = Int(jpc.genAC[j, GEN_BUS])
                if gen_bus == bus_id && jpc.genAC[j, GEN_STATUS] == 1
                    has_power_source = true
                    break
                end
            end
            
            if has_power_source
                break
            end
        end
        
        # If the island has no power source, add its buses to isolated and skip further processing
        if !has_power_source
            append!(isolated, b_external)
            continue
        end
        
        # Record buses in valid islands
        union!(all_valid_island_buses, b_external)
        
        # Convert external bus numbers to internal indices
        b_internal = Int[]
        for bus_id in b_external
            if haskey(e2i, bus_id)
                push!(b_internal, e2i[bus_id])
            else
                @warn "Bus number $bus_id is not in the system"
            end
        end
        
        # Find branches with both ends in island i
        ibr = Int[]
        for j in 1:nl
            f_bus = Int(jpc.branchAC[j, F_BUS])
            t_bus = Int(jpc.branchAC[j, T_BUS])
            if (f_bus in b_external) && (t_bus in b_external)
                push!(ibr, j)
            end
        end
        
        # Find generators connected to buses in island i
        ig = Int[]
        for j in 1:ng
            gen_bus = Int(jpc.genAC[j, GEN_BUS])
            if gen_bus in b_external
                push!(ig, j)
            end
        end
        
        # Find loads connected to buses in island i
        ild = Int[]
        for j in 1:nld
            load_bus = Int(jpc.loadAC[j, LOAD_CND])  # Use LOAD_CND constant for load connection bus
            if load_bus in b_external
                push!(ild, j)
            end
        end
        
        # Find flexible loads connected to buses in island i
        ild_flex = Int[]
        if size(jpc.loadAC_flex, 1) > 0
            for j in 1:size(jpc.loadAC_flex, 1)
                load_bus = Int(jpc.loadAC_flex[j, 1])  # Assume first column is bus number
                if load_bus in b_external
                    push!(ild_flex, j)
                end
            end
        end
        
        # Find asymmetric loads connected to buses in island i
        ild_asymm = Int[]
        if size(jpc.loadAC_asymm, 1) > 0
            for j in 1:size(jpc.loadAC_asymm, 1)
                load_bus = Int(jpc.loadAC_asymm[j, 1])  # Assume first column is bus number
                if load_bus in b_external
                    push!(ild_asymm, j)
                end
            end
        end
        
        # Find three-phase branches with both ends in island i
        ibr3ph = Int[]
        if size(jpc.branch3ph, 1) > 0
            for j in 1:size(jpc.branch3ph, 1)
                f_bus = Int(jpc.branch3ph[j, 1])  # Assume first column is from bus
                t_bus = Int(jpc.branch3ph[j, 2])  # Assume second column is to bus
                if (f_bus in b_external) && (t_bus in b_external)
                    push!(ibr3ph, j)
                end
            end
        end
        
        # Find DC buses in island i
        bdc_external = Int[]
        bdc_internal = Int[]
        if size(jpc.busDC, 1) > 0
            for j in 1:size(jpc.busDC, 1)
                bus_id = Int(jpc.busDC[j, 1])  # Assume first column is DC bus number
                # Need to determine relationship between DC buses and AC buses
                if bus_id in b_external
                    push!(bdc_external, bus_id)
                    push!(bdc_internal, j)
                end
            end
        end
        
        # Find DC branches with both ends in island i
        ibrdc = Int[]
        if size(jpc.branchDC, 1) > 0
            for j in 1:size(jpc.branchDC, 1)
                f_bus = Int(jpc.branchDC[j, 1])  # Assume first column is from bus
                t_bus = Int(jpc.branchDC[j, 2])  # Assume second column is to bus
                if (f_bus in bdc_external) && (t_bus in bdc_external)
                    push!(ibrdc, j)
                end
            end
        end
        
        # Find AC distributed generators connected to buses in island i
        isgen = Int[]
        if size(jpc.sgenAC, 1) > 0
            for j in 1:size(jpc.sgenAC, 1)
                sgen_bus = Int(jpc.sgenAC[j, 1])  # Assume first column is bus number
                if sgen_bus in b_external
                    push!(isgen, j)
                end
            end
        end
        
        # Find storage systems connected to buses in island i
        istorage = Int[]
        if size(jpc.storage, 1) > 0
            for j in 1:size(jpc.storage, 1)
                storage_bus = Int(jpc.storage[j, 1])  # Assume first column is bus number
                if storage_bus in b_external
                    push!(istorage, j)
                end
            end
        end
        
        # Find DC distributed generators connected to DC buses in island i
        isgendc = Int[]
        if size(jpc.sgenDC, 1) > 0
            for j in 1:size(jpc.sgenDC, 1)
                sgen_bus = Int(jpc.sgenDC[j, 1])  # Assume first column is DC bus number
                if sgen_bus in bdc_external
                    push!(isgendc, j)
                end
            end
        end
        
        # Find converters connecting buses in island i
        iconv = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])  # Assume first column is AC bus number
                dc_bus = Int(jpc.converter[j, 2])  # Assume second column is DC bus number
                if (ac_bus in b_external) && (dc_bus in bdc_external)
                    push!(iconv, j)
                end
            end
        end
        
        # Find external grids connected to buses in island i
        iext = Int[]
        if size(jpc.ext_grid, 1) > 0
            for j in 1:size(jpc.ext_grid, 1)
                ext_bus = Int(jpc.ext_grid[j, 1])  # Assume first column is bus number
                if ext_bus in b_external
                    push!(iext, j)
                end
            end
        end
        
        # Find high voltage circuit breakers with both ends in island i
        ihvcb = Int[]
        if size(jpc.hvcb, 1) > 0
            for j in 1:size(jpc.hvcb, 1)
                f_bus = Int(jpc.hvcb[j, 1])  # Assume first column is from bus
                t_bus = Int(jpc.hvcb[j, 2])  # Assume second column is to bus
                if (f_bus in b_external) && (t_bus in b_external)
                    push!(ihvcb, j)
                end
            end
        end
        
        # Find microgrids in island i
        img = Int[]
        if size(jpc.microgrid, 1) > 0
            for j in 1:size(jpc.microgrid, 1)
                mg_bus = Int(jpc.microgrid[j, 1])  # Assume first column is bus number
                if mg_bus in b_external
                    push!(img, j)
                end
            end
        end
        
        # Find PV systems connected to buses in island i
        ipv = Int[]
        if isdefined(jpc, :pv) && size(jpc.pv, 1) > 0
            for j in 1:size(jpc.pv, 1)
                pv_bus = Int(jpc.pv[j, 2])  # Assume second column is bus number
                if pv_bus in b_external
                    # If status field exists, check if device is in service
                    is_in_service = size(jpc.pv, 2) >= 9 ? jpc.pv[j, 9] == 1 : true
                    
                    if is_in_service
                        push!(ipv, j)
                    end
                end
            end
        end
        
        # Find AC PV systems connected to buses in island i
        ipv_ac = Int[]
        if isdefined(jpc, :pv_acsystem) && size(jpc.pv_acsystem, 1) > 0
            for j in 1:size(jpc.pv_acsystem, 1)
                pv_bus = Int(jpc.pv_acsystem[j, 2])  # Assume second column is bus number
                if pv_bus in b_external
                    # If status field exists (assume column 9), check if device is in service
                    is_in_service = size(jpc.pv_acsystem, 2) >= 9 ? jpc.pv_acsystem[j, 9] == 1 : true
                    
                    if is_in_service
                        push!(ipv_ac, j)
                    end
                end
            end
        end
        
        # Create a copy of JPC for this island
        jpck = JPC(jpc.version, jpc.baseMVA)
        
        # Copy relevant data
        if !isempty(b_internal)
            jpck.busAC = jpc.busAC[b_internal, :]
        end
        
        if !isempty(ibr)
            jpck.branchAC = jpc.branchAC[ibr, :]
        end
        
        if !isempty(ig)
            jpck.genAC = jpc.genAC[ig, :]
        end
        
        if !isempty(ild)
            jpck.loadAC = jpc.loadAC[ild, :]
        end
        
        if !isempty(ild_flex)
            jpck.loadAC_flex = jpc.loadAC_flex[ild_flex, :]
        end
        
        if !isempty(ild_asymm)
            jpck.loadAC_asymm = jpc.loadAC_asymm[ild_asymm, :]
        end
        
        if !isempty(ibr3ph)
            jpck.branch3ph = jpc.branch3ph[ibr3ph, :]
        end
        
        if !isempty(bdc_internal)
            jpck.busDC = jpc.busDC[bdc_internal, :]
        end
        
        if !isempty(ibrdc)
            jpck.branchDC = jpc.branchDC[ibrdc, :]
        end
        
        if !isempty(isgen)
            jpck.sgenAC = jpc.sgenAC[isgen, :]
        end
        
        if !isempty(istorage)
            jpck.storage = jpc.storage[istorage, :]
        end
        
        if !isempty(isgendc)
            jpck.sgenDC = jpc.sgenDC[isgendc, :]
        end
        
        if !isempty(iconv)
            jpck.converter = jpc.converter[iconv, :]
        end
        
        if !isempty(iext)
            jpck.ext_grid = jpc.ext_grid[iext, :]
        end
        
        if !isempty(ihvcb)
            jpck.hvcb = jpc.hvcb[ihvcb, :]
        end
        
        if !isempty(img)
            jpck.microgrid = jpc.microgrid[img, :]
        end
        
        # Copy PV data
        if !isempty(ipv) && isa(jpc.pv, Array)
            jpck.pv = copy(jpc.pv[ipv, :])
        else
            # If no PV devices are connected to this island, create an empty array
            if isa(jpc.pv, Array)
                jpck.pv = similar(jpc.pv, 0, size(jpc.pv, 2))
            else
                jpck.pv = deepcopy(jpc.pv)  # If pv is not an array, copy directly
            end
        end
        
        # Copy AC PV system data
        if !isempty(ipv_ac) && isdefined(jpc, :pv_acsystem)
            jpck.pv_acsystem = jpc.pv_acsystem[ipv_ac, :]
        else
            # If no AC PV systems are connected to this island, create an empty array
            if isdefined(jpc, :pv_acsystem) && isa(jpc.pv_acsystem, Array)
                jpck.pv_acsystem = similar(jpc.pv_acsystem, 0, size(jpc.pv_acsystem, 2))
            elseif isdefined(jpc, :pv_acsystem)
                jpck.pv_acsystem = deepcopy(jpc.pv_acsystem)  # If pv_acsystem is not an array, copy directly
            end
        end
        
        push!(jpc_list, jpck)
    end
    
    # Check if there are buses not in any valid island (i.e., not in all_valid_island_buses and not in isolated)
    for i in 1:nb
        bus_id = Int(jpc.busAC[i, BUS_I])
        if !(bus_id in all_valid_island_buses) && !(bus_id in isolated)
            push!(isolated, bus_id)
        end
    end
    
    return jpc_list, isolated
end


"""
    extract_islands_acdc(jpc::JPC)

Extract electrically isolated islands from a hybrid AC-DC power system case.

This function identifies separate electrical islands in a hybrid AC-DC power system
and creates individual JPC objects for each valid island. It handles the complexity
of interconnected AC and DC subsystems through converters.

# Arguments
- `jpc::JPC`: Power system case data structure containing both AC and DC components

# Returns
- `Vector{JPC}`: List of JPC objects, each representing an energized island
- `Vector{Int}`: List of isolated AC buses (buses not connected to any energized island)

# Notes
- An island is considered valid/energized if it contains at least one generator, reference bus,
  or has a power source in its connected DC subsystem
- Islands with only PQ buses and no generation capability are considered isolated
- The function checks for potential issues like multiple reference nodes or batteries in
  constant Vdc mode that might cause power flow calculation errors
"""
function extract_islands_acdc(jpc::JPC)
    # Set up connection matrices
    nb_ac = size(jpc.busAC, 1)     # Number of AC buses
    nb_dc = size(jpc.busDC, 1)     # Number of DC buses
    nl_ac = size(jpc.branchAC, 1)  # Number of AC branches
    nl_dc = size(jpc.branchDC, 1)  # Number of DC branches
    ng = size(jpc.genAC, 1)        # Number of generators
    nld = size(jpc.loadAC, 1)      # Number of loads
    
    # Create mapping from external bus numbers to internal indices
    i2e_ac = jpc.busAC[:, BUS_I]
    e2i_ac = Dict{Int, Int}()
    for i in 1:nb_ac
        e2i_ac[Int(i2e_ac[i])] = i
    end
    
    i2e_dc = nb_dc > 0 ? jpc.busDC[:, 1] : Int[]  # Assume first column is DC bus number
    e2i_dc = Dict{Int, Int}()
    for i in 1:nb_dc
        e2i_dc[Int(i2e_dc[i])] = i
    end
    
    # Find all islands
    groups, isolated = find_islands_acdc(jpc)
    
    # Extract each island
    jpc_list = JPC[]
    
    # Create a set to record all bus numbers in valid islands
    all_valid_island_buses = Set{Int}()
    
    # Process each island
    for i in eachindex(groups)
        # Get external AC bus numbers in island i
        b_external_ac = groups[i]
        
        # Check if the island has generators or reference nodes (i.e., if it is energized)
        has_power_source = false
        for bus_id in b_external_ac
            # Find the corresponding bus row in jpc.busAC
            bus_row = findfirst(x -> Int(x) == bus_id, jpc.busAC[:, BUS_I])
            if bus_row !== nothing
                bus_type = Int(jpc.busAC[bus_row, BUS_TYPE])
                if bus_type == PV || bus_type == REF
                    has_power_source = true
                    break
                end
            end
            
            # Check if there are generators connected to this bus
            for j in 1:ng
                gen_bus = Int(jpc.genAC[j, GEN_BUS])
                if gen_bus == bus_id && jpc.genAC[j, GEN_STATUS] == 1
                    has_power_source = true
                    break
                end
            end
            
            if has_power_source
                break
            end
        end
        
        # Find DC buses connected to these AC buses
        b_external_dc = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])
                dc_bus = Int(jpc.converter[j, 2])
                # Assume column 3 is status field, if no status field, assume it's in service
                is_active = size(jpc.converter, 2) >= 3 ? jpc.converter[j, 3] == 1 : true
                
                if ac_bus in b_external_ac && is_active
                    push!(b_external_dc, dc_bus)
                end
            end
        end
        
        # Expand DC bus set through DC branches
        if !isempty(b_external_dc) && nl_dc > 0
            # Use breadth-first search to find all connected DC buses
            visited_dc = Set(b_external_dc)
            queue = copy(b_external_dc)
            
            while !isempty(queue)
                current_bus = popfirst!(queue)
                
                for j in 1:nl_dc
                    branch = jpc.branchDC[j, :]
                    # Assume DC branch structure is similar to AC branch
                    if branch[11] != 1  # Check status
                        continue
                    end
                    
                    f_bus = Int(branch[1])
                    t_bus = Int(branch[2])
                    
                    if f_bus == current_bus && !(t_bus in visited_dc)
                        push!(visited_dc, t_bus)
                        push!(queue, t_bus)
                    elseif t_bus == current_bus && !(f_bus in visited_dc)
                        push!(visited_dc, f_bus)
                        push!(queue, f_bus)
                    end
                end
            end
            
            b_external_dc = collect(visited_dc)
        end
        
        # Check if DC system has power sources
        if !has_power_source && !isempty(b_external_dc)
            # Check if there are reference nodes (type 2) in DC buses
            for bus_id in b_external_dc
                # Find the corresponding bus row in jpc.busDC
                bus_row = findfirst(x -> Int(x) == bus_id, jpc.busDC[:, 1])
                if bus_row !== nothing && size(jpc.busDC, 2) >= 2
                    bus_type = Int(jpc.busDC[bus_row, 2])  # Assume column 2 is bus type
                    if bus_type == 2  # DC reference node type is 2
                        has_power_source = true
                        break
                    end
                end
            end
            
            # Check DC distributed generation
            if !has_power_source && size(jpc.sgenDC, 1) > 0
                for j in 1:size(jpc.sgenDC, 1)
                    sgen_bus = Int(jpc.sgenDC[j, 1])
                    # Assume column 3 is status field
                    is_active = size(jpc.sgenDC, 2) >= 3 ? jpc.sgenDC[j, 3] == 1 : true
                    
                    if sgen_bus in b_external_dc && is_active
                        has_power_source = true
                        break
                    end
                end
            end
            
            # Check DC generators
            if !has_power_source && isdefined(jpc, :genDC) && size(jpc.genDC, 1) > 0
                for j in 1:size(jpc.genDC, 1)
                    gen_bus = Int(jpc.genDC[j, 1])
                    is_active = size(jpc.genDC, 2) >= 3 ? jpc.genDC[j, 3] == 1 : true
                    
                    if gen_bus in b_external_dc && is_active
                        has_power_source = true
                        break
                    end
                end
            end
        end
        
        # If the island has no power source, add its buses to isolated and skip further processing
        if !has_power_source
            append!(isolated, b_external_ac)
            continue
        end
        
        # Record buses in valid islands
        union!(all_valid_island_buses, b_external_ac)
        
        # Convert external bus numbers to internal indices
        b_internal_ac = Int[]
        for bus_id in b_external_ac
            if haskey(e2i_ac, bus_id)
                push!(b_internal_ac, e2i_ac[bus_id])
            else
                @warn "AC bus number $bus_id is not in the system"
            end
        end
        
        b_internal_dc = Int[]
        for bus_id in b_external_dc
            if haskey(e2i_dc, bus_id)
                push!(b_internal_dc, e2i_dc[bus_id])
            else
                @warn "DC bus number $bus_id is not in the system"
            end
        end
        
        # Find AC branches with both ends in island i
        ibr_ac = Int[]
        for j in 1:nl_ac
            f_bus = Int(jpc.branchAC[j, F_BUS])
            t_bus = Int(jpc.branchAC[j, T_BUS])
            if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                push!(ibr_ac, j)
            end
        end
        
        # Find DC branches with both ends in island i
        ibr_dc = Int[]
        for j in 1:nl_dc
            f_bus = Int(jpc.branchDC[j, 1])  # Assume first column is from bus
            t_bus = Int(jpc.branchDC[j, 2])  # Assume second column is to bus
            if (f_bus in b_external_dc) && (t_bus in b_external_dc)
                push!(ibr_dc, j)
            end
        end
        
        # Find converters connecting AC and DC buses in island i
        iconv = Int[]
        if size(jpc.converter, 1) > 0
            for j in 1:size(jpc.converter, 1)
                ac_bus = Int(jpc.converter[j, 1])
                dc_bus = Int(jpc.converter[j, 2])
                if (ac_bus in b_external_ac) && (dc_bus in b_external_dc)
                    push!(iconv, j)
                end
            end
        end
        
        # Find generators connected to buses in island i
        ig = Int[]
        for j in 1:ng
            gen_bus = Int(jpc.genAC[j, GEN_BUS])
            if gen_bus in b_external_ac
                push!(ig, j)
            end
        end
        
        # Find DC generators connected to DC buses in island i
        igen_dc = Int[]
        if isdefined(jpc, :genDC) && size(jpc.genDC, 1) > 0
            for j in 1:size(jpc.genDC, 1)
                gen_bus = Int(jpc.genDC[j, 1])  # Assume first column is bus number
                if gen_bus in b_external_dc
                    push!(igen_dc, j)
                end
            end
        end
        
        # Find loads connected to buses in island i
        ild = Int[]
        for j in 1:nld
            load_bus = Int(jpc.loadAC[j, LOAD_CND])  # Use LOAD_CND constant for load connection bus
            if load_bus in b_external_ac
                push!(ild, j)
            end
        end
        
        # Find PV devices connected to buses in island i
        ipv = Int[]
        if isa(jpc.pv, Array) && size(jpc.pv, 1) > 0
            for j in 1:size(jpc.pv, 1)
                pv_bus = Int(jpc.pv[j, 2])  # Use column 2 as bus number
                
                # Check if PV device is connected to DC bus
                if pv_bus in b_external_dc
                    # If status field exists (assume column 8 or other position), check if device is in service
                    is_in_service = size(jpc.pv, 2) >= 8 ? jpc.pv[j, 8] == 1 : true
                    
                    if is_in_service
                        push!(ipv, j)
                    end
                end
            end
        end
        
        # Find AC PV systems connected to buses in island i
        ipv_ac = Int[]
        if isdefined(jpc, :pv_acsystem) && size(jpc.pv_acsystem, 1) > 0
            for j in 1:size(jpc.pv_acsystem, 1)
                pv_bus = Int(jpc.pv_acsystem[j, PV_AC_BUS])  # Use PV_AC_BUS constant to get bus number
                if pv_bus in b_external_ac
                    # If status field exists, check if device is in service
                    is_in_service = jpc.pv_acsystem[j, PV_AC_IN_SERVICE] == 1
                    
                    if is_in_service
                        push!(ipv_ac, j)
                    end
                end
            end
        end
        
                # Find flexible loads connected to buses in island i
        ild_flex = Int[]
        if size(jpc.loadAC_flex, 1) > 0
            for j in 1:size(jpc.loadAC_flex, 1)
                load_bus = Int(jpc.loadAC_flex[j, 1])  # Assume first column is bus number
                if load_bus in b_external_ac
                    push!(ild_flex, j)
                end
            end
        end
        
        # Find asymmetric loads connected to buses in island i
        ild_asymm = Int[]
        if size(jpc.loadAC_asymm, 1) > 0
            for j in 1:size(jpc.loadAC_asymm, 1)
                load_bus = Int(jpc.loadAC_asymm[j, 1])  # Assume first column is bus number
                if load_bus in b_external_ac
                    push!(ild_asymm, j)
                end
            end
        end
        
        # Find three-phase branches with both ends in island i
        ibr3ph = Int[]
        if size(jpc.branch3ph, 1) > 0
            for j in 1:size(jpc.branch3ph, 1)
                f_bus = Int(jpc.branch3ph[j, 1])  # Assume first column is from bus
                t_bus = Int(jpc.branch3ph[j, 2])  # Assume second column is to bus
                if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                    push!(ibr3ph, j)
                end
            end
        end
        
        # Find DC loads connected to DC buses in island i
        ild_dc = Int[]
        if size(jpc.loadDC, 1) > 0
            for j in 1:size(jpc.loadDC, 1)
                load_bus = Int(jpc.loadDC[j, 2])  # Assume second column is bus number
                if load_bus in b_external_dc
                    push!(ild_dc, j)
                end
            end
        end
        
        # Find DC distributed generators connected to DC buses in island i
        isgen_dc = Int[]
        if size(jpc.sgenDC, 1) > 0
            for j in 1:size(jpc.sgenDC, 1)
                sgen_bus = Int(jpc.sgenDC[j, 1])  # Assume first column is bus number
                if sgen_bus in b_external_dc
                    push!(isgen_dc, j)
                end
            end
        end
        
        # Find AC distributed generators connected to buses in island i
        isgen_ac = Int[]
        if size(jpc.sgenAC, 1) > 0
            for j in 1:size(jpc.sgenAC, 1)
                sgen_bus = Int(jpc.sgenAC[j, 1])  # Assume first column is bus number
                if sgen_bus in b_external_ac
                    push!(isgen_ac, j)
                end
            end
        end
        
        # Find storage systems connected to DC buses in island i
        istorage = Int[]
        if size(jpc.storage, 1) > 0
            # Determine the column index for DC bus number
            ess_bus_col = -1
            # Try to find column named ESS_BUS (if column names exist)
            if isdefined(jpc.storage, :colnames) && :ESS_BUS in jpc.storage.colnames
                ess_bus_col = findfirst(x -> x == :ESS_BUS, jpc.storage.colnames)
            else
                # If no column names, assume column 1 is the DC bus number
                ess_bus_col = 1
            end
            
            if ess_bus_col > 0 && ess_bus_col <= size(jpc.storage, 2)
                for j in 1:size(jpc.storage, 1)
                    storage_bus = Int(jpc.storage[j, ess_bus_col])
                    if storage_bus in b_external_dc
                        push!(istorage, j)
                    end
                end
            end
        end
        
        # Find external grids connected to buses in island i
        iext = Int[]
        if size(jpc.ext_grid, 1) > 0
            for j in 1:size(jpc.ext_grid, 1)
                ext_bus = Int(jpc.ext_grid[j, 1])  # Assume first column is bus number
                if ext_bus in b_external_ac
                    push!(iext, j)
                end
            end
        end
        
        # Find high voltage circuit breakers with both ends in island i
        ihvcb = Int[]
        if size(jpc.hvcb, 1) > 0
            for j in 1:size(jpc.hvcb, 1)
                f_bus = Int(jpc.hvcb[j, 1])  # Assume first column is from bus
                t_bus = Int(jpc.hvcb[j, 2])  # Assume second column is to bus
                if (f_bus in b_external_ac) && (t_bus in b_external_ac)
                    push!(ihvcb, j)
                end
            end
        end
        
        # Find microgrids in island i
        img = Int[]
        if size(jpc.microgrid, 1) > 0
            for j in 1:size(jpc.microgrid, 1)
                mg_bus = Int(jpc.microgrid[j, 1])  # Assume first column is bus number
                if mg_bus in b_external_ac
                    push!(img, j)
                end
            end
        end
        
        # Create a copy of JPC for this island
        jpck = JPC(jpc.version, jpc.baseMVA)
        
        # Copy relevant data
        if !isempty(b_internal_ac)
            jpck.busAC = deepcopy(jpc.busAC[b_internal_ac, :])
        end
        
        if !isempty(b_internal_dc)
            jpck.busDC = deepcopy(jpc.busDC[b_internal_dc, :])
            
            # Ensure reference node type in DC system remains 2 instead of 3
            if size(jpck.busDC, 1) > 0 && size(jpck.busDC, 2) >= 2
                for j in 1:size(jpck.busDC, 1)
                    if jpck.busDC[j, 2] == 3  # If incorrectly set to 3
                        jpck.busDC[j, 2] = 2  # Change it back to 2
                    end
                end
            end
        end
        
        if !isempty(ibr_ac)
            jpck.branchAC = deepcopy(jpc.branchAC[ibr_ac, :])
        end
        
        if !isempty(ibr_dc)
            jpck.branchDC = deepcopy(jpc.branchDC[ibr_dc, :])
        end
        
        if !isempty(iconv)
            jpck.converter = deepcopy(jpc.converter[iconv, :])
        end
        
        if !isempty(ig)
            jpck.genAC = deepcopy(jpc.genAC[ig, :])
        end
        
        if !isempty(igen_dc) && isdefined(jpc, :genDC)
            jpck.genDC = deepcopy(jpc.genDC[igen_dc, :])
        end
        
        if !isempty(ild)
            jpck.loadAC = deepcopy(jpc.loadAC[ild, :])
        end
        
        if !isempty(ild_dc)
            jpck.loadDC = deepcopy(jpc.loadDC[ild_dc, :])
        end
        
        # Copy PV data
        if !isempty(ipv) && isa(jpc.pv, Array)
            jpck.pv = deepcopy(jpc.pv[ipv, :])
        else
            # If no PV devices are connected to this island, create an empty array
            if isa(jpc.pv, Array)
                jpck.pv = similar(jpc.pv, 0, size(jpc.pv, 2))
            else
                jpck.pv = deepcopy(jpc.pv)  # If pv is not an array, copy directly
            end
        end
        
        # Copy AC PV system data
        if !isempty(ipv_ac) && isdefined(jpc, :pv_acsystem)
            jpck.pv_acsystem = deepcopy(jpc.pv_acsystem[ipv_ac, :])
        else
            # If no AC PV systems are connected to this island, create an empty array
            if isdefined(jpc, :pv_acsystem) && isa(jpc.pv_acsystem, Array)
                jpck.pv_acsystem = similar(jpc.pv_acsystem, 0, size(jpc.pv_acsystem, 2))
            elseif isdefined(jpc, :pv_acsystem)
                jpck.pv_acsystem = deepcopy(jpc.pv_acsystem)  # If pv_acsystem is not an array, copy directly
            end
        end
        
        if !isempty(ild_flex)
            jpck.loadAC_flex = deepcopy(jpc.loadAC_flex[ild_flex, :])
        end
        
        if !isempty(ild_asymm)
            jpck.loadAC_asymm = deepcopy(jpc.loadAC_asymm[ild_asymm, :])
        end
        
        if !isempty(ibr3ph)
            jpck.branch3ph = deepcopy(jpc.branch3ph[ibr3ph, :])
        end
        
        if !isempty(isgen_ac)
            jpck.sgenAC = deepcopy(jpc.sgenAC[isgen_ac, :])
        end
        
        if !isempty(isgen_dc)
            jpck.sgenDC = deepcopy(jpc.sgenDC[isgen_dc, :])
        end
        
        if !isempty(istorage)
            jpck.storage = deepcopy(jpc.storage[istorage, :])
        end
        
        if !isempty(iext)
            jpck.ext_grid = deepcopy(jpc.ext_grid[iext, :])
        end
        
        if !isempty(ihvcb)
            jpck.hvcb = deepcopy(jpc.hvcb[ihvcb, :])
        end
        
        if !isempty(img)
            jpck.microgrid = deepcopy(jpc.microgrid[img, :])
        end
        
        push!(jpc_list, jpck)
    end
    
    # Check if there are buses not in any valid island (i.e., not in all_valid_island_buses and not in isolated)
    for i in 1:nb_ac
        bus_id = Int(jpc.busAC[i, BUS_I])
        if !(bus_id in all_valid_island_buses) && !(bus_id in isolated)
            push!(isolated, bus_id)
        end
    end
    
    # Detect if there are batteries in constant Vdc mode
    for i in 1:length(jpc_list)
        if !isempty(jpc_list[i].storage)
            # Use vectorized logical OR operation .||
            Vdc_converter = findall((jpc_list[i].converter[:,CONV_MODE] .== 4) .|| (jpc_list[i].converter[:,CONV_MODE] .== 5))
            if !isempty(Vdc_converter)
                @error "Batteries in constant Vdc mode detected in island $(i), which may cause power flow calculation errors. Please check JPC data."
            end
        end

        slack_bus_indices = findall(jpc_list[i].busAC[:, BUS_TYPE] .== REF)
        if length(slack_bus_indices) > 1
            # Check if there are multiple reference nodes
            @warn "Multiple reference nodes detected in island $(i), which may cause power flow calculation errors. Please check JPC data."
        end
    end

    return jpc_list, isolated
end

