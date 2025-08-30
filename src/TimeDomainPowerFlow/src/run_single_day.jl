"""
    run_single_day(old_jpc, opt, day_load_matrix, day_price_line, day_irradiance_line, day_storage_line = nothing)

Run a day-ahead simulation for a hybrid AC-DC power system with renewable generation and energy storage.

# Arguments
- `old_jpc`: Original power system data structure
- `opt`: Power flow options
- `day_load_matrix`: Matrix containing hourly load data (hour, bus_id, active_power, reactive_power)
- `day_price_line`: Vector containing hourly electricity prices
- `day_irradiance_line`: Matrix containing hourly solar irradiance data (hour, irradiance)

# Returns
- `results_pf`: Array of power flow results for each hour of the day

# Description
This function performs a day-ahead economic dispatch and power flow analysis for a hybrid AC-DC power system
with renewable generation and energy storage. The process involves:

1. Preprocessing the power system data:
   - Adjusting bus loads to account for converter power flows
   - Renumbering the hybrid system for consistent indexing
   - Extracting load, generator, converter, PV system, and energy storage parameters

2. Running the dynamic economic dispatch optimization:
   - Minimizing generation costs while satisfying system constraints
   - Determining optimal converter operation, generator outputs, and storage charging/discharging

3. Performing hourly power flow calculations:
   - Updating loads, generators, converters, and PV outputs based on optimization results
   - Running hybrid power flow for each hour of the day
   - Handling special converter control modes (e.g., droop control)

The function handles various components including AC and DC loads, generators, converters between AC and DC
subsystems, PV systems (both AC-connected and DC-connected), and energy storage systems. It accounts for
time-varying load profiles, electricity prices, and solar irradiance patterns throughout the day.
"""

function run_single_day(old_jpc, opt, day_load_matrix, day_price_line, day_irradiance_line, day_storage_line = nothing)

    jpc = deepcopy(old_jpc)

    if day_storage_line == nothing
        converter_ac_bus = Int.(jpc.converter[:, CONV_ACBUS])
        converter_ac_p_mw = jpc.converter[:, CONV_P_AC]
        converter_ac_q_mvar = jpc.converter[:, CONV_Q_AC]
        converter_dc_bus = Int.(jpc.converter[:, CONV_DCBUS])
        converter_dc_p_mw = jpc.converter[:, CONV_P_DC]

        # Create mapping from bus ID to index
        ac_bus_map = Dict(jpc.busAC[i, BUS_I] => i for i in 1:size(jpc.busAC, 1))
        dc_bus_map = Dict(jpc.busDC[i, BUS_I] => i for i in 1:size(jpc.busDC, 1))

        # AC side processing - subtract converter power
        for i in 1:length(converter_ac_bus)
            ac_bus_id = converter_ac_bus[i]
            if haskey(ac_bus_map, ac_bus_id)
                idx = ac_bus_map[ac_bus_id]
                if jpc.converter[i, CONV_MODE] != 6 && jpc.converter[i, CONV_MODE] != 7
                jpc.busAC[idx, PD] -= converter_ac_p_mw[i]
                end
                if jpc.converter[i, CONV_MODE] != 7
                    jpc.busAC[idx, QD] -= converter_ac_q_mvar[i]
                end
            else
                @warn "Cannot find AC bus with ID $ac_bus_id"
            end
        end

        # DC side processing - subtract converter power
        for i in 1:length(converter_dc_bus)
            dc_bus_id = converter_dc_bus[i]
            if haskey(dc_bus_map, dc_bus_id)
                idx = dc_bus_map[dc_bus_id]
                if jpc.converter[i, CONV_MODE] != 6 && jpc.converter[i, CONV_MODE] != 7
                    jpc.busDC[idx, PD] -= converter_dc_p_mw[i]
                end
            else
                @warn "Cannot find DC bus with ID $dc_bus_id"
            end
        end

        # Delete or subtract the corresponding load
        for i in 1:size(jpc.loadAC, 1)
            load_bus = jpc.loadAC[i, LOAD_CND]
            for j in 1:length(converter_ac_bus)
                if load_bus == converter_ac_bus[j]
                    # Only subtract the converter power part, not delete the entire load
                    if jpc.converter[j, CONV_MODE] != 6 && jpc.converter[j, CONV_MODE] != 7
                        jpc.loadAC[i, LOAD_PD] -= converter_ac_p_mw[j]
                    end
                    if jpc.converter[j, CONV_MODE] != 7
                        # Only subtract the converter reactive power part, not delete the entire load
                        jpc.loadAC[i, LOAD_QD] -= converter_ac_q_mvar[j]
                    end
                end
            end
        end

        for i in 1:size(jpc.loadDC, 1)
            load_bus = jpc.loadDC[i, LOAD_CND]
            for j in 1:length(converter_dc_bus)
                if load_bus == converter_dc_bus[j]
                    if jpc.converter[j, CONV_MODE] != 6 && jpc.converter[j, CONV_MODE] != 7
                        # Only subtract the converter power part, not delete the entire load
                        jpc.loadDC[i, LOAD_PD] -= converter_dc_p_mw[j]
                    end
                end
            end
        end
    end

   new_jpc, ac_node_mapping, dc_node_mapping = TimeDomainPowerFlow.renumber_hybrid_system(jpc)
    

    ac_node_reverse_mapping = Dict{Int, Float64}()
    dc_node_reverse_mapping = Dict{Int, Float64}()
    
    # Fill reverse mapping dictionaries
    for (old_num, new_num) in ac_node_mapping
        ac_node_reverse_mapping[new_num] = old_num
    end
    
    for (old_num, new_num) in dc_node_mapping
        dc_node_reverse_mapping[new_num] = old_num
    end

    # Use the new JPC object for subsequent operations
    # Extract load information
    loadAC = new_jpc.loadAC
    loadDC = new_jpc.loadDC

    inserviced_line_AC = findall(loadAC[:, LOAD_STATUS] .== 1)
    inserviced_line_DC = findall(loadDC[:, LOAD_STATUS] .== 1)

    # Filter out unused loads
    loadAC = loadAC[inserviced_line_AC, :]
    loadDC = loadDC[inserviced_line_DC, :]

    num_load_ac = size(loadAC, 1)
    num_load_dc = size(loadDC, 1)

    # Extract load data
    loadAC_PD = zeros(num_load_ac, 24)
    loadAC_QD = zeros(num_load_ac, 24)
    # loadDC_PD = zeros(num_load_dc*24, 2)
    
    for i in 1:24
        # Find all data for the current hour
        idx = findall(day_load_matrix[:, 1] .== i)
        # Assign values for AC loads
        loadAC_PD[:, i] = day_load_matrix[idx, 3]
        loadAC_QD[:, i] = day_load_matrix[idx, 4]
        
    end
    
    loadDC_PD = zeros(num_load_dc, 24)
    for i in 1:24
        # Find all data for the current hour
        # idx = findall(day_load_matrix[:, 1] .== i)
        
        # Assign values for DC loads
        loadDC_PD[:, i] = loadDC[:, LOAD_PD]
    end
    # TODO: day_load_matrix needs to be adjusted according to actual conditions
    # loadDC_PD = loadDC[:, LOAD_PD]

    # Extract load location information - using new numbering
    ac_bus_ids = loadAC[:, LOAD_CND]
    dc_bus_ids = loadDC[:, LOAD_CND]

    # Build load connection matrix - using total number of nodes with new numbering
    total_buses = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)
    
    # AC load connection matrix
    Cld_ac = zeros(total_buses, size(loadAC, 1))
    for i in eachindex(ac_bus_ids)
        bus_id = Int(ac_bus_ids[i])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cld_ac[bus_id, i] = 1
        else
            @warn "Cannot find AC bus with ID $bus_id"
        end
    end
    
    # DC load connection matrix
    Cld_dc = zeros(total_buses, size(loadDC, 1))
    for i in eachindex(dc_bus_ids)
        bus_id = Int(dc_bus_ids[i])
        if bus_id in new_jpc.busDC[:, BUS_I]
            Cld_dc[bus_id, i] = 1
        else
            @warn "Cannot find DC bus with ID $bus_id"
        end
    end
    
    # Extract generator information - using new JPC
    genAC = new_jpc.genAC
    genDC = new_jpc.genDC
    genAC = genAC[genAC[:, GEN_STATUS] .== 1, :]  # Only keep generators with status 1
    genDC = genDC[genDC[:, GEN_STATUS] .== 1, :]  # Only keep generators with status 1

    # Find slack bus - using new numbering
    # slack_busAC = new_jpc.busAC[findfirst(new_jpc.busAC[:, BUS_TYPE] .== REF), BUS_I]
    # keep_indices = genAC[:, GEN_BUS] .!= slack_busAC
    # genAC = genAC[keep_indices, :]

    # Find converter nodes - using new numbering
    if !isempty(new_jpc.converter)
        conv_ac_bus = new_jpc.converter[:, CONV_ACBUS]
        keep_indices = genAC[:, GEN_BUS] .!= conv_ac_bus
        genAC = genAC[keep_indices, :]
    end
    
    # Extract generator power
    genAC_PG = genAC[:, PG]
    # Build generator connection matrix - using new numbering
    Cgen_ac = zeros(total_buses, size(genAC, 1))
    for i in eachindex(genAC[:, GEN_BUS])
        bus_id = Int(genAC[i, GEN_BUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cgen_ac[bus_id, i] = 1
        else
            @warn "Cannot find AC bus with ID $bus_id"
        end
    end

    # Extract converter information - using new JPC
    converters = new_jpc.converter
    Pij_inv = zeros(size(converters, 1))
    Pij_rec = zeros(size(converters, 1))
    Qij_inv = zeros(size(converters, 1))
    η_rec = converters[:,CONV_EFF]
    η_inv = converters[:, CONV_EFF]
    
    for i in eachindex(converters[:, 1])
        Pac = converters[i, CONV_P_AC]
        Pdc = converters[i, CONV_P_DC]
        Qac = converters[i, CONV_Q_AC]

        if Pac <= 0
            Pij_inv[i] = -Pac
        else
            Pij_inv[i] = 0
        end

        if Pdc <= 0
            Pij_rec[i] = -Pdc
        else
            Pij_rec[i] = 0
        end

        if Qac < 0
            Qij_inv[i] = -Qac
        else
            Qij_inv[i] = 0
        end
    end

    # Build converter connection matrix - using new numbering
    Cconv_ac = zeros(total_buses, size(converters, 1))
    for i in eachindex(converters[:, CONV_ACBUS])
        bus_id = Int(converters[i, CONV_ACBUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cconv_ac[bus_id, i] = 1
        else
            @warn "Cannot find AC bus with ID $bus_id"
        end
    end

    Cconv_dc = zeros(total_buses, size(converters, 1))
    for i in eachindex(converters[:, CONV_DCBUS])
        bus_id = Int(converters[i, CONV_DCBUS])
        if bus_id in new_jpc.busDC[:, BUS_I]

            Cconv_dc[ bus_id, i] = 1
        else
            @warn "Cannot find DC bus with ID $bus_id"
        end
    end

    # Extract PV AC system information - using new JPC
    pv_ac = new_jpc.pv_acsystem

    # Extract base irradiance for PV AC systems
    pv_ac_base_irradiance = pv_ac[:, PV_AC_IRRADIANCE]
    # Extract PV output upper limit
    pv_ac_p_mw_ratio = zeros(size(pv_ac, 1), 24)
    pv_ac_p_mw_ratio = inv.(pv_ac_base_irradiance) * day_irradiance_line[:, 2]'
    pv_ac_p_mw = pv_ac[:, PV_AC_INVERTER_PAC] 

    pv_ac_inverter_pac = pv_ac[:, PV_AC_INVERTER_PAC]

    # Extract PV AC connection matrix - using new numbering
    Cpv_ac = zeros(total_buses, size(pv_ac, 1))
    for i in eachindex(pv_ac[:, PV_AC_BUS])
        bus_id = Int(pv_ac[i, PV_AC_BUS])
        if bus_id in new_jpc.busAC[:, BUS_I]
            Cpv_ac[bus_id, i] = 1
        else
            @warn "Cannot find AC bus with ID $bus_id"
        end
    end

    # Extract PV information - using new JPC
    pvarray = new_jpc.pv

    # Extract base irradiance for PV
    pv_base_irradiance = pvarray[:, PV_IRRADIANCE]

    # Extract PV output upper limit
    # pv_max_p_mw_ratio = day_irradiance_line[:, 2] ./ pv_base_irradiance
    pv_max_p_mw_ratio = inv.(pv_base_irradiance)* day_irradiance_line[:, 2]'
    pv_max_p_mw = pvarray[:, PV_VMPP].* pvarray[:, PV_IMPP]./1000000

    pv_isc = pvarray[:, PV_ISC]
    pv_impp = pvarray[:, PV_IMPP]
    
    # Extract PV connection matrix - using new numbering
    Cpv_dc = zeros(total_buses, size(pvarray, 1))
    for i in eachindex(pvarray[:, PV_BUS])
        bus_id = Int(pvarray[i, PV_BUS])
        if bus_id in new_jpc.busDC[:, BUS_I]
            Cpv_dc[bus_id, i] = 1
        else
            @warn "Cannot find DC bus with ID $bus_id"
        end
    end

    # Extract storage information - using new JPC
    storage = new_jpc.storage
    # Extract storage output information
    ess_power_capacity_mw = storage[:, ESS_POWER_CAPACITY]
    ess_energy_capacity_mwh = storage[:, ESS_ENERGY_CAPACITY]
    ess_initial_soc = storage[:, ESS_SOC_INIT]

    # Extract storage SOC information
    ess_min_soc = storage[:, ESS_SOC_MIN]
    ess_max_soc = storage[:, ESS_SOC_MAX]
    ess_efficiency = storage[:, ESS_EFFICIENCY]

    # Extract storage connection matrix - using new numbering
    # Note: The dimension in the original code seems problematic, should be total_buses instead of just busDC rows
    Cstorage_ac = zeros(total_buses, size(storage, 1))
    for i in eachindex(storage[:, ESS_BUS])
        bus_id = Int(storage[i, ESS_BUS])
        if bus_id in new_jpc.busDC[:, BUS_I]  # May need to check if this should be busAC
            Cstorage_ac[bus_id, i] = 1
        else
            @warn "Cannot find DC bus with ID $bus_id"
        end
    end

    if day_storage_line == nothing
        result = run_dynamic_dispatch(new_jpc,
            Cld_ac, Cld_dc, 
            loadAC_PD, loadAC_QD,
            loadDC_PD,
            genAC_PG,
            Cgen_ac, 
            Cconv_ac, Cconv_dc, 
            η_rec, η_inv,
            Cpv_ac,
            Cpv_dc,
            pv_ac_p_mw_ratio,
            pv_ac_p_mw,
            pv_max_p_mw,
            pv_max_p_mw_ratio,
            Cstorage_ac,
            ess_initial_soc,
            ess_max_soc,
            ess_min_soc,
            ess_power_capacity_mw,
            ess_energy_capacity_mwh,
            ess_efficiency,
            day_price_line
            )
    end
    # Full day power flow calculation
    results_pf = Array{Any}(undef, 24)
    # Assign values
    for t in 1:24
        # Load amplitudes
        loadAC_PD_t = loadAC_PD[:, t]
        loadAC_QD_t = loadAC_QD[:, t]
        loadDC_PD_t = loadDC_PD[:, t]

        # Update loads
        for i in eachindex(jpc.loadAC[:, LOAD_CND])
            bus_id = jpc.loadAC[i, LOAD_CND]
            if bus_id in jpc.busAC[:, BUS_I]
                bus_index = findfirst(x -> x == bus_id, jpc.busAC[:, BUS_I])
                jpc.busAC[bus_index, PD] = loadAC_PD_t[i]
                jpc.busAC[bus_index, QD] = loadAC_QD_t[i]

                # Update load
                idx = findfirst(jpc.loadAC[:, LOAD_CND] .== bus_id)
                if idx !== nothing
                    jpc.loadAC[idx, LOAD_PD] = loadAC_PD_t[i]
                    jpc.loadAC[idx, LOAD_QD] = loadAC_QD_t[i]
                else
                    @warn "Cannot find AC bus load with ID $bus_id"
                end
            else
                @warn "Cannot find AC bus with ID $bus_id"
            end
        end

        for i in eachindex(jpc.loadDC[:, LOAD_CND])
            bus_id = jpc.loadDC[i, LOAD_CND]
            if bus_id in jpc.busDC[:, BUS_I]
                bus_index = findfirst(x -> x == bus_id, jpc.busDC[:, BUS_I])
                jpc.busDC[bus_index, PD] = loadDC_PD_t[i]

                # Update load
                idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_id)
                if idx !== nothing
                    jpc.loadDC[idx, LOAD_PD] = loadDC_PD_t[i]
                else
                    @warn "Cannot find DC bus load with ID $bus_id"
                end
            else
                @warn "Cannot find DC bus with ID $bus_id"
            end
        end

        if day_storage_line == nothing
            # Reassign results to jpc.converter
            for i in eachindex(converters[:, 1])
                if result["Pij_inv"][i,t] > 0
                    jpc.converter[i, CONV_P_AC] = result["Pij_inv"][i,t]
                else
                    jpc.converter[i, CONV_P_AC] = -result["Pij_rec"][i,t] * jpc.converter[i, CONV_EFF]
                end

                if result["Pij_rec"][i,t] > 0
                    jpc.converter[i, CONV_P_DC] = result["Pij_rec"][i,t]
                else
                    jpc.converter[i, CONV_P_DC] = -result["Pij_inv"][i,t] * jpc.converter[i, CONV_EFF]
                end
                # jpc.converter[i, CONV_Q_AC] = Qij_inv_result[i]
            end

            # Modify load or generator power according to control mode
            for i in eachindex(converters[:, 1])
                if jpc.converter[i, CONV_MODE] == 6
                    # If in Droop_Udc_Qs mode, directly set AC power to 0
                    genDC_idx = findfirst(jpc.genDC[:, GEN_BUS] .== jpc.converter[i, CONV_DCBUS])
                    jpc.genDC[genDC_idx, PG] = -jpc.converter[i, CONV_P_DC]
                    # Step 1: Find the corresponding index
                    idx = findfirst(jpc.busAC[:, BUS_I] .== jpc.converter[i, CONV_ACBUS])
                    
                    # Step 2: If index found, update the active power of the bus
                    # if idx !== nothing
                    #     jpc.busAC[idx, PD] += jpc.converter[i, CONV_P_AC]
                    # end
                    # # Step 1: Find the corresponding index
                    # idx = findfirst(jpc.loadAC[:, LOAD_CND] .== jpc.converter[i, CONV_ACBUS])

                    # # Step 2: If index found, update the active power of the load
                    # if idx !== nothing
                    #     jpc.loadAC[idx, LOAD_PD] += jpc.converter[i, CONV_P_AC]
                    # end

                elseif jpc.converter[i, CONV_MODE] == 7
                    # If in Droop_Udc_Us mode, directly set reactive power to 0
                    genAC_idx = findfirst(jpc.genAC[:, GEN_BUS] .== jpc.converter[i, CONV_ACBUS])
                    jpc.genAC[genAC_idx, PG] = -jpc.converter[i, CONV_P_AC]

                    genDC_idx = findfirst(jpc.genDC[:, GEN_BUS] .== jpc.converter[i, CONV_DCBUS])
                    jpc.genDC[genDC_idx, PG] = -jpc.converter[i, CONV_P_DC]

                end
            end
        

            # Assign values for storage loads
            for i in eachindex(storage[:, ESS_BUS])
                bus_id = storage[i, ESS_BUS]
                if bus_id in new_jpc.busDC[:, BUS_I]
                    bus_id = dc_node_reverse_mapping[bus_id]  # Use new DC bus ID
                    # ac_bus_count = size(new_jpc.busAC, 1)
                    line = findfirst(jpc.busDC[:,BUS_I] .== bus_id)
                    jpc.busDC[line, PD] = result["ess_charge"][i,t]- result["ess_discharge"][i,t]

                    # Update storage load
                    idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_id)
                    if idx !== nothing
                        jpc.loadDC[idx, LOAD_PD] = result["ess_charge"][i,t]- result["ess_discharge"][i,t]
                    else
                        @warn "Cannot find DC bus load with ID $bus_id"
                    end
                else
                    @warn "Cannot find DC bus with ID $bus_id"
                end
            end

        else
            # Assign values for storage loads
            for i in eachindex(storage[:, ESS_BUS])
                bus_id = storage[i, ESS_BUS]
                if bus_id in new_jpc.busDC[:, BUS_I]
                    bus_id = dc_node_reverse_mapping[bus_id]  # Use new DC bus ID
                    # ac_bus_count = size(new_jpc.busAC, 1)
                    line = findfirst(jpc.busDC[:,BUS_I] .== bus_id)
                    jpc.busDC[line, PD] = day_storage_line[t, 2]

                    # Update storage load
                    idx = findfirst(jpc.loadDC[:, LOAD_CND] .== bus_id)
                    if idx !== nothing
                        jpc.loadDC[idx, LOAD_PD] =day_storage_line[t, 2]
                    else
                        @warn "Cannot find DC bus load with ID $bus_id"
                    end
                else
                    @warn "Cannot find DC bus with ID $bus_id"
                end
            end
        end

        # Assign values for PV based on irradiance
        jpc.pv[:, PV_ISC] = pv_isc.*pv_max_p_mw_ratio[:,t]
        jpc.pv[:, PV_IMPP] = pv_impp.*pv_max_p_mw_ratio[:,t]

        # Update PV AC system power
        jpc.pv_acsystem[:, PV_AC_INVERTER_PAC] = pv_ac_inverter_pac.*pv_ac_p_mw_ratio[:, t]

        results_pf[t] = PowerFlow.runhpf(jpc, opt)
    end

    return results_pf
end
