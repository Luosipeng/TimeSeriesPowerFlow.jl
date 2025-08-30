"""
    runhpf(jpc, opt)

Run hybrid power flow calculation for integrated AC/DC systems.
Iteratively solves AC and DC power flow while updating converter power exchanges.
Supports different converter control modes.

Parameters:
- jpc: JPC object containing both AC and DC system data
- opt: Options for power flow calculation

Returns:
- Updated JPC object with power flow results
"""
function runhpf(jpc, opt)
    
    # Set iteration parameters
    max_iterations = 20  # Maximum number of iterations
    convergence_tolerance = 1e-4  # Convergence tolerance
    
    # Create a deep copy of the result
    result_jpc = deepcopy(jpc)
    
    # Save initial load data to avoid infinite accumulation
    initial_busAC = deepcopy(jpc.busAC)
    initial_busDC = deepcopy(jpc.busDC)
    initial_loadAC = deepcopy(jpc.loadAC)
    initial_loadDC = deepcopy(jpc.loadDC)
    
    # Save power values from previous iteration for convergence check
    prev_power_values = Dict()
    
    # Main iteration loop
    for iter in 1:max_iterations
        # Reset to initial load state to avoid cumulative addition
        result_jpc.busAC = deepcopy(initial_busAC)
        result_jpc.busDC = deepcopy(initial_busDC)
        result_jpc.loadAC = deepcopy(initial_loadAC)
        result_jpc.loadDC = deepcopy(initial_loadDC)
        
        # Prepare AC system data
        jpc1 = PowerFlow.JPC()
        jpc1.baseMVA = deepcopy(result_jpc.baseMVA)
        jpc1.busAC = deepcopy(result_jpc.busAC)
        jpc1.genAC = deepcopy(result_jpc.genAC)
        jpc1.loadAC = deepcopy(result_jpc.loadAC)
        jpc1.branchAC = deepcopy(result_jpc.branchAC)
        jpc1.pv_acsystem = deepcopy(result_jpc.pv_acsystem)
        jpc1.version = "2"
        
        # Prepare DC system data
        jpc2 = PowerFlow.JPC()
        jpc2.baseMVA = deepcopy(result_jpc.baseMVA)
        jpc2.busDC = deepcopy(result_jpc.busDC)
        jpc2.genDC = deepcopy(result_jpc.genDC)
        jpc2.loadDC = deepcopy(result_jpc.loadDC)
        jpc2.branchDC = deepcopy(result_jpc.branchDC)
        jpc2.pv = deepcopy(result_jpc.pv)
        jpc2.version = "2"
        
        # 2. Create power values dictionary for current iteration to check convergence
        current_power_values = Dict()
        
        # 3. Update power based on converter modes
        if !isempty(result_jpc.converter)
            # Process converters with mode constant δs, Us (CONV_MODE==1)
            PowerFlow.update_mode_1_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # Process converters with mode constant Ps, Qs (CONV_MODE==2)
            # Default mode, no special handling required
            
            # Process converters with mode constant Ps, Us (CONV_MODE==3)
            PowerFlow.update_mode_3_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # Process converters with mode constant Udc, Qs (CONV_MODE==4)
            PowerFlow.update_mode_4_converters!(result_jpc, jpc1, jpc2, current_power_values)
            
            # Process converters with mode constant Udc, Us (CONV_MODE==5)
            PowerFlow.update_mode_5_converters!(result_jpc, jpc1, jpc2, current_power_values)

            # Process converters with mode Droop Udc, Constant Qs (CONV_MODE==6)
            PowerFlow.update_mode_6_converters!(result_jpc, jpc1, jpc2, current_power_values)

            # Process converters with mode Droop Udc, Constant Us (CONV_MODE==7)
            PowerFlow.update_mode_7_converters!(result_jpc, jpc1, jpc2, current_power_values)

        end
        
        # 1. Calculate AC and DC power flow independently (after updating loads)
        if !isempty(jpc1.busAC)
            jpc1 = runpf(jpc1, opt)
        end
        
        if !isempty(jpc2.busDC)
            jpc2 = rundcpf(jpc2, opt)
        end
        
        # Check if power flow calculation succeeded
        if !isempty(jpc2.busDC)
            if !(jpc1.success && jpc2.success)
                @warn "Power flow calculation failed at iteration $iter"
                break
            end
        elseif !jpc1.success
            @warn "AC power flow calculation failed at iteration $iter"
            break
        end
        
        # 4. Update results
        result_jpc.busAC = deepcopy(jpc1.busAC)
        result_jpc.genAC = deepcopy(jpc1.genAC)
        result_jpc.branchAC = deepcopy(jpc1.branchAC)
        # result_jpc.loadAC = deepcopy(jpc1.loadAC)  # Keep commented as we want to preserve updated loads
        result_jpc.busDC = deepcopy(jpc2.busDC)
        result_jpc.genDC = deepcopy(jpc2.genDC)
        result_jpc.branchDC = deepcopy(jpc2.branchDC)
        result_jpc.loadDC = deepcopy(jpc2.loadDC)
        result_jpc.pv = deepcopy(jpc2.pv)
        result_jpc.pv_acsystem = deepcopy(jpc1.pv_acsystem)
        result_jpc.iterationsAC = deepcopy(jpc1.iterationsAC)
        result_jpc.iterationsDC = deepcopy(jpc2.iterationsDC)
        
        # 5. Check convergence
        if iter > 1
            max_diff = 0.0
            for (key, value) in current_power_values
                if haskey(prev_power_values, key)
                    diff = abs(value - prev_power_values[key])
                    max_diff = max(max_diff, diff)
                end
            end
            
            if max_diff < convergence_tolerance
                @info "Iteration converged at iteration $iter, maximum power difference: $max_diff"
                break
            elseif iter == max_iterations
                @warn "Reached maximum iterations $max_iterations without convergence. Maximum power difference: $max_diff"
            end
        end
        
        prev_power_values = deepcopy(current_power_values)
    end
    
    # Set solution status
    if !isempty(result_jpc.busDC)
        if result_jpc.iterationsAC >= 0 && result_jpc.iterationsDC >= 0
            result_jpc.success = true
        else
            result_jpc.success = false
        end
    else
        result_jpc.success = result_jpc.iterationsAC >= 0
    end
    
    return result_jpc
end

"""
    update_mode_1_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode constant δs, Us (CONV_MODE==1).
These converters have fixed AC voltage angle and magnitude, and the power is calculated from AC side.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_1_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_delta_us = filter(row -> row[CONV_MODE] == 1, eachrow(result_jpc.converter))
    if !isempty(converters_delta_us)
        converters_index = findall(row -> row[CONV_MODE] == 1, eachrow(result_jpc.converter))
        for i in eachindex(converters_delta_us[:, 1])
            conv = converters_delta_us[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            P_ac = -jpc1.genAC[gen_row, PG]
            Q_ac = -jpc1.genAC[gen_row, QG]
            
            if P_ac < 0
                P_dc = -P_ac/conv[CONV_EFF]
            else
                P_dc = -P_ac * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_Q_ac"] = Q_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_Q_AC] = Q_ac
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # Update DC side load
            update_dc_load!(result_jpc, jpc2, conv[CONV_DCBUS], P_dc)
        end
    end
end

"""
    update_mode_3_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode constant Ps, Us (CONV_MODE==3).
These converters have fixed active power and AC voltage magnitude.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_3_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_ps_us = filter(row -> row[CONV_MODE] == 3, eachrow(result_jpc.converter))
    if !isempty(converters_ps_us)
        converters_index = findall(row -> row[CONV_MODE] == 3, eachrow(result_jpc.converter))
        for i in eachindex(converters_ps_us[:, 1])
            conv = converters_ps_us[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            P_ac = -jpc1.genAC[gen_row, PG]
            Q_ac = -jpc1.genAC[gen_row, QG]
            
            if P_ac < 0
                P_dc = -P_ac/conv[CONV_EFF]
            else
                P_dc = -P_ac * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_Q_ac"] = Q_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_Q_AC] = Q_ac
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # Update DC side load
            update_dc_load!(result_jpc, jpc2, conv[CONV_DCBUS], P_dc)
        end
    end
end

"""
    update_mode_4_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode constant Udc, Qs (CONV_MODE==4).
These converters have fixed DC voltage and reactive power.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_4_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_udc_qs = filter(row -> row[CONV_MODE] == 4, eachrow(result_jpc.converter))
    if !isempty(converters_udc_qs)
        converters_index = findall(row -> row[CONV_MODE] == 4, eachrow(result_jpc.converter))
        for i in eachindex(converters_udc_qs[:, 1])
            conv = converters_udc_qs[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # Update AC side load
            update_ac_load!(result_jpc, jpc1, conv[CONV_ACBUS], P_ac)
        end
    end
end

"""
    update_mode_5_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode constant Udc, Us (CONV_MODE==5).
These converters have fixed DC voltage and AC voltage magnitude.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_5_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_udc_us = filter(row -> row[CONV_MODE] == 5, eachrow(result_jpc.converter))
    if !isempty(converters_udc_us)
        converters_index = findall(row -> row[CONV_MODE] == 5, eachrow(result_jpc.converter))
        for i in eachindex(converters_udc_us[:, 1])
            conv = converters_udc_us[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # Update AC side generator
            ac_gen_rows = findall(x -> x == Int(conv[CONV_ACBUS]), result_jpc.genAC[:, 1])
            ac_gen_rows_jpc1 = findall(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            if !isempty(ac_gen_rows)
                result_jpc.genAC[ac_gen_rows, PG] .= -P_ac
                jpc1.genAC[ac_gen_rows_jpc1, PG] .= -P_ac  # Update AC generator in jpc1
            end
        end
    end
end

"""
    update_mode_6_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode Droop Udc, Constant Qs (CONV_MODE==6).
These converters use DC voltage droop control and maintain constant reactive power.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_6_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_droop_udc_qs = filter(row -> row[CONV_MODE] == 6, eachrow(result_jpc.converter))
    if !isempty(converters_droop_udc_qs)
        converters_index = findall(row -> row[CONV_MODE] == 6, eachrow(result_jpc.converter))
        for i in eachindex(converters_droop_udc_qs[:, 1])
            conv = converters_droop_udc_qs[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            U_dc_limited = calculate_droop_voltage(-P_dc, conv[CONV_DROOP_KP],jpc2.genDC[gen_row,VG])
            jpc2.genDC[gen_row, VG] = U_dc_limited  # Update DC side voltage
            jpc2.busDC[findfirst(jpc2.busDC[:,BUS_I].==jpc2.genDC[gen_row,GEN_BUS]), VM] = U_dc_limited  # Update DC bus voltage
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac


            # Update AC side load
            update_ac_load!(result_jpc, jpc1, conv[CONV_ACBUS], P_ac)
        end
    end

end

"""
    update_mode_7_converters!(result_jpc, jpc1, jpc2, current_power_values)

Update converters with mode Droop Udc, Constant Us (CONV_MODE==7).
These converters use DC voltage droop control and maintain constant AC voltage magnitude.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- jpc2: JPC object for DC system
- current_power_values: Dictionary to store current power values for convergence check
"""
function update_mode_7_converters!(result_jpc, jpc1, jpc2, current_power_values)
    converters_droop_udc_us = filter(row -> row[CONV_MODE] == 7, eachrow(result_jpc.converter))
    if !isempty(converters_droop_udc_us)
        converters_index = findall(row -> row[CONV_MODE] == 7, eachrow(result_jpc.converter))
        for i in eachindex(converters_droop_udc_us[:, 1])
            conv = converters_droop_udc_us[i]  # Get current converter from row i
            gen_row = findfirst(x -> x == Int(conv[CONV_DCBUS]), jpc2.genDC[:, 1])
            P_dc = -jpc2.genDC[gen_row, PG]
            U_dc_limited = calculate_droop_voltage(-P_dc, conv[CONV_DROOP_KP],jpc2.genDC[gen_row,VG])
            jpc2.genDC[gen_row, VG] = U_dc_limited  # Update DC side voltage
            jpc2.busDC[findfirst(jpc2.busDC[:,BUS_I].==jpc2.genDC[gen_row,GEN_BUS]), VM] = U_dc_limited  # Update DC bus voltage
            
            if P_dc < 0
                P_ac = -P_dc/conv[CONV_EFF]
            else
                P_ac = -P_dc * conv[CONV_EFF]
            end

            # Save current power values for convergence check, using row index as converter ID
            conv_id = converters_index[i]
            current_power_values["conv_$(conv_id)_P_ac"] = P_ac
            current_power_values["conv_$(conv_id)_P_dc"] = P_dc

            # Update converter power
            result_jpc.converter[converters_index[i], CONV_P_DC] = P_dc
            result_jpc.converter[converters_index[i], CONV_P_AC] = P_ac

            # Update AC side generator
            ac_gen_rows = findall(x -> x == Int(conv[CONV_ACBUS]), result_jpc.genAC[:, 1])
            ac_gen_rows_jpc1 = findall(x -> x == Int(conv[CONV_ACBUS]), jpc1.genAC[:, 1])
            if !isempty(ac_gen_rows)
                result_jpc.genAC[ac_gen_rows, PG] .= -P_ac
                jpc1.genAC[ac_gen_rows_jpc1, PG] .= -P_ac  # Update AC generator in jpc1
            end
        end
    end
    
end

"""
    update_dc_load!(result_jpc, jpc2, dc_bus_id, P_dc)

Helper function to update DC side load with power from converter.
Creates a new load if none exists at the specified bus.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc2: JPC object for DC system
- dc_bus_id: ID of the DC bus to update
- P_dc: DC power to add to the load
"""
function update_dc_load!(result_jpc, jpc2, dc_bus_id, P_dc)
    dc_bus_id = Int(dc_bus_id)
    
    # Update DC bus power
    dc_bus_rows = findall(x -> x == dc_bus_id, result_jpc.busDC[:, 1])
    dc_bus_rows_jpc2 = findall(x -> x == dc_bus_id, jpc2.busDC[:, 1])
    if !isempty(dc_bus_rows)
        result_jpc.busDC[dc_bus_rows, PD] .+= P_dc
    end
    if !isempty(dc_bus_rows_jpc2)
        jpc2.busDC[dc_bus_rows_jpc2, PD] .+= P_dc
    end
    
    # Find corresponding DC load
    dc_load_rows = findall(x -> x == dc_bus_id, result_jpc.loadDC[:, 2])
    dc_load_rows_jpc2 = findall(x -> x == dc_bus_id, jpc2.loadDC[:, 2])
    
    if !isempty(dc_load_rows)
        # If load exists, update power
        result_jpc.loadDC[dc_load_rows, LOAD_PD] .+= P_dc
    else
        # If no DC load exists, create a virtual load
        new_load_id = isempty(result_jpc.loadDC) ? 1 : maximum(result_jpc.loadDC[:, 1]) + 1
        
        # Create new load row, initialize all values to 0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # Load ID
        new_load_row[LOAD_CND] = dc_bus_id  # Bus ID
        new_load_row[LOAD_STATUS] = 1       # Operating status
        new_load_row[LOAD_PD] = P_dc        # Active power
        new_load_row[LOADP_PERCENT] = 1.0   # Active power percentage
        
        # Add new load row to result_jpc's loadDC
        result_jpc.loadDC = vcat(result_jpc.loadDC, reshape(new_load_row, 1, :))
    end
    if !isempty(dc_load_rows_jpc2)
        # If load exists, update power
        jpc2.loadDC[dc_load_rows_jpc2, LOAD_PD] .+= P_dc
    else
        # If no DC load exists, create a virtual load
        new_load_id = isempty(jpc2.loadDC) ? 1 : maximum(jpc2.loadDC[:, 1]) + 1
        
        # Create new load row, initialize all values to 0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # Load ID
        new_load_row[LOAD_CND] = dc_bus_id  # Bus ID
        new_load_row[LOAD_STATUS] = 1       # Operating status
        new_load_row[LOAD_PD] = P_dc        # Active power
        new_load_row[LOADP_PERCENT] = 1.0   # Active power percentage
        
        # Add new load row to jpc2's loadDC
        jpc2.loadDC = vcat(jpc2.loadDC, reshape(new_load_row, 1, :))
    end
end

"""
    update_ac_load!(result_jpc, jpc1, ac_bus_id, P_ac)

Helper function to update AC side load with power from converter.
Creates a new load if none exists at the specified bus.

Parameters:
- result_jpc: Main JPC object containing all system data
- jpc1: JPC object for AC system
- ac_bus_id: ID of the AC bus to update
- P_ac: AC power to add to the load
"""
function update_ac_load!(result_jpc, jpc1, ac_bus_id, P_ac)
    ac_bus_id = Int(ac_bus_id)
    
    # Update AC bus power
    ac_bus_rows = findall(x -> x == ac_bus_id, result_jpc.busAC[:, 1])
    ac_bus_rows_jpc1 = findall(x -> x == ac_bus_id, jpc1.busAC[:, 1])
    if !isempty(ac_bus_rows)
        result_jpc.busAC[ac_bus_rows, PD] .+= P_ac
    end
    if !isempty(ac_bus_rows_jpc1)
        jpc1.busAC[ac_bus_rows_jpc1, PD] .+= P_ac
    end
    
    # Find corresponding AC load
    ac_load_rows = findall(x -> x == ac_bus_id, result_jpc.loadAC[:, 2])
    ac_load_rows_jpc1 = findall(x -> x == ac_bus_id, jpc1.loadAC[:, 2])
    
    if !isempty(ac_load_rows)
        # If load exists, update power
        result_jpc.loadAC[ac_load_rows, LOAD_PD] .+= P_ac
    else
        # If no AC load exists, create a virtual load
        new_load_id = isempty(result_jpc.loadAC) ? 1 : maximum(result_jpc.loadAC[:, 1]) + 1
        
        # Create new load row, initialize all values to 0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # Load ID
        new_load_row[LOAD_CND] = ac_bus_id  # Bus ID
        new_load_row[LOAD_STATUS] = 1       # Operating status
        new_load_row[LOAD_PD] = P_ac        # Active power
        new_load_row[LOADP_PERCENT] = 1.0   # Active power percentage
        
                # Add new load row to result_jpc's loadAC
        result_jpc.loadAC = vcat(result_jpc.loadAC, reshape(new_load_row, 1, :))
    end
    if !isempty(ac_load_rows_jpc1)
        # If load exists, update power
        jpc1.loadAC[ac_load_rows_jpc1, LOAD_PD] .+= P_ac
    else
        # If no AC load exists, create a virtual load
        new_load_id = isempty(jpc1.loadAC) ? 1 : maximum(jpc1.loadAC[:, 1]) + 1
        
        # Create new load row, initialize all values to 0
        new_load_row = zeros(1, 8)
        new_load_row[LOAD_I] = new_load_id  # Load ID
        new_load_row[LOAD_CND] = ac_bus_id  # Bus ID
        new_load_row[LOAD_STATUS] = 1       # Operating status
        new_load_row[LOAD_PD] = P_ac        # Active power
        new_load_row[LOADP_PERCENT] = 1.0   # Active power percentage
        
        # Add new load row to jpc1's loadAC
        jpc1.loadAC = vcat(jpc1.loadAC, reshape(new_load_row, 1, :))
    end
end


"""
    calculate_droop_voltage(P_dc, k_p, U_dc_ref=1.0, U_dc_min=0.95, U_dc_max=1.05)

Calculate DC voltage with droop control using the equation: U_dc = U_dc_ref - k_p * P_dc
Limits the output voltage to be within specified minimum and maximum values.

Parameters:
- P_dc: DC power
- k_p: Droop coefficient
- U_dc_ref: Reference DC voltage (default: 1.0 p.u.)
- U_dc_min: Minimum allowed DC voltage (default: 0.95 p.u.)
- U_dc_max: Maximum allowed DC voltage (default: 1.05 p.u.)

Returns:
- Limited DC voltage value
"""
function calculate_droop_voltage(P_dc, k_p, U_dc_ref=1.0, U_dc_min=0.95, U_dc_max=1.05)
    # Calculate DC voltage with droop control
    U_dc_droop = U_dc_ref - k_p * P_dc
    
    # Voltage limit
    U_dc_limited = max(U_dc_min, min(U_dc_max, U_dc_droop))
    
    return U_dc_limited
end

