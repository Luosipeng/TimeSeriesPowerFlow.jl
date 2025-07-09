"""
Write system summary section
"""
function write_system_summary(f::IOStream, mpc::JPC, area, isolated)
    # Extract necessary data
    baseMVA = mpc.baseMVA
    
    # Calculate basic statistics
    if hasproperty(mpc, :busAC)
        n_buses = size(mpc.busAC, 1)
    else
        # Infer number of buses from other data
        bus_ids = Set()
        if hasproperty(mpc, :branchAC)
            for i in 1:size(mpc.branchAC, 1)
                push!(bus_ids, mpc.branchAC[i, 1])
                push!(bus_ids, mpc.branchAC[i, 2])
            end
        end
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                push!(bus_ids, mpc.genAC[i, 1])
            end
        end
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                push!(bus_ids, mpc.loadAC[i, 1])
            end
        end
        n_buses = length(bus_ids)
    end
    n_isolated = length(isolated)
    n_buses = n_buses + n_isolated
    
    # Get number of generators
    n_gens = hasproperty(mpc, :genAC) ? size(mpc.genAC, 1) : 0
    
    # Get number of loads
    n_loads = hasproperty(mpc, :loadAC) ? size(mpc.loadAC, 1) : 0
    
    # Get number of branches
    n_branches = hasproperty(mpc, :branchAC) ? size(mpc.branchAC, 1) : 0
    
    # Calculate number of transformers (assuming branch matrix has tap column)
    n_transformers = 0
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 9
        for i in 1:size(mpc.branchAC, 1)
            if abs(mpc.branchAC[i, 9]) > 1e-6
                n_transformers += 1
            end
        end
    end
    
    # Calculate total generation and load
    total_gen_p = 0.0
    total_gen_q = 0.0
    if hasproperty(mpc, :genAC)
        # Assuming gen matrix column 2 is active power, column 3 is reactive power
        total_gen_p = sum(mpc.genAC[:, 2]) 
        total_gen_q = sum(mpc.genAC[:, 3]) 
    end
    
    total_load_p = 0.0
    total_load_q = 0.0
    if hasproperty(mpc, :busAC)
        # Assuming busAC matrix column 3 is active load, column 4 is reactive load
        if size(mpc.busAC, 2) >= 4
            total_load_p = sum(mpc.busAC[:, 3]) 
            total_load_q = sum(mpc.busAC[:, 4]) 
        end
    elseif hasproperty(mpc, :loadAC)
        # Get load data from loadAC
        if size(mpc.loadAC, 2) >= 4
            total_load_p = sum(mpc.loadAC[:, 3])
            total_load_q = sum(mpc.loadAC[:, 4])
        end
    end
    
    # Calculate losses
    total_p_loss = total_gen_p - total_load_p
    total_q_loss = 0.0
    charging_q = 0.0
    
    if hasproperty(mpc, :branchAC)
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
                
            # Get line parameters
            r = mpc.branchAC[i, 3] 
            x = mpc.branchAC[i, 4]
            b = mpc.branchAC[i, 5]  # Line half charging susceptance
                
            # Get actual bus voltage values and angles
            v_from = 1.0  # Default value if actual voltage not found
            v_to = 1.0    # Default value if actual voltage not found
            ang_from = 0.0
            ang_to = 0.0
                
            # Look up actual voltage values and angles from bus data
            if hasproperty(mpc, :busAC)
                for j in 1:size(mpc.busAC, 1)
                    if Int(mpc.busAC[j, 1]) == from_bus
                        v_from = mpc.busAC[j, 8]  # Use actual voltage magnitude
                        ang_from = mpc.busAC[j, 9] * pi/180  # Convert to radians
                    elseif Int(mpc.busAC[j, 1]) == to_bus
                        v_to = mpc.busAC[j, 8]    # Use actual voltage magnitude
                        ang_to = mpc.busAC[j, 9] * pi/180  # Convert to radians
                    end
                end
            end
                
            # Angle difference
            angle_diff = ang_from - ang_to
                
            # Calculate line admittance
            y = 1 / complex(r, x)
            y_abs = abs(y)
                
            # Directly calculate line current magnitude squared
            i_mag_squared = (v_from^2 + v_to^2 - 2*v_from*v_to*cos(angle_diff)) * y_abs^2
                
            # Calculate reactive losses - using reactance and current squared
            q_loss = x * i_mag_squared * baseMVA

            charging_from = 0.5 * b * v_from^2 * baseMVA
            charging_to = 0.5 * b * v_to^2 * baseMVA
            charging_q += charging_from + charging_to
            total_q_loss += q_loss
        end
    end
    
    # Calculate generator capacity
    total_gen_pmax = 0.0
    total_gen_qmin = 0.0
    total_gen_qmax = 0.0
    if hasproperty(mpc, :genAC) && size(mpc.genAC, 2) >= 9
        total_gen_pmax = sum(mpc.genAC[:, 9]) 
        total_gen_qmin = sum(mpc.genAC[:, 5])
        total_gen_qmax = sum(mpc.genAC[:, 4])
    end
    
    # Write system summary
    write(f, "================================================================================\n")
    write(f, "|     System Summary                                                           |\n")
    write(f, "================================================================================\n\n")
    
    @printf(f, "How many?                How much?              P (MW)            Q (MVAr)\n")
    @printf(f, "---------------------    -------------------  -------------  -----------------\n")
    @printf(f, "Buses            %3d     Total Gen Capacity    %7.1f       %7.1f to %7.1f\n", 
            n_buses, total_gen_pmax, total_gen_qmin, total_gen_qmax)
    @printf(f, "Generators       %3d     On-line Capacity      %7.1f       %7.1f to %7.1f\n", 
            n_gens, total_gen_pmax, total_gen_qmin, total_gen_qmax)
    @printf(f, "Committed Gens   %3d     Generation (actual)   %7.1f             %7.1f\n", 
            n_gens, total_gen_p, total_gen_q)
    @printf(f, "Loads            %3d     Load                  %7.1f            %7.1f\n", 
            n_loads, total_load_p, total_load_q)
    @printf(f, "  Fixed          %3d       Fixed               %7.1f            %7.1f\n", 
            n_loads, total_load_p, total_load_q)
    @printf(f, "  Dispatchable    %2d       Dispatchable          %4.1f of %4.1f      %5.1f\n", 
            0, 0.0, 0.0, 0.0)
    @printf(f, "Shunts           %3d     Shunt (inj)             %5.1f              %5.1f\n", 
            0, 0.0, 0.0)  # Assuming no shunt elements
    @printf(f, "Branches         %3d     Losses (I^2 * Z)       %6.2f            %6.2f\n", 
            n_branches, total_p_loss, total_q_loss)  # Reactive losses need to be calculated
    @printf(f, "Transformers     %3d     Branch Charging (inj)     -             %6.1f\n", 
            n_transformers, charging_q)  # Need to calculate actual values
    @printf(f, "Inter-ties        %2d     Total Inter-tie Flow     %4.1f               %4.1f\n", 
            0, 0.0, 0.0)
    @printf(f, "Areas             %2d\n\n", area)
    
    # Voltage and angle min/max values
    # This part needs to be extracted from bus data, if complete bus data is not available,
    # it can be omitted or estimated values can be used
    if hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 9
        min_vm = Inf
        max_vm = -Inf
        min_vm_bus = 0
        max_vm_bus = 0
        min_va = Inf
        max_va = -Inf
        min_va_bus = 0
        max_va_bus = 0
        
        for i in 1:size(mpc.busAC, 1)
            vm = mpc.busAC[i, 8]
            va = mpc.busAC[i, 9]
            
            if vm < min_vm
                min_vm = vm
                min_vm_bus = Int(mpc.busAC[i, 1])
            end
            if vm > max_vm
                max_vm = vm
                max_vm_bus = Int(mpc.busAC[i, 1])
            end
            
            if va < min_va
                min_va = va
                min_va_bus = Int(mpc.busAC[i, 1])
            end
            if va > max_va
                max_va = va
                max_va_bus = Int(mpc.busAC[i, 1])
            end
        end
        
        @printf(f, "                          Minimum                      Maximum\n")
        @printf(f, "                 -------------------------  --------------------------------\n")
        @printf(f, "Voltage Magnitude   %5.3f p.u. @ bus %3d         %5.3f p.u. @ bus %3d  \n", 
                min_vm, min_vm_bus, max_vm, max_vm_bus)
        @printf(f, "Voltage Angle      %6.2f deg   @ bus %3d        %6.2f deg   @ bus %3d  \n", 
                min_va, min_va_bus, max_va, max_va_bus)
    end
    
    # Line loss information
    # If detailed branch loss data is available, this part can be added
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 18
        # Process branch loss data
        max_p_loss= -Inf
        max_q_loss = -Inf
        max_p_line = 0
        max_q_line = 0
        p_loss = mpc.branchAC[:, 15] + mpc.branchAC[:, 17]
        q_loss = mpc.branchAC[:, 16] + mpc.branchAC[:, 18]
        for i in 1:size(mpc.branchAC, 1)
            ploss= p_loss[i]
            qloss= q_loss[i]
            if ploss > max_p_loss
                max_p_loss = ploss
                max_p_line = i
            end
            if qloss > max_q_loss
                max_q_loss = qloss
                max_q_line = i
            end
        end
        # Use estimated values or omit
        @printf(f, "P Losses (I^2*R)             -                  %5.2f MW    @ line %s\n", 
        max_p_loss, string(Int(mpc.branchAC[max_p_line, 1]), "-", Int(mpc.branchAC[max_p_line, 2])))
        @printf(f, "Q Losses (I^2*X)             -                 %5.2f MVAr  @ line %s\n", 
        max_q_loss, string(Int(mpc.branchAC[max_q_line, 1]), "-", Int(mpc.branchAC[max_q_line, 2])))
    else
        # Use estimated values or omit
        @printf(f, "P Losses (I^2*R)             -                  %5.2f MW    @ line %s\n", 
                0.0, "X-X")
        @printf(f, "Q Losses (I^2*X)             -                 %5.2f MVAr  @ line %s\n", 
                0.0, "X-X")
    end
    
    write(f, "\n")
end

"""
Write bus data section
"""
function write_bus_data(f::IOStream, mpc::JPC, isolated)
    baseMVA = mpc.baseMVA
    
    write(f, "================================================================================\n")
    write(f, "|     Bus Data                                                                 |\n")
    write(f, "================================================================================\n")
    write(f, " Bus      Voltage          Generation             Load        \n")
    write(f, "  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n")
    write(f, "----- ------- --------  --------  --------  --------  --------\n")
    
    # Assuming we have complete bus data
    if hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 9
        # Create generator and load lookup tables
        gen_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                bus_id = Int(mpc.genAC[i, 1])
                pg = mpc.genAC[i, 2] 
                qg = mpc.genAC[i, 3] 
                gen_lookup[bus_id] = (pg, qg)
            end
        end
        
        load_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                bus_id = Int(mpc.loadAC[i, 1])
                pd = mpc.loadAC[i, 3] 
                qd = mpc.loadAC[i, 4] 
                load_lookup[bus_id] = (pd, qd)
            end
        elseif hasproperty(mpc, :busAC) && size(mpc.busAC, 2) >= 4
            for i in 1:size(mpc.busAC, 1)
                bus_id = Int(mpc.busAC[i, 1])
                pd = mpc.busAC[i, 3] 
                qd = mpc.busAC[i, 4] 
                if pd != 0.0 || qd != 0.0
                    load_lookup[bus_id] = (pd, qd)
                end
            end
        end
        
        # Add isolated buses
        bus_data = copy(mpc.busAC)
        for i in eachindex(isolated)
            isolated_bus = zeros(1, size(bus_data, 2))
            isolated_bus[1] = isolated[i]
            isolated_bus[2] = 1.0  # Type (PQ node)
            isolated_bus[8] = 0.0  # Default voltage magnitude
            isolated_bus[9] = 0.0  # Default angle
            bus_data = vcat(bus_data, isolated_bus)
        end
        
        # Sort bus_data by bus number
        bus_ids = bus_data[:, 1]
        sorted_indices = sortperm(bus_ids)
        bus_data = bus_data[sorted_indices, :]

        # Totals
        total_pg = 0.0
        total_qg = 0.0
        total_pd = 0.0
        total_qd = 0.0
        
        # Iterate through all buses
        for i in 1:size(bus_data, 1)
            bus_id = Int(bus_data[i, 1])
            vm = bus_data[i, 8]
            va = bus_data[i, 9]
            
            # Get generation data
            pg_str = "-"
            qg_str = "-"
            if haskey(gen_lookup, bus_id)
                pg, qg = gen_lookup[bus_id]
                pg_str = @sprintf("%.2f", pg)
                qg_str = @sprintf("%.2f", qg)
                total_pg += pg
                total_qg += qg
            end
            
            # Get load data
            pd = 0.0
            qd = 0.0
            if haskey(load_lookup, bus_id)
                pd, qd = load_lookup[bus_id]
                total_pd += pd
                total_qd += qd
            end
            
            # Print bus data
            @printf(f, "%5d  %5.3f   %6.3f   %8s   %8s   %7.2f   %7.2f \n", 
                    bus_id, vm, va, pg_str, qg_str, pd, qd)
        end
        
        # Print totals
        @printf(f, "                        --------  --------  --------  --------\n")
        @printf(f, "               Total:   %7.2f   %7.2f   %7.2f   %7.2f\n", 
                total_pg, total_qg, total_pd, total_qd)
    else
        # If complete bus data is not available, try to construct from other data
        bus_ids = Set{Int}()
        
        # Collect all bus IDs from gen, load, and branch data
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                push!(bus_ids, Int(mpc.genAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                push!(bus_ids, Int(mpc.loadAC[i, 1]))
            end
        end
        
        if hasproperty(mpc, :branchAC)
            for i in 1:size(mpc.branchAC, 1)
                push!(bus_ids, Int(mpc.branchAC[i, 1]))
                push!(bus_ids, Int(mpc.branchAC[i, 2]))
            end
        end
        
        # Create generator and load lookup tables
        gen_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :genAC)
            for i in 1:size(mpc.genAC, 1)
                bus_id = Int(mpc.genAC[i, 1])
                pg = mpc.genAC[i, 2] * baseMVA
                qg = mpc.genAC[i, 3] * baseMVA
                gen_lookup[bus_id] = (pg, qg)
            end
        end
        
        load_lookup = Dict{Int, Tuple{Float64, Float64}}()
        if hasproperty(mpc, :loadAC)
            for i in 1:size(mpc.loadAC, 1)
                bus_id = Int(mpc.loadAC[i, 1])
                pd = mpc.loadAC[i, 3] * baseMVA
                qd = mpc.loadAC[i, 4] * baseMVA
                load_lookup[bus_id] = (pd, qd)
            end
        end
        
        # Totals
        total_pg = 0.0
        total_qg = 0.0
        total_pd = 0.0
        total_qd = 0.0
        
        # Iterate through all collected bus IDs
        for bus_id in sort(collect(bus_ids))
            # Assume voltage data
            vm = 1.0  # Default value
            va = 0.0  # Default value
            
            # Get generation data
            pg_str = "-"
            qg_str = "-"
            if haskey(gen_lookup, bus_id)
                pg, qg = gen_lookup[bus_id]
                pg_str = @sprintf("%.2f", pg)
                qg_str = @sprintf("%.2f", qg)
                total_pg += pg
                total_qg += qg
            end
            
            # Get load data
            pd = 0.0
            qd = 0.0
            if haskey(load_lookup, bus_id)
                pd, qd = load_lookup[bus_id]
                total_pd += pd
                total_qd += qd
            end
            
            # Print bus data
            @printf(f, "%5d  %5.3f   %6.3f   %8s   %8s   %7.2f   %7.2f \n", 
                    bus_id, vm, va, pg_str, qg_str, pd, qd)
        end
        
        # Print totals
        @printf(f, "                        --------  --------  --------  --------\n")
        @printf(f, "               Total:   %7.2f   %7.2f   %7.2f   %7.2f\n", 
                total_pg, total_qg, total_pd, total_qd)
    end
    
    write(f, "\n")
end

"""
Write branch data section
"""
function write_branch_data(f::IOStream, mpc::JPC)
    baseMVA = mpc.baseMVA
    
    write(f, "================================================================================\n")
    write(f, "|     Branch Data                                                              |\n")
    write(f, "================================================================================\n")
    write(f, "Brnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  \n")
    write(f, "  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n")
    write(f, "-----  -----  -----  --------  --------  --------  --------  --------  --------\n")
    
    # Check if branch data is available
    if !hasproperty(mpc, :branchAC) || size(mpc.branchAC, 1) == 0
        @printf(f, "No branch data available.\n")
        return
    end
    
    # Check if branch power flow data is available
    has_flow_data = false
    if hasproperty(mpc, :branchAC) && size(mpc.branchAC, 2) >= 18
        has_flow_data = true
    end
    
    total_p_loss = 0.0
    total_q_loss = 0.0
    
    if has_flow_data
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
            
            # Get power flow data
            pf = mpc.branchAC[i, 15] 
            qf = mpc.branchAC[i, 16] 
            pt = mpc.branchAC[i, 17] 
            qt = mpc.branchAC[i, 18] 
            
            # Calculate losses
            p_loss = pf + pt
            
            # Get line parameters
            r = mpc.branchAC[i, 3] 
            x = mpc.branchAC[i, 4]
            
            # Get actual bus voltage values and angles
            v_from = 1.0  # Default value if actual voltage not found
            v_to = 1.0    # Default value if actual voltage not found
            ang_from = 0.0
            ang_to = 0.0
            
            # Look up actual voltage values and angles from bus data
            if hasproperty(mpc, :busAC)
                for j in 1:size(mpc.busAC, 1)
                    if Int(mpc.busAC[j, 1]) == from_bus
                        v_from = mpc.busAC[j, 8]  # Use actual voltage magnitude
                        ang_from = mpc.busAC[j, 9] * pi/180  # Convert to radians
                    elseif Int(mpc.busAC[j, 1]) == to_bus
                        v_to = mpc.busAC[j, 8]    # Use actual voltage magnitude
                        ang_to = mpc.busAC[j, 9] * pi/180  # Convert to radians
                    end
                end
            end
            
            # Angle difference
            angle_diff = ang_from - ang_to
            
            # Calculate line admittance
            y = 1 / complex(r, x)
            y_abs = abs(y)
            
            # Directly calculate line current magnitude squared
            i_mag_squared = (v_from^2 + v_to^2 - 2*v_from*v_to*cos(angle_diff)) * y_abs^2
            
            # Calculate reactive losses - using reactance and current squared
            q_loss = x * i_mag_squared * baseMVA
            
            total_p_loss += p_loss
            total_q_loss += q_loss
            
            # Print branch data
            @printf(f, "%5d  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n", 
                    branch_id, from_bus, to_bus, pf, qf, pt, qt, p_loss, q_loss)
        end
    else
        # If power flow data is not available, only print branch topology information
        for i in 1:size(mpc.branchAC, 1)
            branch_id = i
            from_bus = Int(mpc.branchAC[i, 1])
            to_bus = Int(mpc.branchAC[i, 2])
            
            # Print branch data (without power flow information)
            @printf(f, "%5d  %5d  %5d  %8s  %8s  %8s  %8s  %8s  %8s\n", 
                    branch_id, from_bus, to_bus, "-", "-", "-", "-", "-", "-")
        end
    end
    
    # Print totals
    @printf(f, "                                                             --------  --------\n")
    @printf(f, "                                                    Total:   %8.2f  %8.2f\n", 
            total_p_loss, total_q_loss)
    
    write(f, "\n")
end

"""
Extract bus data from your data structure
"""
function extract_bus_data(mpc::JPC)
    # If bus data already exists, return it directly
    if hasproperty(mpc, :busAC)
        return mpc.busAC
    end
    
    # Otherwise, try to construct basic bus data from gen, load, and branch data
    bus_ids = Set{Int}()
    
    # Collect all bus IDs from gen, load, and branch data
    if hasproperty(mpc, :genAC)
        for i in 1:size(mpc.genAC, 1)
            push!(bus_ids, Int(mpc.genAC[i, 1]))
        end
    end
    
    if hasproperty(mpc, :loadAC)
        for i in 1:size(mpc.loadAC, 1)
            push!(bus_ids, Int(mpc.loadAC[i, 1]))
        end
    end
    
    if hasproperty(mpc, :branchAC)
        for i in 1:size(mpc.branchAC, 1)
            push!(bus_ids, Int(mpc.branchAC[i, 1]))
            push!(bus_ids, Int(mpc.branchAC[i, 2]))
        end
    end
    
    # Create basic bus data matrix
    # Columns: [bus_id, Vm, Va]
    bus_data = zeros(length(bus_ids), 3)
    
    for (i, bus_id) in enumerate(sort(collect(bus_ids)))
        bus_data[i, 1] = bus_id
        bus_data[i, 2] = 1.0  # Default voltage magnitude
        bus_data[i, 3] = 0.0  # Default angle
    end
    
    return bus_data
end

"""
Format power flow calculation results as MATPOWER-style report and save to a text file
"""
function generate_matpower_report(mpc::JPC, area, execution_time, isolated, output_file::String="PowerFlow_report.txt")
    # Open file for writing
    open(output_file, "w") do f
        # Write report header
        write(f, "JUPOWER Version 0.01, $(Dates.format(now(), "dd-u-yyyy"))\n")
        write(f, "Power Flow -- AC-polar-power formulation\n\n")
        
        # Write convergence information
        write(f, "Newton's method converged in $(mpc.iterationsAC) iterations.\n")
        if mpc.success
            write(f, "PF successful\n\n")
        else
            write(f, "PF NOT successful\n\n")
        end
        
        # Assumed calculation time
        write(f, "Converged in $(execution_time) seconds\n")
        
        # System summary
        write_system_summary(f, mpc, area, isolated)
        
        # Bus data
        write_bus_data(f, mpc, isolated)
        
        # Branch data
        write_branch_data(f, mpc)
    end
    
    println("Report saved to $output_file")
end

