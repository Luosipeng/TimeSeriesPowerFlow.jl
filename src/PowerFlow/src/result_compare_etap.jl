"""
    compare_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String; 
                           tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3)

Compare power flow calculation results in JPC format with reference voltage values, using JuliaPowerCase for node mapping.
Supports both original format and new format with mixed AC/DC data.
Parameters:
- results: NamedTuple containing results in JPC format (typically the return value of power flow calculation)
- case: Original JuliaPowerCase object, used to get node names and ID mappings
- reference_file: Path to Excel file containing reference voltage values
- tolerance_mag: Tolerance for voltage magnitude comparison (default: 1e-4)
- tolerance_ang: Tolerance for voltage angle comparison (default: 1e-3)

Returns:
- NamedTuple containing comparison results, including ac and dc DataFrames
"""
function compare_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String; 
                                tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3)
    # Check if reference file exists
    if !isfile(reference_file)
        error("Reference file $reference_file does not exist")
    end
    
    # Read reference file
    reference_data = DataFrame(XLSX.readtable(reference_file, "Sheet1"))
    
    # Check reference file format, determine if it contains type column (AC/DC mixed format)
    has_type_column = "type" in lowercase.(names(reference_data))
    
    # Create mapping from AC bus ID to name
    ac_bus_id_to_name = Dict{Int, String}()
    ac_bus_name_to_id = Dict{String, Int}()
    for bus in case.busesAC
        if bus.in_service
            ac_bus_id_to_name[bus.bus_id] = bus.name
            ac_bus_name_to_id[bus.name] = bus.bus_id
        end
    end
    
    # Create mapping from DC bus ID to name (if case has DC buses)
    dc_bus_id_to_name = Dict{Int, String}()
    dc_bus_name_to_id = Dict{String, Int}()
    if hasproperty(case, :busesDC)
        for bus in case.busesDC
            if bus.in_service
                dc_bus_id_to_name[bus.bus_id] = bus.name
                dc_bus_name_to_id[bus.name] = bus.bus_id
            end
        end
    end
    
    # Create AC results comparison dataframe
    ac_comparison_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Calc_Volt_Mag = Float64[],
        Ref_Volt_Mag = Float64[],
        Mag_Diff = Float64[],
        Mag_Error_Percent = Float64[],
        Mag_Within_Tolerance = Bool[],
        Calc_Volt_Ang = Float64[],
        Ref_Volt_Ang = Float64[],
        Ang_Diff = Float64[],
        Ang_Within_Tolerance = Bool[]
    )
    
    # Create DC results comparison dataframe
    dc_comparison_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Calc_Volt_Mag = Float64[],
        Ref_Volt_Mag = Float64[],
        Mag_Diff = Float64[],
        Mag_Error_Percent = Float64[],
        Mag_Within_Tolerance = Bool[]
    )
    
    # Get voltage values from calculation results
    voltage_results = get_bus_voltage_results_acdc(results, case)
    ac_results = voltage_results.ac
    dc_results = voltage_results.dc
    
    # Create mapping from AC bus name to results
    ac_name_to_result = Dict{String, NamedTuple}()
    for row in eachrow(ac_results)
        ac_name_to_result[row.Bus_Name] = (vm = row.Volt_Mag, va = row.Volt_Ang)
    end
    
    # Create mapping from DC bus name to results
    dc_name_to_result = Dict{String, Float64}()
    for row in eachrow(dc_results)
        dc_name_to_result[row.Bus_Name] = row.Volt_Mag
    end
    
    # Iterate through each bus in the reference file
    for i in 1:size(reference_data, 1)
        bus_name = reference_data.bus_ID[i]
        
        # Determine bus type
        if has_type_column
            # New format: use type column to determine bus type
            bus_type = uppercase(reference_data.type[i])
        else
            # Original format: assume all buses are AC
            bus_type = "AC"
        end
        
        ref_mag = reference_data.volt_mag[i]
        
        if bus_type == "AC"
            # Process AC bus
            ref_ang = reference_data.volt_ang[i]
            
            # Find corresponding bus in calculation results
            if haskey(ac_name_to_result, bus_name)
                result = ac_name_to_result[bus_name]
                vm = result.vm
                va = result.va
                
                # Calculate differences
                mag_diff = abs(vm - ref_mag)
                mag_error_percent = ref_mag != 0 ? (mag_diff / ref_mag) * 100 : 0.0
                mag_within_tol = mag_diff <= tolerance_mag
                
                ang_diff = abs(va - ref_ang)
                ang_within_tol = ang_diff <= tolerance_ang
                
                # Get bus ID
                bus_id = get(ac_bus_name_to_id, bus_name, -1)
                
                push!(ac_comparison_df, [
                    bus_id,
                    bus_name,
                    vm,
                    ref_mag,
                    mag_diff,
                    mag_error_percent,
                    mag_within_tol,
                    va,
                    ref_ang,
                    ang_diff,
                    ang_within_tol
                ])
            else
                @warn "AC bus $bus_name not found in calculation results, skipping comparison"
            end
        elseif bus_type == "DC"
            # Process DC bus
            # Find corresponding bus in calculation results
            if haskey(dc_name_to_result, bus_name)
                vm = dc_name_to_result[bus_name]
                
                # Calculate differences
                mag_diff = abs(vm - ref_mag)
                mag_error_percent = ref_mag != 0 ? (mag_diff / ref_mag) * 100 : 0.0
                mag_within_tol = mag_diff <= tolerance_mag
                
                # Get bus ID
                bus_id = get(dc_bus_name_to_id, bus_name, -1)
                
                push!(dc_comparison_df, [
                    bus_id,
                    bus_name,
                    vm,
                    ref_mag,
                    mag_diff,
                    mag_error_percent,
                    mag_within_tol
                ])
            else
                @warn "DC bus $bus_name not found in calculation results, skipping comparison"
            end
        else
            @warn "Unknown bus type '$bus_type' for bus $bus_name, skipping comparison"
        end
    end
    
    # Add statistics
    # AC statistics
    total_ac_buses = nrow(ac_comparison_df)
    if total_ac_buses > 0
        ac_mag_match_count = count(ac_comparison_df.Mag_Within_Tolerance)
        ac_ang_match_count = count(ac_comparison_df.Ang_Within_Tolerance)
        
        println("\nAC Buses Statistics:")
        println("Total AC buses: $total_ac_buses")
        println("AC buses with voltage magnitude within tolerance: $ac_mag_match_count ($(round(ac_mag_match_count/total_ac_buses*100, digits=2))%)")
        println("AC buses with voltage angle within tolerance: $ac_ang_match_count ($(round(ac_ang_match_count/total_ac_buses*100, digits=2))%)")
        
        # Output AC buses with mismatches
        if ac_mag_match_count < total_ac_buses
            println("\nAC buses with voltage magnitude mismatch:")
            for row in eachrow(filter(row -> !row.Mag_Within_Tolerance, ac_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Mag), Ref=$(row.Ref_Volt_Mag), Diff=$(row.Mag_Diff), Error=$(round(row.Mag_Error_Percent, digits=2))%")
            end
        end
        
        if ac_ang_match_count < total_ac_buses
            println("\nAC buses with voltage angle mismatch:")
            for row in eachrow(filter(row -> !row.Ang_Within_Tolerance, ac_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Ang), Ref=$(row.Ref_Volt_Ang), Diff=$(row.Ang_Diff)")
            end
        end
    else
        println("No AC buses found for comparison!")
    end
    
    # DC statistics
    total_dc_buses = nrow(dc_comparison_df)
    if total_dc_buses > 0
        dc_mag_match_count = count(dc_comparison_df.Mag_Within_Tolerance)
        
        println("\nDC Buses Statistics:")
        println("Total DC buses: $total_dc_buses")
        println("DC buses with voltage magnitude within tolerance: $dc_mag_match_count ($(round(dc_mag_match_count/total_dc_buses*100, digits=2))%)")
        
        # Output DC buses with mismatches
        if dc_mag_match_count < total_dc_buses
            println("\nDC buses with voltage magnitude mismatch:")
            for row in eachrow(filter(row -> !row.Mag_Within_Tolerance, dc_comparison_df))
                println("  $(row.Bus_Name) (ID=$(row.Bus_ID)): Calc=$(row.Calc_Volt_Mag), Ref=$(row.Ref_Volt_Mag), Diff=$(row.Mag_Diff), Error=$(round(row.Mag_Error_Percent, digits=2))%")
            end
        end
    end
    
    return (ac = ac_comparison_df, dc = dc_comparison_df)
end

"""
    save_comparison_results(comparison_results, output_file::String)

Save comparison results to an Excel file. Supports both original format and AC/DC mixed format.
Parameters:
- comparison_results: Can be a DataFrame or a NamedTuple containing ac and dc DataFrames
- output_file: Output file path
"""
function save_comparison_results(comparison_results, output_file::String)
    # Delete the file if it already exists
    if isfile(output_file)
        rm(output_file)
    end
    
    # Determine result type
    if isa(comparison_results, DataFrame)
        # Original format: single DataFrame
        XLSX.writetable(output_file, 
            Comparison = (collect(eachcol(comparison_results)), names(comparison_results))
        )
    else
        # AC/DC mixed format: NamedTuple containing ac and dc
        # Create parameters for worksheets to write
        sheets = Dict()
        
        # Add AC results (if any)
        if !isempty(comparison_results.ac)
            sheets[:AC_Comparison] = (collect(eachcol(comparison_results.ac)), names(comparison_results.ac))
        end
        
        # Add DC results (if any)
        if !isempty(comparison_results.dc)
            sheets[:DC_Comparison] = (collect(eachcol(comparison_results.dc)), names(comparison_results.dc))
        end
        
        # Write to Excel file
        XLSX.writetable(output_file; sheets...)
    end
    
    println("Comparison results saved to $output_file")
end



"""
    get_bus_voltage_results(results::NamedTuple, case::JuliaPowerCase)

Extract bus voltage information from power flow calculation results in JPC format, and map to node names in JuliaPowerCase.
Parameters:
- results: NamedTuple containing results in JPC format (typically the return value of power flow calculation)
- case: Original JuliaPowerCase object, used to get node names and ID mappings

Returns:
- DataFrame containing bus voltage results
"""
function get_bus_voltage_results(results::NamedTuple, case::JuliaPowerCase)
    # Extract JPC object from results
    jpc = results.value[1]
    
    # Create mapping from bus ID to name
    bus_id_to_name = Dict{Int, String}()
    for bus in case.busesAC
        if bus.in_service
            bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # Create results dataframe
    voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Volt_Ang = Float64[],
        Bus_Type = Int[]
    )
    
    # Iterate through all buses in JPC
    for i in 1:size(jpc.busAC, 1)
        bus_id = Int(jpc.busAC[i, 1])
        vm = jpc.busAC[i, 8]  # Voltage magnitude is typically in column 8
        va = jpc.busAC[i, 9]  # Voltage angle is typically in column 9
        bus_type = Int(jpc.busAC[i, 2])  # Bus type
        
        # Get bus name
        bus_name = get(bus_id_to_name, bus_id, "Unknown_Bus_$bus_id")
        
        push!(voltage_df, [
            bus_id,
            bus_name,
            vm,
            va,
            bus_type
        ])
    end
    
    # Sort by bus ID
    sort!(voltage_df, :Bus_ID)
    
    return voltage_df
end

"""
    get_bus_voltage_results_acdc(results::NamedTuple, case::JuliaPowerCase)

Extract both AC and DC bus voltage information from power flow calculation results in JPC format, and map to node names in JuliaPowerCase.
Parameters:
- results: NamedTuple containing results in JPC format (typically the return value of power flow calculation)
- case: Original JuliaPowerCase object, used to get node names and ID mappings

Returns:
- NamedTuple containing AC and DC bus voltage results as two DataFrames: ac and dc
"""
function get_bus_voltage_results_acdc(results::NamedTuple, case::JuliaPowerCase)
    # Create mapping from AC bus ID to name
    ac_bus_id_to_name = Dict{Int, String}()
    for bus in case.busesAC
        if bus.in_service
            ac_bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # Create mapping from DC bus ID to name
    dc_bus_id_to_name = Dict{Int, String}()
    for bus in case.busesDC
        if bus.in_service
            dc_bus_id_to_name[bus.bus_id] = bus.name
        end
    end
    
    # Create AC results dataframe
    ac_voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Volt_Ang = Float64[],
        Bus_Type = Int[],
        Island = Int[]  # Add island identifier
    )
    
    # Create DC results dataframe
    dc_voltage_df = DataFrame(
        Bus_ID = Int[],
        Bus_Name = String[],
        Volt_Mag = Float64[],
        Bus_Type = Int[],
        Island = Int[]  # Add island identifier
    )
    
    # Iterate through all islands
    n_islands = length(results.value)
    for island_idx in 1:n_islands
        jpc = results.value[island_idx]
        
        # Process AC buses
        if size(jpc.busAC, 1) > 0
            # Iterate through all AC buses in JPC
            for i in 1:size(jpc.busAC, 1)
                bus_id = Int(jpc.busAC[i, 1])
                vm = jpc.busAC[i, 8]  # Voltage magnitude is typically in column 8
                va = jpc.busAC[i, 9]  # Voltage angle is typically in column 9
                bus_type = Int(jpc.busAC[i, 2])  # Bus type
                
                # Get bus name
                bus_name = get(ac_bus_id_to_name, bus_id, "Unknown_AC_Bus_$bus_id")
                
                push!(ac_voltage_df, [
                    bus_id,
                    bus_name,
                    vm,
                    va,
                    bus_type,
                    island_idx  # Add island number
                ])
            end
        end
        
        # Process DC buses
        if size(jpc.busDC, 1) > 0
            # Iterate through all DC buses in JPC
            for i in 1:size(jpc.busDC, 1)
                bus_id = Int(jpc.busDC[i, 1])
                vm = jpc.busDC[i, 8]  # DC voltage value is typically in column 8
                bus_type = Int(jpc.busDC[i, 2])  # Bus type
                
                # Get bus name
                bus_name = get(dc_bus_id_to_name, bus_id, "Unknown_DC_Bus_$bus_id")
                
                push!(dc_voltage_df, [
                    bus_id,
                    bus_name,
                    vm,
                    bus_type,
                    island_idx  # Add island number
                ])
            end
        end
    end
    
    # Sort by Bus_ID only
    sort!(ac_voltage_df, :Bus_ID)
    sort!(dc_voltage_df, :Bus_ID)
    
    return (ac = ac_voltage_df, dc = dc_voltage_df)
end



"""
    plot_voltage_errors(comparison_results, output_file::String="voltage_errors.png";
                       show_plot::Bool=true)

Plot voltage magnitude and angle error curves. Supports both original format and AC/DC mixed format.
Parameters:
- comparison_results: Can be a DataFrame or a NamedTuple containing ac and dc DataFrames
- output_file: Path to save the chart file (default: "voltage_errors.png")
- show_plot: Whether to display the chart (default: true)

Returns:
- Generated chart object or NamedTuple containing chart objects
"""
function plot_voltage_errors(comparison_results, output_file::String="voltage_errors.png";
                            show_plot::Bool=true)
    # Set font, using system default font
    default_font = "Arial"  # Use a commonly available Western font
    
    # Determine result type
    if isa(comparison_results, DataFrame)
        # Original format: single DataFrame
        # Sort buses by ID
        sorted_df = sort(comparison_results, :Bus_ID)
        
        # Create a chart layout with 2 subplots
        plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                   fontfamily=default_font)
        
        # 1. Voltage magnitude error percentage curve
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Mag_Error_Percent, 
              label="Magnitude Error %", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Magnitude Error Percentage", 
              xlabel="Bus ID", ylabel="Error Percentage (%)")
        # Add zero error reference line
        hline!(plt[1], [0], linestyle=:dash, color=:black, label="Zero Error")
        
        # 2. Voltage angle error curve
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Ang_Diff, 
              label="Angle Error", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Angle Error", 
              xlabel="Bus ID", ylabel="Error (degrees)")
        # Add zero error reference line
        hline!(plt[2], [0], linestyle=:dash, color=:black, label="Zero Error")
        
        # Add overall title
        plot!(plt, title="Voltage Calculation Error Analysis", titlefontsize=14)
        
        # Save chart
        savefig(plt, output_file)
        
        if show_plot
            display(plt)
        end
        
        return plt
    else
        # AC/DC mixed format: NamedTuple containing ac and dc
        # Ensure output directory exists
        output_dir = dirname(output_file)
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        
        results = Dict()
        
        # Process AC part
        if !isempty(comparison_results.ac)
            # Sort buses by ID
            sorted_ac_df = sort(comparison_results.ac, :Bus_ID)
            
            # Create a chart layout with 2 subplots
            ac_plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 1. Voltage magnitude error percentage curve
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Mag_Error_Percent, 
                  label="Magnitude Error %", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Magnitude Error Percentage", 
                  xlabel="Bus ID", ylabel="Error Percentage (%)")
            # Add zero error reference line
            hline!(ac_plt[1], [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # 2. Voltage angle error curve
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Ang_Diff, 
                  label="Angle Error", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Angle Error", 
                  xlabel="Bus ID", ylabel="Error (degrees)")
            # Add zero error reference line
            hline!(ac_plt[2], [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # Add overall title
            plot!(ac_plt, title="AC Voltage Calculation Error Analysis", titlefontsize=14)
            
            # Save chart
            ac_error_plot_file = joinpath(dirname(output_file), "ac_voltage_errors.png")
            savefig(ac_plt, ac_error_plot_file)
            
            if show_plot
                display(ac_plt)
            end
            
            results[:ac] = (plot = ac_plt, file = ac_error_plot_file)
        end
        
        # Process DC part
        if !isempty(comparison_results.dc)
            # Sort buses by ID
            sorted_dc_df = sort(comparison_results.dc, :Bus_ID)
            
            # Create a chart
            dc_plt = plot(size=(800, 500), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # Voltage magnitude error percentage curve
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Mag_Error_Percent, 
                  label="Magnitude Error %", marker=:circle, markersize=4, 
                  linewidth=2, title="DC Voltage Magnitude Error Percentage", 
                  xlabel="Bus ID", ylabel="Error Percentage (%)")
            # Add zero error reference line
            hline!(dc_plt, [0], linestyle=:dash, color=:black, label="Zero Error")
            
            # Save chart
            dc_error_plot_file = joinpath(dirname(output_file), "dc_voltage_errors.png")
            savefig(dc_plt, dc_error_plot_file)
            
            if show_plot
                display(dc_plt)
            end
            
            results[:dc] = (plot = dc_plt, file = dc_error_plot_file)
        end
        
        return (ac = get(results, :ac, nothing), dc = get(results, :dc, nothing))
    end
end


"""
    plot_voltage_comparison(comparison_results, output_file::String="voltage_comparison.png";
                           show_plot::Bool=true)

Plot voltage magnitude and angle comparison curves between calculated and reference values. Supports both original format and AC/DC mixed format.
Parameters:
- comparison_results: Can be a DataFrame or a NamedTuple containing ac and dc DataFrames
- output_file: Path to save the chart file (default: "voltage_comparison.png")
- show_plot: Whether to display the chart (default: true)

Returns:
- Generated chart object or NamedTuple containing chart objects
"""
function plot_voltage_comparison(comparison_results, output_file::String="voltage_comparison.png";
                                show_plot::Bool=true)
    # Set font, using system default font
    default_font = "Arial"  # Use a commonly available Western font
    
    # Determine result type
    if isa(comparison_results, DataFrame)
        # Original format: single DataFrame
        # Sort buses by ID
        sorted_df = sort(comparison_results, :Bus_ID)
        
        # Create a chart layout with 2 subplots
        plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                   fontfamily=default_font)
        
        # 1. Voltage magnitude comparison (calculated vs reference)
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Calc_Volt_Mag, 
              label="Calculated", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Magnitude Comparison", 
              xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
        plot!(plt[1], sorted_df.Bus_ID, sorted_df.Ref_Volt_Mag, 
              label="Reference", marker=:square, markersize=4, 
              linewidth=2, linestyle=:dash)
        
        # 2. Voltage angle comparison (calculated vs reference)
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Calc_Volt_Ang, 
              label="Calculated", marker=:circle, markersize=4, 
              linewidth=2, title="Voltage Angle Comparison", 
              xlabel="Bus ID", ylabel="Voltage Angle (degrees)")
        plot!(plt[2], sorted_df.Bus_ID, sorted_df.Ref_Volt_Ang, 
              label="Reference", marker=:square, markersize=4, 
              linewidth=2, linestyle=:dash)
        
        # Add overall title
        plot!(plt, title="Voltage Calculation vs Reference Comparison", titlefontsize=14)
        
        # Save chart
        savefig(plt, output_file)
        
        if show_plot
            display(plt)
        end
        
        return plt
    else
        # AC/DC mixed format: NamedTuple containing ac and dc
        # Ensure output directory exists
        output_dir = dirname(output_file)
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        
        results = Dict()
        
        # Process AC part
        if !isempty(comparison_results.ac)
            # Sort buses by ID
            sorted_ac_df = sort(comparison_results.ac, :Bus_ID)
            
            # Create a chart layout with 2 subplots
            ac_plt = plot(layout=(2,1), size=(1000, 800), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
            # 1. Voltage magnitude comparison (calculated vs reference)
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Calc_Volt_Mag, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Magnitude Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
            plot!(ac_plt[1], sorted_ac_df.Bus_ID, sorted_ac_df.Ref_Volt_Mag, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # 2. Voltage angle comparison (calculated vs reference)
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Calc_Volt_Ang, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="AC Voltage Angle Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Angle (degrees)")
            plot!(ac_plt[2], sorted_ac_df.Bus_ID, sorted_ac_df.Ref_Volt_Ang, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # Add overall title
            plot!(ac_plt, title="AC Voltage Calculation vs Reference Comparison", titlefontsize=14)
            
            # Save chart
            ac_comparison_plot_file = joinpath(dirname(output_file), "ac_voltage_comparison.png")
            savefig(ac_plt, ac_comparison_plot_file)
            
            if show_plot
                display(ac_plt)
            end
            
            results[:ac] = (plot = ac_plt, file = ac_comparison_plot_file)
        end
        
        # Process DC part
        if !isempty(comparison_results.dc)
            # Sort buses by ID
            sorted_dc_df = sort(comparison_results.dc, :Bus_ID)
            
            # Create a chart
            dc_plt = plot(size=(800, 500), dpi=300, legend=:outertopright,
                         fontfamily=default_font)
            
                        # Voltage magnitude comparison (calculated vs reference)
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Calc_Volt_Mag, 
                  label="Calculated", marker=:circle, markersize=4, 
                  linewidth=2, title="DC Voltage Magnitude Comparison", 
                  xlabel="Bus ID", ylabel="Voltage Magnitude (p.u.)")
            plot!(dc_plt, sorted_dc_df.Bus_ID, sorted_dc_df.Ref_Volt_Mag, 
                  label="Reference", marker=:square, markersize=4, 
                  linewidth=2, linestyle=:dash)
            
            # Save chart
            dc_comparison_plot_file = joinpath(dirname(output_file), "dc_voltage_comparison.png")
            savefig(dc_plt, dc_comparison_plot_file)
            
            if show_plot
                display(dc_plt)
            end
            
            results[:dc] = (plot = dc_plt, file = dc_comparison_plot_file)
        end
        
        return (ac = get(results, :ac, nothing), dc = get(results, :dc, nothing))
    end
end



"""
    analyze_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String;
                           tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3,
                           output_dir::String="./results", save_pdf::Bool=true)

Analyze voltage differences between power flow calculation results and reference file, generate comparison reports and charts.
Supports both original format and AC/DC mixed format.
Parameters:
- results: NamedTuple containing results in JPC format (typically the return value of power flow calculation)
- case: Original JuliaPowerCase object, used to get node names and ID mappings
- reference_file: Path to Excel file containing reference voltage values
- tolerance_mag: Tolerance for voltage magnitude comparison (default: 1e-4)
- tolerance_ang: Tolerance for voltage angle comparison (default: 1e-3)
- output_dir: Output directory (default: "./results")
- save_pdf: Whether to save plots in PDF format in addition to PNG (default: true)

Returns:
- DataFrame or NamedTuple containing comparison results
"""
function analyze_voltage_results(results::NamedTuple, case::JuliaPowerCase, reference_file::String;
                                tolerance_mag::Float64=1e-4, tolerance_ang::Float64=1e-3,
                                output_dir::String="./results", save_pdf::Bool=true)
    # Ensure output directory exists
    mkpath(output_dir)
    
    # Compare voltage results
    comparison_results = compare_voltage_results(results, case, reference_file, 
                                               tolerance_mag=tolerance_mag, 
                                               tolerance_ang=tolerance_ang)
    
    # Save comparison results to Excel file
    excel_output = joinpath(output_dir, "voltage_comparison_results.xlsx")
    save_comparison_results(comparison_results, excel_output)
    
    # Plot voltage error curves
    error_plot_file = joinpath(output_dir, "voltage_errors.png")
    error_plots = plot_voltage_errors(comparison_results, error_plot_file)
    
    # Plot voltage calculated vs reference comparison curves
    comparison_plot_file = joinpath(output_dir, "voltage_comparison.png")
    comparison_plots = plot_voltage_comparison(comparison_results, comparison_plot_file)
    
    # Save plots in PDF format if requested
    if save_pdf
        if isa(comparison_results, DataFrame)
            # Original format
            pdf_error_plot_file = joinpath(output_dir, "voltage_errors.pdf")
            savefig(error_plots, pdf_error_plot_file)
            
            pdf_comparison_plot_file = joinpath(output_dir, "voltage_comparison.pdf")
            savefig(comparison_plots, pdf_comparison_plot_file)
        else
            # AC/DC mixed format
            if !isempty(comparison_results.ac)
                pdf_ac_error_plot_file = joinpath(output_dir, "ac_voltage_errors.pdf")
                savefig(error_plots.ac.plot, pdf_ac_error_plot_file)
                
                pdf_ac_comparison_plot_file = joinpath(output_dir, "ac_voltage_comparison.pdf")
                savefig(comparison_plots.ac.plot, pdf_ac_comparison_plot_file)
            end
            
            if !isempty(comparison_results.dc)
                pdf_dc_error_plot_file = joinpath(output_dir, "dc_voltage_errors.pdf")
                savefig(error_plots.dc.plot, pdf_dc_error_plot_file)
                
                pdf_dc_comparison_plot_file = joinpath(output_dir, "dc_voltage_comparison.pdf")
                savefig(comparison_plots.dc.plot, pdf_dc_comparison_plot_file)
            end
        end
    end
    
    println("\nAnalysis complete!")
    println("Comparison results saved to: $excel_output")
    
    # Output different information based on result type
    if isa(comparison_results, DataFrame)
        # Original format
        println("Voltage error plots saved to: $error_plot_file")
        println("Voltage comparison plots saved to: $comparison_plot_file")
        if save_pdf
            println("PDF voltage error plots saved to: $(joinpath(output_dir, "voltage_errors.pdf"))")
            println("PDF voltage comparison plots saved to: $(joinpath(output_dir, "voltage_comparison.pdf"))")
        end
    else
        # AC/DC mixed format
        if !isempty(comparison_results.ac)
            println("AC voltage error plots saved to: $(error_plots.ac.file)")
            println("AC voltage comparison plots saved to: $(comparison_plots.ac.file)")
            if save_pdf
                println("PDF AC voltage error plots saved to: $(joinpath(output_dir, "ac_voltage_errors.pdf"))")
                println("PDF AC voltage comparison plots saved to: $(joinpath(output_dir, "ac_voltage_comparison.pdf"))")
            end
        end
        
        if !isempty(comparison_results.dc)
            println("DC voltage error plots saved to: $(error_plots.dc.file)")
            println("DC voltage comparison plots saved to: $(comparison_plots.dc.file)")
            if save_pdf
                println("PDF DC voltage error plots saved to: $(joinpath(output_dir, "dc_voltage_errors.pdf"))")
                println("PDF DC voltage comparison plots saved to: $(joinpath(output_dir, "dc_voltage_comparison.pdf"))")
            end
        end
    end
    
    return comparison_results
end

