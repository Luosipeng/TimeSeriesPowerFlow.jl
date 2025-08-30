"""
    load_julia_power_data(file_path::String)

Load power system data from an Excel file and convert it to a JuliaPowerCase structure.

Parameters:
- `file_path::String`: Path to the Excel file

Returns:
- `JuliaPowerCase`: Structure containing all power system components
"""
function load_julia_power_data(file_path::String)
    # Verify if the file exists
    if !isfile(file_path)
        error("File does not exist: $file_path")
    end
    
    # Create empty JuliaPowerCase structure
    case = JuliaPowerCase()
    
    try
        # Read Excel file
        @info "Reading Excel file: $file_path"
        xf = XLSX.readxlsx(file_path)
        
        # Get all worksheet names
        sheet_names = XLSX.sheetnames(xf)
        @info "Found worksheets: $(join(sheet_names, ", "))"
        
        # Define component loaders and corresponding worksheet names
        component_loaders = [
            ("bus", load_buses!),
            ("dcbus", load_buses!),
            ("cable", load_lines!),
            ("xline", load_lines!),
            ("dcimpedance", load_dclines!),
            ("sgen", load_static_generators!),
            ("lumpedload", load_loads!),
            ("dclumpload", load_dcloads!),
            ("xform2w", load_trafo!),
            ("xform3w", load_trafo3ws!),
            ("gen", load_generators!),
            ("battery", load_storages!),
            ("inverter", load_converters!),
            ("vpp", load_virtual_power_plants!),
            ("util", load_ext_grids!),
            ("hvcb", load_switches! ),
            ("pvarray", load_pv_arrays!),
            ("pvarray", load_ac_pv_system!),
            # ("switch", load_switches!), # To do
            # ("equipment_carbon", load_carbons!), # To do 
            # ("carbon_time_series", load_carbon_time_series!), # To do
            # ("carbon_scenario", load_carbon_scenarios!), # To do 
            # ("charging_station", load_charging_stations!), # To do
            # ("charger", load_chargers!), # To do
            # ("v2g_service", load_v2g_services!), # To do
            # ("mobile_storage", load_mess!), # To do 
            # ("load_time_series", load_time_series!), # To do
            # ("res_time_series", res_time_series!), # To do
        ]
        
        # Load component data
        for (sheet_name, loader_func) in component_loaders
            if sheet_name in sheet_names
                try
                    @info "Loading $sheet_name data..."
                    loader_func(case, file_path, sheet_name)
                catch e
                    @error "Error loading $sheet_name data" exception=(e, catch_backtrace())
                end
            else
                @debug "$sheet_name worksheet does not exist, skipping"
            end
        end
        
        # Validate loaded data
        validate_case(case)
        
    catch e
        @error "Error loading power system data" exception=(e, catch_backtrace())
        rethrow(e)
    end
    
    # Return filled JuliaPowerCase structure
    return case
end

"""
    validate_case(case::JuliaPowerCase)

Validate the integrity and consistency of the loaded power system case data.
"""
function validate_case(case::JuliaPowerCase)
    # Check if bus data exists
    if isempty(case.busesAC)
        @warn "Warning: No bus data loaded"
    end
    
    # Check line and bus associations
    if !isempty(case.branchesAC)
        bus_indices = Set(bus.index for bus in case.busesAC)
        for line in case.branchesAC
            if !(line.from_bus in bus_indices) || !(line.to_bus in bus_indices)
                @warn "Warning: Line $(line.name) (ID: $(line.index)) connects to non-existent buses"
            end
        end
    end
    
    # More validation logic can be added...
    
    @info "Power system case data validation complete"
end

"""
    safe_get_value(cell, default_value, type_converter=identity)

Safely get value from a cell, providing type conversion and default value.

Parameters:
- `cell`: Excel cell value
- `default_value`: Default value if cell is empty or conversion fails
- `type_converter`: Type conversion function or target type

Returns:
- Converted value or default value
"""
function safe_get_value(cell, default_value, type_converter=identity)
    if ismissing(cell) || cell === nothing || (typeof(cell) <: AbstractString && isempty(strip(string(cell))))
        return default_value
    else
        try
            # Check if type_converter is a type or function
            if type_converter isa DataType
                # If it's a type, create an appropriate conversion function
                if type_converter <: Number && typeof(cell) <: AbstractString
                    return parse(type_converter, cell)
                else
                    return convert(type_converter, cell)
                end
            else
                # If it's a function, use it directly
                return type_converter(cell)
            end
        catch e
            @debug "Value conversion failed: $cell to $(typeof(default_value)) type" exception=e
            return default_value
        end
    end
end



"""
    parse_bool(value)

Parse boolean values from various types.

Parameters:
- `value`: Value to parse

Returns:
- Parsed boolean value
"""
function parse_bool(value)
    if typeof(value) <: Bool
        return value
    elseif typeof(value) <: AbstractString
        lowercase_value = lowercase(strip(value))
        if lowercase_value in ["true", "yes", "1", "t", "y"]
            return true
        elseif lowercase_value in ["false", "no", "0", "f", "n"]
            return false
        else
            @debug "Cannot parse boolean value: $value defaulting to false"
            return false
        end
    elseif typeof(value) <: Number
        return value != 0
    else
        @debug "Cannot parse boolean type: $(typeof(value)) defaulting to false"
        return false
    end
end

"""
    load_buses!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load bus data from Excel file and add to power system case.
Also creates mappings from bus names to integer IDs and from zone names to integer IDs.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing bus data
"""
function load_buses!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        # Ensure data is not empty
        if isempty(df)
            @info "Bus table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        if sheet_name == "bus"
            required_columns = [:index, :id, :nominalkv]
        elseif sheet_name == "dcbus"
            required_columns = [:index, :id, :nominalv]
        end
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Bus table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # (1) Map df[:name] to consecutive integers starting from 1
        name_to_id = Dict{String, Int}()
        for (i, name) in enumerate(df[:, :id])
            if !ismissing(name) && !isempty(strip(string(name)))
                name_str = string(name)
                if !haskey(name_to_id, name_str)
                    name_to_id[name_str] = i
                end
            end
        end
        
        
        # (2) Map df[:zone] to consecutive integer vector starting from 1
        zone_to_id = Dict{String, Int}()
        
        # Check if zone column exists
        has_zone_column = any(col -> lowercase(string(col)) == "zone", names(df))
        
        if has_zone_column
            # First assign consecutive IDs for each unique zone
            unique_zones = []
            for zone in df[:, :zone]
                zone_str = safe_get_value(zone, "", String)
                if !isempty(zone_str) && !(zone_str in unique_zones)
                    push!(unique_zones, zone_str)
                end
            end
            
            # Assign consecutive IDs starting from 1 for unique zones
            for (i, zone) in enumerate(unique_zones)
                zone_to_id[zone] = i
            end
        else
            @info "No zone column in bus table, skipping zone mapping creation"
        end
        
        #(3) Map df[:area] to integer IDs
        area_to_id = Dict{String, Int}()
        has_area_column = any(col -> lowercase(string(col)) == "area", names(df))
        if has_area_column
            # First assign consecutive IDs for each unique area
            unique_areas = []
            for area in df[:, :area]
                area_str = safe_get_value(area, "", String)
                if !isempty(area_str) && !(area_str in unique_areas)
                    push!(unique_areas, area_str)
                end
            end
            
            # Assign consecutive IDs starting from 1 for unique areas
            for (i, area) in enumerate(unique_areas)
                area_to_id[area] = i
            end
        else
            @info "No area column in bus table, skipping area mapping creation"
        end

        # Save mappings to case
        if sheet_name == "bus"
            case.bus_name_to_id = name_to_id
        elseif sheet_name == "dcbus"
            # For DC buses, use different mapping
            case.busdc_name_to_id = name_to_id
        end
        case.zone_to_id = zone_to_id
        case.area_to_id = area_to_id
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Based on the index information to check whether it exists or not
                if index <= 0
                    @warn "Row $i: Invalid bus index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                bus_id = name_to_id[name]
                zone = haskey(row, :zone) ? safe_get_value(row[:zone], "", String) : ""
                area = haskey(row, :area) ? safe_get_value(row[:area], "", String) : ""
                if sheet_name == "bus"
                    vn_kv = safe_get_value(row[:nominalkv], 0.0, Float64)
                elseif sheet_name == "dcbus"
                    vn_kv = safe_get_value(row[:nominalv], 0.0, Float64)/1000.0
                end

                if has_zone_column
                    zone_id = zone_to_id[zone]
                else
                    zone_id = 1
                end
                if has_area_column
                    area_id = area_to_id[area]
                else
                    area_id = 1
                end
                # Validate voltage level is reasonable
                if vn_kv <= 0.0
                    @warn "Row $i: Bus $name (ID: $index) has invalid voltage level ($vn_kv kV), using default value 0.4 kV"
                    vn_kv = 0.4
                end
                
                max_vm_pu = haskey(row, :max_vm_pu) ? safe_get_value(row[:max_vm_pu], 1.05, Float64) : 1.05
                min_vm_pu = haskey(row, :min_vm_pu) ? safe_get_value(row[:min_vm_pu], 0.95, Float64) : 0.95
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Validate voltage limits are reasonable
                if min_vm_pu >= max_vm_pu
                    @warn "Row $i: Bus $name (ID: $index) has invalid voltage limits (min: $min_vm_pu, max: $max_vm_pu), using default values"
                    min_vm_pu = 0.95
                    max_vm_pu = 1.05
                end
                
                # Create Bus object and add to case
                if sheet_name == "bus"
                    push!(case.busesAC, Bus(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id))
                elseif sheet_name == "dcbus"
                    push!(case.busesDC, BusDC(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id))
                    
                end
                processed_rows += 1

            catch e
                @error "Error processing bus data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Bus data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        if sheet_name == "bus"
            @info "Created $(length(case.bus_name_to_id)) bus name to ID mappings"
        else
            @info "Created $(length(case.busdc_name_to_id)) bus name to ID mappings"
        end
        if has_zone_column
            @info "Created $(length(case.zone_to_id)) zone name to ID mappings"
        end
        
    catch e
        @error "Error loading bus data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end



"""
    load_lines!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load line data from Excel file and add to power system case.
Uses case.bus_name_to_id to map bus names to integer IDs.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing line data
"""
function load_lines!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading line data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Line table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        if sheet_name == "cable"
            required_columns = [:index, :frombus, :tobus, :lengthvalue, :cablelengthunit, :ohmsperlengthunit, :ohmsperlengthvalue,:rposvalue, :xposvalue,:rzerovalue, :xzerovalue]
        elseif sheet_name == "xline"
            required_columns = [:index, :frombus, :tobus, :length, :lengthunit, :perlength, :perlengthunit, :rpos, :xpos,:rzero, :xzero]
        end
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Line table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Check if bus_name_to_id mapping exists
        if !isdefined(case, :bus_name_to_id) || isempty(case.bus_name_to_id)
            @warn "case.bus_name_to_id mapping does not exist or is empty, cannot map bus names to IDs"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid line index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Get bus names from data
                from_bus_name = safe_get_value(row[:frombus], "", String)
                to_bus_name = safe_get_value(row[:tobus], "", String)
                
                # Use bus_name_to_id mapping to convert bus names to integer IDs
                from_bus = 0
                to_bus = 0
                
                if haskey(case.bus_name_to_id, from_bus_name)
                    from_bus = case.bus_name_to_id[from_bus_name]
                else
                    @warn "Row $i: Line $name (ID: $index) has from-bus name '$from_bus_name' that doesn't exist in mapping, skipping this row"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.bus_name_to_id, to_bus_name)
                    to_bus = case.bus_name_to_id[to_bus_name]
                else
                    @warn "Row $i: Line $name (ID: $index) has to-bus name '$to_bus_name' that doesn't exist in mapping, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus indices are valid
                if from_bus <= 0 || to_bus <= 0
                    @warn "Row $i: Line $name (ID: $index) connects to invalid bus IDs (from: $from_bus, to: $to_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate buses exist
                bus_indices = Set(bus.bus_id for bus in case.busesAC)
                if !(from_bus in bus_indices) || !(to_bus in bus_indices)
                    @warn "Row $i: Line $name (ID: $index) connects to non-existent bus IDs (from: $from_bus, to: $to_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate line is not a self-loop
                if from_bus == to_bus
                    @warn "Row $i: Line $name (ID: $index) connects to the same bus ($from_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                
                if sheet_name == "cable"
                    length_value = haskey(row, :lengthvalue) ? safe_get_value(row[:lengthvalue], 0.0, Float64) : 0.0
                    length_unit = safe_get_value(row[:cablelengthunit], "", String)
                else sheet_name == "xline"
                    length_value = haskey(row, :length) ? safe_get_value(row[:length], 0.0, Float64) : 0.0
                    length_unit = safe_get_value(row[:lengthunit], "", String)
                end

                if length_unit =="0"
                    length_km = length_value * 0.0003048 # Convert feet to kilometers
                elseif length_unit == "1"
                    length_km = length_value * 1.60934 # Convert miles to kilometers
                elseif length_unit == "2"
                    length_km = length_value * 0.001 # Convert meters to kilometers
                elseif length_unit == "3"
                    length_km = length_value * 1.0 # Convert kilometers to kilometers
                end
                
                # Validate length is reasonable
                if length_km < 0.0
                    @warn "Row $i: Line $name (ID: $index) has invalid length ($length_km km), setting to 0"
                    length_km = 0.0
                end
                if sheet_name == "cable"
                    OhmsPerLengthUnit = safe_get_value(row[:ohmsperlengthunit], "", String)
                    OhmsPerLengthValue = safe_get_value(row[:ohmsperlengthvalue], 0.0, Float64)
                    RPosValue = safe_get_value(row[:rposvalue], 0.0, Float64)
                    XPosValue = safe_get_value(row[:xposvalue], 0.0, Float64)
                    CPosValue = safe_get_value(row[:yposvalue], 0.0, Float64)
                    RZeroValue = safe_get_value(row[:rzerovalue], 0.0, Float64)
                    XZeroValue = safe_get_value(row[:xzerovalue], 0.0, Float64)
                    CZeroValue = safe_get_value(row[:yzerovalue], 0.0, Float64)
                elseif sheet_name == "xline"
                    OhmsPerLengthUnit = safe_get_value(row[:perlengthunit], "", String)
                    OhmsPerLengthValue = safe_get_value(row[:perlength], 0.0, Float64)
                    RPosValue = safe_get_value(row[:rpos], 0.0, Float64)
                    XPosValue = safe_get_value(row[:xpos], 0.0, Float64)
                    CPosValue = safe_get_value(row[:ypos], 0.0, Float64)
                    RZeroValue = safe_get_value(row[:rzero], 0.0, Float64)
                    XZeroValue = safe_get_value(row[:xzero], 0.0, Float64)
                    CZeroValue = safe_get_value(row[:yzero], 0.0, Float64)
                end

                if OhmsPerLengthUnit == "0"
                    ohmsperkm = OhmsPerLengthValue * 0.0003048 # Convert feet to kilometers
                elseif OhmsPerLengthUnit == "1"
                    ohmsperkm = OhmsPerLengthValue * 1.60934 # Convert miles to kilometers
                elseif OhmsPerLengthUnit == "2"
                    ohmsperkm = OhmsPerLengthValue * 0.001 # Convert meters to kilometers
                elseif OhmsPerLengthUnit == "3"
                    ohmsperkm = OhmsPerLengthValue * 1.0 # Convert kilometers to kilometers
                end

                r_ohm_per_km = RPosValue/ohmsperkm
                x_ohm_per_km = XPosValue/ohmsperkm
                c_nf_per_km = CPosValue/ohmsperkm
                r0_ohm_per_km = RZeroValue/ohmsperkm
                x0_ohm_per_km = XZeroValue/ohmsperkm
                c0_nf_per_km = CZeroValue/ohmsperkm
                g_us_per_km = CZeroValue/ohmsperkm
                
                
                # Validate impedance is reasonable
                if r_ohm_per_km < 0.0
                    @warn "Row $i: Line $name (ID: $index) has invalid resistance value ($r_ohm_per_km Ω/km), using default value 0.1"
                    r_ohm_per_km = 0.1
                end
                
                if x_ohm_per_km < 0.0
                    @warn "Row $i: Line $name (ID: $index) has invalid reactance value ($x_ohm_per_km Ω/km), using default value 0.1"
                    x_ohm_per_km = 0.1
                end
                
                max_i_ka = haskey(row, :max_i_ka) ? safe_get_value(row[:max_i_ka], 0.0, Float64) : 0.0
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                max_loading_percent = haskey(row, :max_loading_percent) ? safe_get_value(row[:max_loading_percent], 100.0, Float64) : 100.0
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # Validate parallel count is reasonable
                if parallel <= 0
                    @warn "Row $i: Line $name (ID: $index) has invalid parallel count ($parallel), setting to 1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Collect reliability parameters as named parameters
                reliability_params = Dict{Symbol, Any}()
                
                # Define reliability parameter list and default values
                reliability_fields = [
                    (:mtbf_hours, 0.0),
                    (:mttr_hours, 0.0),
                    (:failure_rate_per_year, 0.0),
                    (:planned_outage_hours_per_year, 0.0),
                    (:forced_outage_rate, 0.0),
                    (:permanent_fault_rate_per_km_year, 0.0),
                    (:temporary_fault_rate_per_km_year, 0.0),
                    (:repair_time_permanent_hours, 0.0),
                    (:auto_reclosing_success_rate, 0.0)
                ]
                
                # Safely extract reliability parameters
                for (field, default_value) in reliability_fields
                    field_str = String(field)
                    if haskey(row, Symbol(field_str))
                        reliability_params[field] = safe_get_value(row[Symbol(field_str)], default_value, Float64)
                    else
                        reliability_params[field] = default_value
                    end
                end
                
                # Create Line object and add to case, using named parameters for reliability parameters
                push!(case.branchesAC, Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                                      c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km,
                                      c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent,
                                      parallel, df, in_service; reliability_params...))
                
                processed_rows += 1
                
            catch e
                @error "Error processing line data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Line data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading line data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_dclines!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load DC line data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing DC line data
"""
function load_dclines!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading DC line data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "DC line table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :id, :frombus, :tobus, :rvalue, :lvalue]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "DC line table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid DC line index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus names to bus IDs
                from_bus_name = safe_get_value(row[:frombus], "", String)
                to_bus_name = safe_get_value(row[:tobus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus names to integer IDs
                from_bus = 0
                to_bus = 0
                
                                if haskey(case.busdc_name_to_id, from_bus_name)
                    from_bus = case.busdc_name_to_id[from_bus_name]
                else
                    @warn "Row $i: DC line $name (ID: $index) has from-bus name '$from_bus_name' that doesn't exist in busdc_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.busdc_name_to_id, to_bus_name)
                    to_bus = case.busdc_name_to_id[to_bus_name]
                else
                    @warn "Row $i: DC line $name (ID: $index) has to-bus name '$to_bus_name' that doesn't exist in busdc_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus indices are valid
                if from_bus <= 0 || to_bus <= 0
                    @warn "Row $i: DC line $name (ID: $index) connects to invalid buses (from: $from_bus, to: $to_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate buses exist
                bus_indices = Set(bus.index for bus in case.busesDC)
                if !(from_bus in bus_indices) || !(to_bus in bus_indices)
                    @warn "Row $i: DC line $name (ID: $index) connects to non-existent buses (from: $from_bus, to: $to_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate line is not a self-loop
                if from_bus == to_bus
                    @warn "Row $i: DC line $name (ID: $index) connects to the same bus ($from_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                if sheet_name == "dcimpedance"
                    r = safe_get_value(row[:rvalue], 0.0, Float64)
                    x = safe_get_value(row[:lvalue], 0.0, Float64)
                end
                
                # if x == 0.0
                #     @warn "Row $i: DC line $name (ID: $index) has invalid reactance value ($x), using default value 1e-6"
                #     x = 1e-6
                # end
                
                # Extract other field values from row data
                if sheet_name == "dcimpedance"
                    length_km = 1.0
                else
                    length_km = haskey(row, :length_km) ? safe_get_value(row[:length_km], 0.0, Float64) : 0.0
                end
                
                # Validate length is reasonable
                if length_km < 0.0
                    @warn "Row $i: DC line $name (ID: $index) has invalid length ($length_km km), setting to 0"
                    length_km = 0.0
                end
                if sheet_name == "dcimpedance"
                    r_ohm_per_km = r
                    x_ohm_per_km = x
                else
                    r_ohm_per_km = safe_get_value(row[:rvalue], 0.0, Float64)
                end
                
                # Validate impedance is reasonable
                if r_ohm_per_km < 0.0
                    @warn "Row $i: DC line $name (ID: $index) has invalid resistance value ($r_ohm_per_km Ω/km), using default value 0.1"
                    r_ohm_per_km = 0.1
                end
                
                g_us_per_km = haskey(row, :g_us_per_km) ? safe_get_value(row[:g_us_per_km], 0.0, Float64) : 0.0
                max_i_ka = haskey(row, :max_i_ka) ? safe_get_value(row[:max_i_ka], 0.0, Float64) : 0.0
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                max_loading_percent = haskey(row, :max_loading_percent) ? safe_get_value(row[:max_loading_percent], 100.0, Float64) : 100.0
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # Validate parallel count is reasonable
                if parallel <= 0
                    @warn "Row $i: DC line $name (ID: $index) has invalid parallel count ($parallel), setting to 1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : false
                
                # Create LineDC object and add to case
                push!(case.branchesDC, LineDC(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                                          g_us_per_km, max_i_ka, type, max_loading_percent,
                                          parallel, df, in_service))
                
                processed_rows += 1
                
            catch e
                @error "Error processing DC line data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "DC line data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading DC line data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_static_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load static generator data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing static generator data
"""
function load_static_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading static generator data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Static generator table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :bus, :p_mw, :q_mvar]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Static generator table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid static generator index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "Row $i: Static generator $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: Static generator $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus exists
                bus_indices = Set(b.index for b in case.buses)
                if !(bus in bus_indices)
                    @warn "Row $i: Static generator $name (ID: $index) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                p_mw = safe_get_value(row[:p_mw], 0.0, Float64)
                q_mvar = safe_get_value(row[:q_mvar], 0.0, Float64)
                scaling = haskey(row, :scaling) ? safe_get_value(row[:scaling], 1.0, Float64) : 1.0
                
                # Validate scaling factor is reasonable
                if scaling < 0.0
                    @warn "Row $i: Static generator $name (ID: $index) has invalid scaling factor ($scaling), setting to 1.0"
                    scaling = 1.0
                end
                
                max_p_mw = haskey(row, :max_p_mw) ? safe_get_value(row[:max_p_mw], 0.0, Float64) : 0.0
                min_p_mw = haskey(row, :min_p_mw) ? safe_get_value(row[:min_p_mw], 0.0, Float64) : 0.0
                max_q_mvar = haskey(row, :max_q_mvar) ? safe_get_value(row[:max_q_mvar], 0.0, Float64) : 0.0
                min_q_mvar = haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0
                
                # Validate power limits are reasonable
                if max_p_mw < min_p_mw
                    @warn "Row $i: Static generator $name (ID: $index) has invalid active power limits (min: $min_p_mw, max: $max_p_mw), swapping values"
                    max_p_mw, min_p_mw = min_p_mw, max_p_mw
                end
                
                if max_q_mvar < min_q_mvar
                    @warn "Row $i: Static generator $name (ID: $index) has invalid reactive power limits (min: $min_q_mvar, max: $max_q_mvar), swapping values"
                    max_q_mvar, min_q_mvar = min_q_mvar, max_q_mvar
                end
                
                k = haskey(row, :k) ? safe_get_value(row[:k], 0.0, Float64) : 0.0
                rx = haskey(row, :rx) ? safe_get_value(row[:rx], 0.0, Float64) : 0.0
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], false)) : false
                
                # Create StaticGenerator object and add to case
                push!(case.static_generators, StaticGenerator(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                                                           max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable))
                
                processed_rows += 1
                
            catch e
                @error "Error processing static generator data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Static generator data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading static generator data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_loads!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load load data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing load data
"""
function load_loads!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading load data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Load table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :bus, :mva, :pf]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Load table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid load index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "Row $i: Load $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: Load $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus exists
                bus_indices = Set(b.index for b in case.busesAC)
                if !(bus in bus_indices)
                    @warn "Row $i: Load $name (ID: $index) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                s_mva = safe_get_value(row[:mva], 0.0, Float64)/1000
                p_mw = s_mva * safe_get_value(row[:pf], 0.0, Float64)/100
                q_mvar = s_mva * sqrt(1 - (safe_get_value(row[:pf], 0.0, Float64)/100)^2)
                
                # Validate load model parameters
                const_p_percent = haskey(row, :mtloadpercent) ? safe_get_value(row[:mtloadpercent], 100.0, Float64) : 100.0
                const_z_percent = 100.0 - const_p_percent
                const_i_percent = haskey(row, :const_i_percent) ? safe_get_value(row[:const_i_percent], 0.0, Float64) : 0.0
                
                
                # Validate percentage sum is 100%
                total_percent = const_z_percent + const_i_percent + const_p_percent
                if abs(total_percent - 100.0) > 1e-6
                    @warn "Row $i: Load $name (ID: $index) has load model percentages that don't sum to 100% ($total_percent%), normalizing"
                    if total_percent > 0
                        const_z_percent = const_z_percent * 100.0 / total_percent
                        const_i_percent = const_i_percent * 100.0 / total_percent
                        const_p_percent = const_p_percent * 100.0 / total_percent
                    else
                        const_z_percent = 0.0
                        const_i_percent = 0.0
                        const_p_percent = 100.0
                    end
                end
                
                scaling = haskey(row, :uniformscalepq) ? safe_get_value(row[:uniformscalepq], 1.0, Float64) : 1.0
                
                # Validate scaling factor is reasonable
                if scaling < 0.0
                    @warn "Row $i: Load $name (ID: $index) has invalid scaling factor ($scaling), setting to 1.0"
                    scaling = 1.0
                end
                
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                type = haskey(row, :LoadType) ? safe_get_value(row[:type], "", String) : "wye"
                
                if type == "0"
                    type = "wye"
                elseif type == "1"
                    type = "delta"
                end
                
                # Create Load object and add to case
                push!(case.loadsAC, Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent,
                                     const_p_percent, scaling, in_service, type))
                
                processed_rows += 1
                
            catch e
                @error "Error processing load data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Load data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading load data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_dcloads!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load DC load data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing DC load data
"""
function load_dcloads!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading load data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Load table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :bus, :kw]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Load table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid load index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus = case.busdc_name_to_id[bus_name]
                else
                    @warn "Row $i: Load $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: Load $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus exists
                bus_indices = Set(b.index for b in case.busesDC)
                if !(bus in bus_indices)
                    @warn "Row $i: Load $name (ID: $index) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                p_mw = safe_get_value(row[:kw], 0.0, Float64)/1000
                
                # Validate load model parameters
                const_p_percent = haskey(row, :mtloadpercent) ? safe_get_value(row[:mtloadpercent], 100.0, Float64) : 100.0
                const_z_percent = haskey(row, :staticloadpercent) ? safe_get_value(row[:staticloadpercent], 100.0, Float64) : 100.0
                const_i_percent = 100.0 - const_z_percent - const_p_percent
                
                
                # Validate percentage sum is 100%
                total_percent = const_z_percent + const_i_percent + const_p_percent
                if abs(total_percent - 100.0) > 1e-6
                    @warn "Row $i: Load $name (ID: $index) has load model percentages that don't sum to 100% ($total_percent%), normalizing"
                    if total_percent > 0
                        const_z_percent = const_z_percent * 100.0 / total_percent
                        const_i_percent = const_i_percent * 100.0 / total_percent
                        const_p_percent = const_p_percent * 100.0 / total_percent
                    else
                        const_z_percent = 0.0
                        const_i_percent = 0.0
                        const_p_percent = 100.0
                    end
                end
                
                # scaling = haskey(row, :uniformscalepq) ? safe_get_value(row[:uniformscalepq], 1.0, Float64) : 1.0
                scaling = 1.0
                
                # Validate scaling factor is reasonable
                if scaling < 0.0
                    @warn "Row $i: Load $name (ID: $index) has invalid scaling factor ($scaling), setting to 1.0"
                    scaling = 1.0
                end
                
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Create Load object and add to case
                push!(case.loadsDC, LoadDC(index, name, bus, p_mw, const_z_percent, const_i_percent,
                                     const_p_percent, scaling, in_service))
                
                processed_rows += 1
                
            catch e
                @error "Error processing load data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Load data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading load data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_trafo!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load two-winding transformer data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing two-winding transformer data
"""
function load_trafo!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading two-winding transformer data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Two-winding transformer table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :frombus, :tobus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Two-winding transformer table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid transformer index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus names to bus IDs
                hv_bus_name = safe_get_value(row[:frombus], "", String)
                lv_bus_name = safe_get_value(row[:tobus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus names to integer IDs
                hv_bus = 0
                lv_bus = 0
                
                if haskey(case.bus_name_to_id, hv_bus_name)
                    hv_bus = case.bus_name_to_id[hv_bus_name]
                else
                    @warn "Row $i: Transformer $name (ID: $index) has HV bus name '$hv_bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.bus_name_to_id, lv_bus_name)
                    lv_bus = case.bus_name_to_id[lv_bus_name]
                else
                    @warn "Row $i: Transformer $name (ID: $index) has LV bus name '$lv_bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus indices are valid
                if hv_bus <= 0 || lv_bus <= 0
                    @warn "Row $i: Transformer $name (ID: $index) connects to invalid buses (HV: $hv_bus, LV: $lv_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate buses exist
                bus_indices = Set(bus.index for bus in case.busesAC)
                if !(hv_bus in bus_indices) || !(lv_bus in bus_indices)
                    @warn "Row $i: Transformer $name (ID: $index) connects to non-existent buses (HV: $hv_bus, LV: $lv_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate HV and LV buses are different
                if hv_bus == lv_bus
                    @warn "Row $i: Transformer $name (ID: $index) has same HV and LV buses ($hv_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                sn_mva = haskey(row, :ansimva) ? safe_get_value(row[:ansimva], 0.0, Float64)/1000 : 0.0
                
                # Validate rated capacity is reasonable
                if sn_mva <= 0.0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid rated capacity ($sn_mva MVA), using default value 10.0 MVA"
                    sn_mva = 10.0
                end
                
                vn_hv_kv = haskey(row, :primkv) ? safe_get_value(row[:primkv], 0.0, Float64) : 0.0
                vn_lv_kv = haskey(row, :seckv) ? safe_get_value(row[:seckv], 0.0, Float64) : 0.0
                
                # Validate rated voltages are reasonable
                                if vn_hv_kv <= 0.0
                    # Try to get rated voltage from bus data
                    hv_bus_obj = findfirst(b -> b.index == hv_bus, case.buses)
                    if hv_bus_obj !== nothing
                        vn_hv_kv = case.buses[hv_bus_obj].vn_kv
                        @warn "Row $i: Transformer $name (ID: $index) has invalid HV rated voltage, using bus voltage $vn_hv_kv kV"
                    else
                        vn_hv_kv = 110.0
                        @warn "Row $i: Transformer $name (ID: $index) has invalid HV rated voltage, using default value $vn_hv_kv kV"
                    end
                end
                
                if vn_lv_kv <= 0.0
                    # Try to get rated voltage from bus data
                    lv_bus_obj = findfirst(b -> b.index == lv_bus, case.buses)
                    if lv_bus_obj !== nothing
                        vn_lv_kv = case.buses[lv_bus_obj].vn_kv
                        @warn "Row $i: Transformer $name (ID: $index) has invalid LV rated voltage, using bus voltage $vn_lv_kv kV"
                    else
                        vn_lv_kv = 10.0
                        @warn "Row $i: Transformer $name (ID: $index) has invalid LV rated voltage, using default value $vn_lv_kv kV"
                    end
                end
                
                # Validate HV and LV voltage relationship
                if vn_hv_kv <= vn_lv_kv
                    @warn "Row $i: Transformer $name (ID: $index) has HV voltage ($vn_hv_kv kV) not greater than LV voltage ($vn_lv_kv kV), marked but continuing"
                end
                
                # Impedance parameters
                z_percent = haskey(row, :ansiposz) ? safe_get_value(row[:ansiposz], 0.0, Float64)/100 : 0.0
                x_r = haskey(row, :ansiposxr) ? safe_get_value(row[:ansiposxr], 0.0, Float64) : 0.0

                z0_percent = haskey(row, :ansizeroz) ? safe_get_value(row[:ansizeroz], 0.0, Float64)/100 : 0.0
                x0_r0 = haskey(row, :ansizeroxoverr) ? safe_get_value(row[:ansizeroxoverr], 0.0, Float64) : 0.0

                # Validate impedance parameters are reasonable
                if z_percent <= 0.0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid impedance percentage ($z_percent%), using default value 10.0%"
                    z_percent = 10.0/100
                end
                
                if x_r < 0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid short circuit resistance ($x_r%), using default value 20.0"
                    x_r = 20.0
                end
                
                if z0_percent <= 0.0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid zero sequence impedance ($z0_percent%), using default value 10.0"
                    z0_percent = 10.0/100
                end

                if x0_r0 < 0.0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid zero sequence resistance ($x0_r0%), using default value 20.0"
                    x0_r0 = 20.0
                end
                
                # Transformer taps
                tap_neutral = haskey(row, :centertap) ? safe_get_value(row[:centertap], 0.0, Float64) : 0.0
                prim_tap = haskey(row, :primpercenttap) ? safe_get_value(row[:primpercenttap], 0.0, Float64) : 0.0
                sec_tap = haskey(row, :secpercenttap) ? safe_get_value(row[:secpercenttap], 0.0, Float64) : 0.0

                prim_tap_min = haskey(row, :primminpercentfixedtap) ? safe_get_value(row[:primminpercentfixedtap], 0.0, Float64) : 0.0
                prim_tap_max = haskey(row, :primmaxpercentfixedtap) ? safe_get_value(row[:primmaxpercentfixedtap], 0.0, Float64) : 0.0

                sec_tap_min = haskey(row, :secminpercentfixedtap) ? safe_get_value(row[:secminpercentfixedtap], 0.0, Float64) : 0.0
                sec_tap_max = haskey(row, :secmaxpercentfixedtap) ? safe_get_value(row[:secmaxpercentfixedtap], 0.0, Float64) : 0.0
                
                
                # Validate transformer taps are reasonable
                if prim_tap > prim_tap_max
                    @warn "Row $i: Transformer $name (ID: $index) has tap exceeding upper limit, invalid ($prim_tap), setting to 0"
                    prim_tap = 0.0
                end
                
                if prim_tap < prim_tap_min
                    @warn "Row $i: Transformer $name (ID: $index) has tap below lower limit, invalid ($prim_tap), setting to 0"
                    prim_tap = 0.0
                end
                
                if sec_tap > sec_tap_max
                    @warn "Row $i: Transformer $name (ID: $index) has tap exceeding upper limit, invalid ($sec_tap), setting to 0"
                    sec_tap = 0.0
                end

                if sec_tap < sec_tap_min
                    @warn "Row $i: Transformer $name (ID: $index) has tap below lower limit, invalid ($sec_tap), setting to 0"
                    sec_tap = 0.0
                end

                # Phase shift angle
                phaseshifthl = haskey(row, :phaseshifthl) ? safe_get_value(row[:phaseshifthl], 0.0, Float64) : 0.0
                phaseshiftps = haskey(row, :phaseshiftps) ? safe_get_value(row[:phaseshiftps], 0.0, Float64) : 0.0

                # Connection type
                vectororwinding = haskey(row, :vectororwinding) ? safe_get_value(row[:vectororwinding], "", String) : ""
                primconnectionbutton = haskey(row, :primconnectionbutton) ? safe_get_value(row[:primconnectionbutton], "", String) : ""
                secconnectionbutton = haskey(row, :secconnectionbutton) ? safe_get_value(row[:secconnectionbutton], "", String) : ""
                primneutralconn = haskey(row, :primneutralconn) ? safe_get_value(row[:primneutralconn], "", String) : ""
                secneutralconn = haskey(row, :secneutralconn) ? safe_get_value(row[:secneutralconn], "", String) : ""

                # First check if it's vector winding mode
                if vectororwinding == "1"
                    # Initialize primary and secondary connection types
                    prim_conn_type = ""
                    sec_conn_type = ""
                    
                    # Determine primary connection type
                    if primconnectionbutton == "1"
                        prim_conn_type = "D"  # D-type connection
                    elseif primconnectionbutton == "0"
                        prim_conn_type = "Y"  # Y-type connection
                        # Check if primary neutral is grounded
                        if primneutralconn == "1"
                            prim_conn_type = "Yn"  # Grounded Y-type connection
                        end
                    end
                    
                    # Determine secondary connection type
                    if secconnectionbutton == "1"
                        sec_conn_type = "d"  # D-type connection
                    elseif secconnectionbutton == "0"
                        sec_conn_type = "y"  # Y-type connection
                        # Check if secondary neutral is grounded
                        if secneutralconn == "1"
                            sec_conn_type = "yn"  # Grounded Y-type connection
                        end
                    end
                    
                    # Combine to form final vector_group
                    vector_group = prim_conn_type * sec_conn_type
                else
                    # If not vector winding mode, can set default value or other processing
                    vector_group = ""
                end

                
                # Other parameters
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # Validate parallel count is reasonable
                if parallel <= 0
                    @warn "Row $i: Transformer $name (ID: $index) has invalid parallel count ($parallel), setting to 1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                
                # Create Transformer2Wetap object and add to case
                # Modify code creating Transformer2Wetap object in load_trafo! function
                push!(case.transformers_2w_etap, Transformer2Wetap(index, name, "", hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,z_percent,
                        x_r, z0_percent, x0_r0, Int(round(tap_neutral)), prim_tap=prim_tap, sec_tap=sec_tap, prim_tap_min=Int(round(prim_tap_min)), prim_tap_max=Int(round(prim_tap_max)), # Convert to Int
                        sec_tap_min=Int(round(sec_tap_min)),sec_tap_max=Int(round(sec_tap_max)), phaseshifthl=phaseshifthl, phaseshiftps=phaseshiftps, vector_group=vector_group, parallel=parallel, 
                        df=df, in_service=in_service))
                
                processed_rows += 1
                
            catch e
                @error "Error processing two-winding transformer data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Two-winding transformer data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading two-winding transformer data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_trafo3ws!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load three-winding transformer data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing three-winding transformer data
"""
function load_trafo3ws!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading three-winding transformer data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Three-winding transformer table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :hv_bus, :mv_bus, :lv_bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Three-winding transformer table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract basic field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid three-winding transformer index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                std_type = haskey(row, :std_type) ? safe_get_value(row[:std_type], "", String) : ""
                
                # Three winding bus connections
                hv_bus = safe_get_value(row[:hv_bus], 0, Int)
                mv_bus = safe_get_value(row[:mv_bus], 0, Int)
                lv_bus = safe_get_value(row[:lv_bus], 0, Int)
                
                # Validate bus indices are valid
                if hv_bus <= 0 || mv_bus <= 0 || lv_bus <= 0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) connects to invalid buses (HV: $hv_bus, MV: $mv_bus, LV: $lv_bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate buses exist
                bus_indices = Set(bus.index for bus in case.buses)
                if !(hv_bus in bus_indices) || !(mv_bus in bus_indices) || !(lv_bus in bus_indices)
                    @warn "Row $i: Three-winding transformer $name (ID: $index) connects to non-existent buses, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate three buses are different
                if hv_bus == mv_bus || hv_bus == lv_bus || mv_bus == lv_bus
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has windings connected to the same bus, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Rated capacity
                sn_hv_mva = haskey(row, :sn_hv_mva) ? safe_get_value(row[:sn_hv_mva], 0.0, Float64) : 0.0
                sn_mv_mva = haskey(row, :sn_mv_mva) ? safe_get_value(row[:sn_mv_mva], 0.0, Float64) : 0.0
                sn_lv_mva = haskey(row, :sn_lv_mva) ? safe_get_value(row[:sn_lv_mva], 0.0, Float64) : 0.0
                
                # Validate rated capacity is reasonable
                if sn_hv_mva <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV side rated capacity ($sn_hv_mva MVA), using default value 10.0 MVA"
                    sn_hv_mva = 10.0
                end
                
                if sn_mv_mva <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV side rated capacity ($sn_mv_mva MVA), using default value 10.0 MVA"
                    sn_mv_mva = 10.0
                end
                
                if sn_lv_mva <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV side rated capacity ($sn_lv_mva MVA), using default value 10.0 MVA"
                    sn_lv_mva = 10.0
                end
                
                # Rated voltage
                vn_hv_kv = haskey(row, :vn_hv_kv) ? safe_get_value(row[:vn_hv_kv], 0.0, Float64) : 0.0
                vn_mv_kv = haskey(row, :vn_mv_kv) ? safe_get_value(row[:vn_mv_kv], 0.0, Float64) : 0.0
                vn_lv_kv = haskey(row, :vn_lv_kv) ? safe_get_value(row[:vn_lv_kv], 0.0, Float64) : 0.0
                
                # Validate rated voltage is reasonable
                if vn_hv_kv <= 0.0
                    # Try to get rated voltage from bus data
                    hv_bus_obj = findfirst(b -> b.index == hv_bus, case.buses)
                    if hv_bus_obj !== nothing
                        vn_hv_kv = case.buses[hv_bus_obj].vn_kv
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV rated voltage, using bus voltage $vn_hv_kv kV"
                    else
                        vn_hv_kv = 110.0
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV rated voltage, using default value $vn_hv_kv kV"
                    end
                end
                
                if vn_mv_kv <= 0.0
                    # Try to get rated voltage from bus data
                    mv_bus_obj = findfirst(b -> b.index == mv_bus, case.buses)
                    if mv_bus_obj !== nothing
                        vn_mv_kv = case.buses[mv_bus_obj].vn_kv
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV rated voltage, using bus voltage $vn_mv_kv kV"
                    else
                        vn_mv_kv = 35.0
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV rated voltage, using default value $vn_mv_kv kV"
                    end
                end
                
                if vn_lv_kv <= 0.0
                    # Try to get rated voltage from bus data
                    lv_bus_obj = findfirst(b -> b.index == lv_bus, case.buses)
                    if lv_bus_obj !== nothing
                        vn_lv_kv = case.buses[lv_bus_obj].vn_kv
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV rated voltage, using bus voltage $vn_lv_kv kV"
                    else
                        vn_lv_kv = 10.0
                        @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV rated voltage, using default value $vn_lv_kv kV"
                    end
                end
                
                # Validate voltage level relationship
                if !(vn_hv_kv > vn_mv_kv && vn_mv_kv > vn_lv_kv)
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has incorrect voltage level relationship (HV: $vn_hv_kv kV, MV: $vn_mv_kv kV, LV: $vn_lv_kv kV), marked but continuing"
                end
                
                # Short circuit voltage percentage
                vk_hv_percent = haskey(row, :vk_hv_percent) ? safe_get_value(row[:vk_hv_percent], 0.0, Float64) : 0.0
                vk_mv_percent = haskey(row, :vk_mv_percent) ? safe_get_value(row[:vk_mv_percent], 0.0, Float64) : 0.0
                vk_lv_percent = haskey(row, :vk_lv_percent) ? safe_get_value(row[:vk_lv_percent], 0.0, Float64) : 0.0
                
                # Validate short circuit voltage is reasonable
                if vk_hv_percent <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV side short circuit voltage ($vk_hv_percent%), using default value 6.0%"
                    vk_hv_percent = 6.0
                end
                
                if vk_mv_percent <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV side short circuit voltage ($vk_mv_percent%), using default value 6.0%"
                    vk_mv_percent = 6.0
                end
                
                if vk_lv_percent <= 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV side short circuit voltage ($vk_lv_percent%), using default value 6.0%"
                    vk_lv_percent = 6.0
                end
                
                # Short circuit loss percentage
                vkr_hv_percent = haskey(row, :vkr_hv_percent) ? safe_get_value(row[:vkr_hv_percent], 0.0, Float64) : 0.0
                vkr_mv_percent = haskey(row, :vkr_mv_percent) ? safe_get_value(row[:vkr_mv_percent], 0.0, Float64) : 0.0
                vkr_lv_percent = haskey(row, :vkr_lv_percent) ? safe_get_value(row[:vkr_lv_percent], 0.0, Float64) : 0.0
                
                # Validate short circuit loss is reasonable
                if vkr_hv_percent < 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV side short circuit loss ($vkr_hv_percent%), setting to 0"
                    vkr_hv_percent = 0.0
                end
                
                if vkr_mv_percent < 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV side short circuit loss ($vkr_mv_percent%), setting to 0"
                    vkr_mv_percent = 0.0
                end
                
                if vkr_lv_percent < 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV side short circuit loss ($vkr_lv_percent%), setting to 0"
                    vkr_lv_percent = 0.0
                end
                
                # Validate short circuit impedance relationship
                if vkr_hv_percent > vk_hv_percent
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has HV side short circuit loss ($vkr_hv_percent%) greater than short circuit voltage ($vk_hv_percent%), adjusting to $(vk_hv_percent * 0.9)%"
                    vkr_hv_percent = vk_hv_percent * 0.9
                end
                
                if vkr_mv_percent > vk_mv_percent
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has MV side short circuit loss ($vkr_mv_percent%) greater than short circuit voltage ($vk_mv_percent%), adjusting to $(vk_mv_percent * 0.9)%"
                    vkr_mv_percent = vk_mv_percent * 0.9
                end
                
                if vkr_lv_percent > vk_lv_percent
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has LV side short circuit loss ($vkr_lv_percent%) greater than short circuit voltage ($vk_lv_percent%), adjusting to $(vk_lv_percent * 0.9)%"
                    vkr_lv_percent = vk_lv_percent * 0.9
                end
                
                # Iron loss and no-load current
                pfe_kw = haskey(row, :pfe_kw) ? safe_get_value(row[:pfe_kw], 0.0, Float64) : 0.0
                i0_percent = haskey(row, :i0_percent) ? safe_get_value(row[:i0_percent], 0.0, Float64) : 0.0
                
                # Validate iron loss and no-load current are reasonable
                if pfe_kw < 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid iron loss ($pfe_kw kW), setting to 0"
                    pfe_kw = 0.0
                end
                
                if i0_percent < 0.0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid no-load current ($i0_percent%), setting to 0"
                    i0_percent = 0.0
                end
                
                # Phase shift angles
                shift_mv_degree = haskey(row, :shift_mv_degree) ? safe_get_value(row[:shift_mv_degree], 0.0, Float64) : 0.0
                shift_lv_degree = haskey(row, :shift_lv_degree) ? safe_get_value(row[:shift_lv_degree], 0.0, Float64) : 0.0
                
                # Tap changer parameters
                tap_side = haskey(row, :tap_side) ? safe_get_value(row[:tap_side], "", String) : ""
                
                # Validate tap side is valid
                if !isempty(tap_side) && !(tap_side in ["hv", "mv", "lv"])
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid tap side ($tap_side), setting to hv"
                    tap_side = "hv"
                end
                
                tap_neutral = haskey(row, :tap_neutral) ? safe_get_value(row[:tap_neutral], 0, Int) : 0
                tap_min = haskey(row, :tap_min) ? safe_get_value(row[:tap_min], 0, Int) : 0
                tap_max = haskey(row, :tap_max) ? safe_get_value(row[:tap_max], 0, Int) : 0
                
                # Validate tap position range
                if tap_min > tap_max && tap_max != 0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid tap position range (min: $tap_min, max: $tap_max), swapping values"
                    tap_min, tap_max = tap_max, tap_min
                end
                
                tap_step_percent = haskey(row, :tap_step_percent) ? safe_get_value(row[:tap_step_percent], 0.0, Float64) : 0.0
                tap_step_degree = haskey(row, :tap_step_degree) ? safe_get_value(row[:tap_step_degree], 0.0, Float64) : 0.0
                tap_at_star_point = haskey(row, :tap_at_star_point) ? parse_bool(safe_get_value(row[:tap_at_star_point], false)) : false
                tap_pos = haskey(row, :tap_pos) ? safe_get_value(row[:tap_pos], 0, Int) : 0
                
                # Validate tap position is within range
                if tap_pos < tap_min && tap_min != 0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has current tap position ($tap_pos) less than minimum ($tap_min), setting to minimum"
                    tap_pos = tap_min
                elseif tap_pos > tap_max && tap_max != 0
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has current tap position ($tap_pos) greater than maximum ($tap_max), setting to maximum"
                    tap_pos = tap_max
                end
                
                # Operational status
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                
                # Technical parameters
                vector_group_hv_mv = haskey(row, :vector_group_hv_mv) ? safe_get_value(row[:vector_group_hv_mv], "", String) : ""
                vector_group_hv_lv = haskey(row, :vector_group_hv_lv) ? safe_get_value(row[:vector_group_hv_lv], "", String) : ""
                vector_group_mv_lv = haskey(row, :vector_group_mv_lv) ? safe_get_value(row[:vector_group_mv_lv], "", String) : ""
                hv_connection = haskey(row, :hv_connection) ? safe_get_value(row[:hv_connection], "", String) : ""
                mv_connection = haskey(row, :mv_connection) ? safe_get_value(row[:mv_connection], "", String) : ""
                lv_connection = haskey(row, :lv_connection) ? safe_get_value(row[:lv_connection], "", String) : ""
                thermal_capacity_mw = haskey(row, :thermal_capacity_mw) ? safe_get_value(row[:thermal_capacity_mw], 0.0, Float64) : 0.0
                cooling_type = haskey(row, :cooling_type) ? safe_get_value(row[:cooling_type], "", String) : ""
                oil_volume_liters = haskey(row, :oil_volume_liters) ? safe_get_value(row[:oil_volume_liters], 0.0, Float64) : 0.0
                winding_material = haskey(row, :winding_material) ? safe_get_value(row[:winding_material], "", String) : ""
                
                # Validate connection types are valid
                valid_connections = ["Y", "D", "Z", "y", "d", "z", ""]
                if !(hv_connection in valid_connections)
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid HV side connection type ($hv_connection), setting to Y"
                    hv_connection = "Y"
                end
                
                if !(mv_connection in valid_connections)
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid MV side connection type ($mv_connection), setting to Y"
                    mv_connection = "Y"
                end
                
                if !(lv_connection in valid_connections)
                    @warn "Row $i: Three-winding transformer $name (ID: $index) has invalid LV side connection type ($lv_connection), setting to Y"
                    lv_connection = "Y"
                end
                
                # Create dictionary for additional parameters
                kwargs = Dict{Symbol, Any}(
                    :tap_side => tap_side,
                    :tap_neutral => tap_neutral,
                    :tap_min => tap_min,
                    :tap_max => tap_max,
                    :tap_step_percent => tap_step_percent,
                    :tap_step_degree => tap_step_degree,
                    :tap_at_star_point => tap_at_star_point,
                    :tap_pos => tap_pos,
                    :in_service => in_service,
                    :vector_group_hv_mv => vector_group_hv_mv,
                    :vector_group_hv_lv => vector_group_hv_lv,
                    :vector_group_mv_lv => vector_group_mv_lv,
                    :hv_connection => hv_connection,
                    :mv_connection => mv_connection,
                    :lv_connection => lv_connection,
                    :thermal_capacity_mw => thermal_capacity_mw,
                    :cooling_type => cooling_type,
                    :oil_volume_liters => oil_volume_liters,
                    :winding_material => winding_material
                )
                
                # Create Transformer3W object and add to case
                try
                    trafo3w = Transformer3W(
                        index, name, std_type, hv_bus, mv_bus, lv_bus, 
                        sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                        vn_hv_kv, vn_mv_kv, vn_lv_kv,
                        vk_hv_percent, vk_mv_percent, vk_lv_percent,
                        vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                        pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; 
                        kwargs...
                    )
                    
                    push!(case.trafo3ws, trafo3w)
                    processed_rows += 1
                catch e
                    @error "Error creating three-winding transformer object" exception=(e, catch_backtrace()) transformer_data=(index=index, name=name)
                    error_rows += 1
                end
                
            catch e
                @error "Error processing three-winding transformer data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Three-winding transformer data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading three-winding transformer data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load generator data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing generator data
"""
function load_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading generator data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Generator table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :name, :bus, :p_mw]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Generator table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract basic field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid generator index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "Row $i: Generator $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: Generator $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if bus exists
                if !any(b -> b.index == bus, case.buses)
                    @warn "Row $i: Generator $name (ID: $index) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                p_mw = safe_get_value(row[:p_mw], 0.0, Float64)
                vm_pu = haskey(row, :vm_pu) ? safe_get_value(row[:vm_pu], 1.0, Float64) : 1.0
                
                # Validate voltage value is reasonable
                if vm_pu <= 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid voltage setpoint ($vm_pu p.u.), setting to default value 1.0 p.u."
                    vm_pu = 1.0
                end
                
                sn_mva = haskey(row, :sn_mva) ? safe_get_value(row[:sn_mva], 0.0, Float64) : 0.0
                
                # Validate rated capacity is reasonable
                if sn_mva < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid rated capacity ($sn_mva MVA), setting to 0 MVA"
                    sn_mva = 0.0
                end
                
                min_p_mw = haskey(row, :min_p_mw) ? safe_get_value(row[:min_p_mw], 0.0, Float64) : 0.0
                max_p_mw = haskey(row, :max_p_mw) ? safe_get_value(row[:max_p_mw], 0.0, Float64) : 0.0
                min_q_mvar = haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0
                max_q_mvar = haskey(row, :max_q_mvar) ? safe_get_value(row[:max_q_mvar], 0.0, Float64) : 0.0
                
                # Validate power limits are reasonable
                if max_p_mw < min_p_mw && max_p_mw != 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid active power limits (min: $min_p_mw MW, max: $max_p_mw MW), swapping values"
                    max_p_mw, min_p_mw = min_p_mw, max_p_mw
                end
                
                if max_q_mvar < min_q_mvar && max_q_mvar != 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid reactive power limits (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar), swapping values"
                    max_q_mvar, min_q_mvar = min_q_mvar, max_q_mvar
                end
                
                scaling = haskey(row, :scaling) ? safe_get_value(row[:scaling], 1.0, Float64) : 1.0
                
                # Validate scaling factor is reasonable
                if scaling < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid scaling factor ($scaling), setting to default value 1.0"
                    scaling = 1.0
                end
                
                slack = haskey(row, :slack) ? parse_bool(safe_get_value(row[:slack], false)) : false
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # Economic parameters
                cost_a = haskey(row, :cost_a) ? safe_get_value(row[:cost_a], 0.0, Float64) : 0.0
                cost_b = haskey(row, :cost_b) ? safe_get_value(row[:cost_b], 0.0, Float64) : 0.0
                                cost_c = haskey(row, :cost_c) ? safe_get_value(row[:cost_c], 0.0, Float64) : 0.0
                startup_cost = haskey(row, :startup_cost) ? safe_get_value(row[:startup_cost], 0.0, Float64) : 0.0
                shutdown_cost = haskey(row, :shutdown_cost) ? safe_get_value(row[:shutdown_cost], 0.0, Float64) : 0.0
                
                # Validate cost coefficients are reasonable
                if cost_a < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid quadratic cost coefficient ($cost_a), setting to 0"
                    cost_a = 0.0
                end
                
                if cost_b < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid linear cost coefficient ($cost_b), setting to 0"
                    cost_b = 0.0
                end
                
                # Technical parameters
                min_up_time = haskey(row, :min_up_time) ? safe_get_value(row[:min_up_time], 0, Int) : 0
                min_down_time = haskey(row, :min_down_time) ? safe_get_value(row[:min_down_time], 0, Int) : 0
                ramp_up_mw_per_min = haskey(row, :ramp_up_mw_per_min) ? safe_get_value(row[:ramp_up_mw_per_min], 0.0, Float64) : 0.0
                ramp_down_mw_per_min = haskey(row, :ramp_down_mw_per_min) ? safe_get_value(row[:ramp_down_mw_per_min], 0.0, Float64) : 0.0
                startup_time = haskey(row, :startup_time) ? safe_get_value(row[:startup_time], 0, Int) : 0
                
                # Validate technical parameters are reasonable
                if min_up_time < 0
                    @warn "Row $i: Generator $name (ID: $index) has invalid minimum up time ($min_up_time), setting to 0"
                    min_up_time = 0
                end
                
                if min_down_time < 0
                    @warn "Row $i: Generator $name (ID: $index) has invalid minimum down time ($min_down_time), setting to 0"
                    min_down_time = 0
                end
                
                if ramp_up_mw_per_min < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid ramp up rate ($ramp_up_mw_per_min MW/min), setting to 0"
                    ramp_up_mw_per_min = 0.0
                end
                
                if ramp_down_mw_per_min < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid ramp down rate ($ramp_down_mw_per_min MW/min), setting to 0"
                    ramp_down_mw_per_min = 0.0
                end
                
                if startup_time < 0
                    @warn "Row $i: Generator $name (ID: $index) has invalid startup time ($startup_time), setting to 0"
                    startup_time = 0
                end
                
                # Reliability parameters
                mtbf_hours = haskey(row, :mtbf_hours) ? safe_get_value(row[:mtbf_hours], 0.0, Float64) : 0.0
                mttr_hours = haskey(row, :mttr_hours) ? safe_get_value(row[:mttr_hours], 0.0, Float64) : 0.0
                failure_rate_per_year = haskey(row, :failure_rate_per_year) ? safe_get_value(row[:failure_rate_per_year], 0.0, Float64) : 0.0
                planned_outage_hours_per_year = haskey(row, :planned_outage_hours_per_year) ? safe_get_value(row[:planned_outage_hours_per_year], 0.0, Float64) : 0.0
                forced_outage_rate = haskey(row, :forced_outage_rate) ? safe_get_value(row[:forced_outage_rate], 0.0, Float64) : 0.0
                
                # Validate reliability parameters are reasonable
                if mtbf_hours < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid mean time between failures ($mtbf_hours hours), setting to 0"
                    mtbf_hours = 0.0
                end
                
                if mttr_hours < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid mean time to repair ($mttr_hours hours), setting to 0"
                    mttr_hours = 0.0
                end
                
                if failure_rate_per_year < 0.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid annual failure rate ($failure_rate_per_year), setting to 0"
                    failure_rate_per_year = 0.0
                end
                
                if forced_outage_rate < 0.0 || forced_outage_rate > 1.0
                    @warn "Row $i: Generator $name (ID: $index) has invalid forced outage rate ($forced_outage_rate), setting to 0"
                    forced_outage_rate = 0.0
                end
                
                # Create Generator object and add to case
                generator = Dict(
                    :index => index,
                    :name => name,
                    :bus => bus,
                    :p_mw => p_mw,
                    :vm_pu => vm_pu,
                    :sn_mva => sn_mva,
                    :min_p_mw => min_p_mw,
                    :max_p_mw => max_p_mw,
                    :min_q_mvar => min_q_mvar,
                    :max_q_mvar => max_q_mvar,
                    :scaling => scaling,
                    :slack => slack,
                    :in_service => in_service,
                    :type => type,
                    :controllable => controllable,
                    :cost_a => cost_a,
                    :cost_b => cost_b,
                    :cost_c => cost_c,
                    :startup_cost => startup_cost,
                    :shutdown_cost => shutdown_cost,
                    :min_up_time => min_up_time,
                    :min_down_time => min_down_time,
                    :ramp_up_mw_per_min => ramp_up_mw_per_min,
                    :ramp_down_mw_per_min => ramp_down_mw_per_min,
                    :startup_time => startup_time,
                    :mtbf_hours => mtbf_hours,
                    :mttr_hours => mttr_hours,
                    :failure_rate_per_year => failure_rate_per_year,
                    :planned_outage_hours_per_year => planned_outage_hours_per_year,
                    :forced_outage_rate => forced_outage_rate
                )
                
                push!(case.generators, generator)
                processed_rows += 1
                
            catch e
                @error "Error processing generator data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Generator data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
        # Check if there are slack buses
        if any(gen -> gen[:slack], case.generators)
            @info "System has $(count(gen -> gen[:slack], case.generators)) slack buses"
        else
            @warn "System has no slack bus, which may cause power flow calculation to fail to converge"
        end
        
    catch e
        @error "Error loading generator data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_storages!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load energy storage device data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing energy storage device data
"""
function load_storages!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading energy storage device data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Storage table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Storage table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid storage device index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus = case.busdc_name_to_id[bus_name]
                else
                    @warn "Row $i: Storage device $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: Storage device $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if bus exists
                if !any(b -> b.index == bus, case.busesDC)
                    @warn "Row $i: Storage device $name (ID: $index) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # ETAP storage device specific parameters
                ra = haskey(row, :ra) ? safe_get_value(row[:ra], 0.0, Float64) : 0.0025
                cell = haskey(row, :nrofcells) ? safe_get_value(row[:nrofcells], 0.0, Float64) : 0.0
                str = haskey(row, :nrofstrings) ? safe_get_value(row[:nrofstrings], 0.0, Float64) : 0.0
                package = haskey(row, :noofpacks) ? safe_get_value(row[:noofpacks], 0.0, Float64) : 0.0
                
                # Validate parameters are reasonable
                if ra < 0.0
                    @warn "Row $i: Storage device $name (ID: $index) has invalid ra parameter ($ra), setting to default value 0.0"
                    ra = 0.0
                end
                
                if cell < 0.0
                    @warn "Row $i: Storage device $name (ID: $index) has invalid cell parameter ($cell), setting to default value 0.0"
                    cell = 0.0
                end
                
                if str < 0.0
                    @warn "Row $i: Storage device $name (ID: $index) has invalid string parameter ($str), setting to default value 0.0"
                    str = 0.0
                end
                
                if package < 0.0
                    @warn "Row $i: Storage device $name (ID: $index) has invalid package parameter ($package), setting to default value 0.0"
                    package = 0.0
                end

                vpc = haskey(row, :vpc) ? safe_get_value(row[:vpc], 0.0, Float64) : 2.06
                voc = vpc * cell *package/1000
                if voc < 0.0
                    @warn "Row $i: Storage device $name (ID: $index) has invalid voc parameter ($voc), setting to default value 0.0"
                    voc = 0.0
                end
                
                # Other parameters
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # Create Storageetap object and add to case
                push!(case.storageetap, Storageetap(
                    index, name, bus, ra, cell, str, package, voc, in_service, type, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "Error processing storage device data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Storage device data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading storage device data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_converters!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load converter data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing converter data
"""
function load_converters!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading converter data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Converter table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :busid, :cznetwork]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Converter table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid converter index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                parent_pv = haskey(row, :parentname) ? safe_get_value(row[:parentname], "", String) : ""
                if parent_pv !== ""
                    @warn "Row $i: Converter $name (ID: $index) is connected to an AC PV system, will be processed during PV system loading"
                    continue
                end

                # Map bus names to bus IDs
                bus_ac_name = safe_get_value(row[:busid], "", String)
                bus_dc_name = safe_get_value(row[:cznetwork], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus names to integer IDs
                bus_ac = 0
                bus_dc = 0
                
                if haskey(case.bus_name_to_id, bus_ac_name)
                    bus_ac = case.bus_name_to_id[bus_ac_name]
                else
                    @warn "Row $i: Converter $name (ID: $index) has AC side bus name '$bus_ac_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.busdc_name_to_id, bus_dc_name)
                    bus_dc = case.busdc_name_to_id[bus_dc_name]
                else
                    @warn "Row $i: Converter $name (ID: $index) has DC side bus name '$bus_dc_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus indices are valid
                if bus_ac <= 0 || bus_dc <= 0
                    @warn "Row $i: Converter $name (ID: $index) connects to invalid buses (AC: $bus_ac, DC: $bus_dc), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract other field values from row data
                p_mw = safe_get_value(row[:gencat0ackw], 0.0, Float64)/1000
                q_mvar = (haskey(row, :gencat0kvar) ? safe_get_value(row[:gencat0kvar], 0.0, Float64) : 0.0)/1000
                vm_ac_pu = haskey(row, :vm_ac_pu) ? safe_get_value(row[:vm_ac_pu], 1.0, Float64) : 1.0
                vm_dc_pu = haskey(row, :vm_dc_pu) ? safe_get_value(row[:vm_dc_pu], 1.0, Float64) : 1.0
                
                # Validate voltage values are reasonable
                if vm_ac_pu <= 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid AC side voltage ($vm_ac_pu p.u.), setting to default value 1.0 p.u."
                    vm_ac_pu = 1.0
                end
                
                if vm_dc_pu <= 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid DC side voltage ($vm_dc_pu p.u.), setting to default value 1.0 p.u."
                    vm_dc_pu = 1.0
                end
                
                loss_percent = 1-(haskey(row, :dcpercenteff) ? safe_get_value(row[:dcpercenteff], 0.0, Float64) : 0.0)/100
                loss_mw = haskey(row, :loss_mw) ? safe_get_value(row[:loss_mw], 0.0, Float64) : 0.0
                
                # Validate losses are reasonable
                if loss_percent < 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid loss percentage ($loss_percent%), setting to 0%"
                    loss_percent = 0.0
                end
                
                if loss_mw < 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid loss power ($loss_mw MW), setting to 0 MW"
                    loss_mw = 0.0
                end
                
                max_p_mw = (haskey(row, :kwmax) ? safe_get_value(row[:kwmax], 0.0, Float64) : 0.0)/1000
                min_p_mw = (haskey(row, :kwmin) ? safe_get_value(row[:kwmin], 0.0, Float64) : 0.0)/1000
                max_q_mvar = (haskey(row, :kvarmax) ? safe_get_value(row[:kvarmax], 0.0, Float64) : 0.0)/1000
                min_q_mvar = (haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0)/1000
                
                # Validate power limits are reasonable
                if min_p_mw > max_p_mw && max_p_mw != 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid active power limits (min: $min_p_mw MW, max: $max_p_mw MW), swapping values"
                    min_p_mw, max_p_mw = max_p_mw, min_p_mw
                end
                
                if min_q_mvar > max_q_mvar && max_q_mvar != 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid reactive power limits (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar), swapping values"
                    min_q_mvar, max_q_mvar = max_q_mvar, min_q_mvar
                end
                
                control_mode = haskey(row, :acoperationmode) ? safe_get_value(row[:acoperationmode], "", String) : ""
                droop_kv = haskey(row, :droop_kv) ? safe_get_value(row[:droop_kv], 0.0, Float64) : 0.05
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # Create Converter object and add to case
                push!(case.converters, Converter(
                    index, name, bus_ac, bus_dc, p_mw, q_mvar, 
                    vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                    max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                    control_mode, droop_kv, in_service, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "Error processing converter data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Converter data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading converter data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_virtual_power_plants!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load virtual power plant (VPP) data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing virtual power plant data
"""
function load_virtual_power_plants!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    # Use DataFrame processing
    df = DataFrame(XLSX.readtable(file_path, sheet_name))
    
    # Ensure data is not empty
    if isempty(df)
        @info "Virtual power plant table is empty"
        return
    end
    
    # Convert column names to lowercase
    rename!(df, lowercase.(names(df)))
    
    # Iterate through each row of data
    for row in eachrow(df)
        try
            # Extract basic field values from row data
            index = safe_get_value(row[:index], 0, Int)
            name = safe_get_value(row[:name], "", String)
            description = haskey(row, :description) ? safe_get_value(row[:description], "", String) : ""
            control_area = haskey(row, :control_area) ? safe_get_value(row[:control_area], "", String) : ""
            
            # Capacity and energy parameters
            capacity_mw = haskey(row, :capacity_mw) ? safe_get_value(row[:capacity_mw], 0.0, Float64) : 0.0
            energy_mwh = haskey(row, :energy_mwh) ? safe_get_value(row[:energy_mwh], 0.0, Float64) : 0.0
            
            # Response and ramp parameters
            response_time_s = haskey(row, :response_time_s) ? safe_get_value(row[:response_time_s], 0.0, Float64) : 0.0
            ramp_rate_mw_per_min = haskey(row, :ramp_rate_mw_per_min) ? safe_get_value(row[:ramp_rate_mw_per_min], 0.0, Float64) : 0.0
            availability_percent = haskey(row, :availability_percent) ? safe_get_value(row[:availability_percent], 100.0, Float64) : 100.0
            
            # Operational information
            operator = haskey(row, :operator) ? safe_get_value(row[:operator], "", String) : ""
            in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
            
            # Collect additional kwargs parameters
            kwargs = Dict{Symbol, Any}()
            
            # Resource information
            if haskey(row, :resource_type)
                kwargs[:resource_type] = safe_get_value(row[:resource_type], "", String)
            end
            
            if haskey(row, :resource_id)
                kwargs[:resource_id] = safe_get_value(row[:resource_id], 0, Int)
            end
            
            if haskey(row, :capacity_share_percent)
                kwargs[:capacity_share_percent] = safe_get_value(row[:capacity_share_percent], 0.0, Float64)
            end
            
            if haskey(row, :control_priority)
                kwargs[:control_priority] = safe_get_value(row[:control_priority], 0, Int)
            end
            
            if haskey(row, :resource_response_time_s)
                kwargs[:resource_response_time_s] = safe_get_value(row[:resource_response_time_s], 0.0, Float64)
            end
            
            if haskey(row, :max_duration_h)
                kwargs[:max_duration_h] = safe_get_value(row[:max_duration_h], 0.0, Float64)
            end
            
            # Load information
            if haskey(row, :timestamp) && !ismissing(row[:timestamp])
                # Process timestamp, adjust according to actual format
                if isa(row[:timestamp], String)
                    kwargs[:timestamp] = DateTime(row[:timestamp], dateformat"yyyy-mm-dd HH:MM:SS")
                elseif isa(row[:timestamp], DateTime)
                    kwargs[:timestamp] = row[:timestamp]
                else
                    kwargs[:timestamp] = DateTime(now())
                end
            end
            
            if haskey(row, :p_mw)
                kwargs[:p_mw] = safe_get_value(row[:p_mw], 0.0, Float64)
            end
            
            if haskey(row, :q_mvar)
                kwargs[:q_mvar] = safe_get_value(row[:q_mvar], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_up_mw)
                kwargs[:flexibility_up_mw] = safe_get_value(row[:flexibility_up_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_down_mw)
                kwargs[:flexibility_down_mw] = safe_get_value(row[:flexibility_down_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_duration_h)
                kwargs[:flexibility_duration_h] = safe_get_value(row[:flexibility_duration_h], 0.0, Float64)
            end
            
            # Create VirtualPowerPlant object and add to case
            vpp = VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                                   response_time_s, ramp_rate_mw_per_min, availability_percent,
                                   operator, in_service; kwargs...)
            
            push!(case.virtual_power_plants, vpp)
            
        catch e
            @warn "Error processing virtual power plant: $e"
            @warn "Problem row: $(row)"
        end
    end
    
    @info "Loaded $(length(case.virtual_power_plants)) virtual power plants"
end


"""
    load_ext_grids!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load external grid data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing external grid data
"""
function load_ext_grids!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading external grid data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "External grid table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "External grid table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid external grid index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "Row $i: External grid $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus <= 0
                    @warn "Row $i: External grid $name (ID: $index) connects to invalid bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if bus exists
                if !any(b -> b.bus_id == bus, case.busesAC)
                    @warn "Row $i: External grid $name (ID: $bus_id) connects to non-existent bus ($bus), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Voltage parameters
                vm_pu = haskey(row, :vm_pu) ? safe_get_value(row[:vm_pu], 1.0, Float64) : 1.0
                va_degree = haskey(row, :va_degree) ? safe_get_value(row[:va_degree], 0.0, Float64) : 0.0
                
                # Validate voltage value is reasonable
                if vm_pu <= 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid voltage setpoint ($vm_pu p.u.), setting to default value 1.0 p.u."
                    vm_pu = 1.0
                end
                
                # Operational status
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Short circuit capacity parameters
                s_sc_max_mva = haskey(row, :s_sc_max_mva) ? safe_get_value(row[:s_sc_max_mva], 1000.0, Float64) : 1000.0
                s_sc_min_mva = haskey(row, :s_sc_min_mva) ? safe_get_value(row[:s_sc_min_mva], 1000.0, Float64) : 1000.0
                
                # Validate short circuit capacity is reasonable
                if s_sc_max_mva < s_sc_min_mva && s_sc_min_mva > 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid short circuit capacity range (min: $s_sc_min_mva MVA, max: $s_sc_max_mva MVA), swapping values"
                    s_sc_max_mva, s_sc_min_mva = s_sc_min_mva, s_sc_max_mva
                end
                
                if s_sc_max_mva <= 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid maximum short circuit capacity ($s_sc_max_mva MVA), setting to default value 1000.0 MVA"
                    s_sc_max_mva = 1000.0
                end
                
                if s_sc_min_mva <= 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid minimum short circuit capacity ($s_sc_min_mva MVA), setting to default value 1000.0 MVA"
                    s_sc_min_mva = 1000.0
                end
                
                # Impedance ratio parameters
                posr = haskey(row, :posr) ? safe_get_value(row[:posr], 0.1, Float64) : 0.1
                posx = haskey(row, :posx) ? safe_get_value(row[:posx], 1.0, Float64) : 1.0
                zeror = haskey(row, :zeror) ? safe_get_value(row[:zeror], 0.1, Float64) : 0.1
                zerox = haskey(row, :zerox) ? safe_get_value(row[:zerox], 1.0, Float64) : 1.0
                rx_max = posr / posx
                rx_min = posr / posx
                r0x0_max = zeror / zerox
                x0x_max = zerox / posx
                
                # Validate impedance ratios are reasonable
                if rx_max < 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid maximum R/X ratio ($rx_max), setting to default value 0.1"
                    rx_max = 0.1
                end
                
                if rx_min < 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid minimum R/X ratio ($rx_min), setting to default value 0.1"
                    rx_min = 0.1
                end
                
                if rx_max < rx_min && rx_min > 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid R/X ratio range (min: $rx_min, max: $rx_max), swapping values"
                    rx_max, rx_min = rx_min, rx_max
                end
                
                if r0x0_max < 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid zero sequence R0/X0 ratio ($r0x0_max), setting to default value 0.1"
                    r0x0_max = 0.1
                end
                
                if x0x_max < 0.0
                    @warn "Row $i: External grid $name (ID: $index) has invalid zero sequence X0/X ratio ($x0x_max), setting to default value 1.0"
                    x0x_max = 1.0
                end
                
                # Controllability
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], false)) : false
                
                # Create ExternalGrid object and add to case
                push!(case.ext_grids, ExternalGrid(
                    index, name, bus, vm_pu, va_degree, in_service,
                    s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "Error processing external grid data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "External grid data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading external grid data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end
"""
    load_switches!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load switch data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing switch data
"""
function load_switches!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading switch data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "Switch table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :fromelement, :toelement]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "Switch table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid switch index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # Map bus name to bus ID - from bus
                bus_from_name = safe_get_value(row[:fromelement], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus_from = 0
                
                if haskey(case.bus_name_to_id, bus_from_name)
                    bus_from = case.bus_name_to_id[bus_from_name]
                else
                    @warn "Row $i: Switch $name (ID: $index) has from bus name '$bus_from_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Map bus name to bus ID - to bus
                bus_to_name = safe_get_value(row[:toelement], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus_to = 0
                
                if haskey(case.bus_name_to_id, bus_to_name)
                    bus_to = case.bus_name_to_id[bus_to_name]
                else
                    @warn "Row $i: Switch $name (ID: $index) has to bus name '$bus_to_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus indices are valid
                if bus_from <= 0
                    @warn "Row $i: Switch $name (ID: $index) connects to invalid from bus ($bus_from), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if bus_to <= 0
                    @warn "Row $i: Switch $name (ID: $index) connects to invalid to bus ($bus_to), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if buses exist
                if !any(b -> b.index == bus_from, case.busesAC)
                    @warn "Row $i: Switch $name (ID: $index) connects to non-existent from bus ($bus_from), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if !any(b -> b.index == bus_to, case.busesAC)
                    @warn "Row $i: Switch $name (ID: $index) connects to non-existent to bus ($bus_to), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate from and to buses are different
                if bus_from == bus_to
                    @warn "Row $i: Switch $name (ID: $index) has the same from and to buses ($bus_from), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Element type and ID
                element_type = haskey(row, :element_type) ? safe_get_value(row[:element_type], "l", String) : "l"
                element_id = haskey(row, :element_id) ? safe_get_value(row[:element_id], 0, Int) : 0
                
                # Validate element type is valid
                valid_element_types = ["l", "t", "t3", "b"]
                if !(lowercase(element_type) in valid_element_types)
                    @warn "Row $i: Switch $name (ID: $index) has invalid element type ($element_type), setting to default value 'l'"
                    element_type = "l"
                end
                
                # Switch status
                closed = haskey(row, :pdestatus) ? parse_bool(safe_get_value(row[:pdestatus], true)=="Closed") : true
                
                # Switch type and parameters
                type = haskey(row, :type) ? safe_get_value(row[:type], "CB", String) : "CB"
                z_ohm = haskey(row, :z_ohm) ? safe_get_value(row[:z_ohm], 0.0, Float64) : 0.0
                
                # Validate impedance value is reasonable
                if z_ohm < 0.0
                    @warn "Row $i: Switch $name (ID: $index) has invalid impedance value ($z_ohm Ω), setting to default value 0.0 Ω"
                    z_ohm = 0.0
                end
                
                # Operational status
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Create Switch object and add to case
                push!(case.hvcbs, HighVoltageCircuitBreaker(
                    index, name, bus_from, bus_to, element_type, element_id,
                    closed, type, z_ohm, in_service
                ))
                
                # Update index mapping
                # case.switch_indices[index] = length(case.hvcbs)
                
                processed_rows += 1
                
            catch e
                @error "Error processing switch data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "Switch data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading switch data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_pv_arrays!(case::JuliaPowerCase, file_path::String, sheet_name::String)
Load PV array data from Excel file and add to power system case.
Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing PV array data
"""
function load_pv_arrays!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading PV array data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # Ensure data is not empty
        if isempty(df)
            @info "PV array table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        
        # Verify necessary columns exist
        required_columns = [:index, :id, :bus, :numpanelseries, :numpanelparallel, :vmpp, :impp, :voc, :isc, :pvanumcells]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "PV array table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid PV array index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)

                inverterinclude = haskey(row, :inverterincluded) ? safe_get_value(row[:inverterincluded], "", String) : "0"
                if  inverterinclude == "1"
                    @warn "Row $i: PV array $name (ID: $index) is used for building AC PV system, skipping this row"
                    continue
                end
                
                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus_id = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus_id = case.busdc_name_to_id[bus_name]
                else
                    @warn "Row $i: PV array $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus_id <= 0
                    @warn "Row $i: PV array $name (ID: $index) connects to invalid bus ($bus_id), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if bus exists
                if !any(b -> b.index == bus_id, case.busesDC)
                    @warn "Row $i: PV array $name (ID: $index) connects to non-existent bus ($bus_id), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Extract PV array parameters
                numpanelseries = safe_get_value(row[:numpanelseries], 0, Int)
                numpanelparallel = safe_get_value(row[:numpanelparallel], 0, Int)
                vmpp = safe_get_value(row[:vmpp], 0.0, Float64)
                impp = safe_get_value(row[:impp], 0.0, Float64)
                voc = safe_get_value(row[:voc], 0.0, Float64)
                isc = safe_get_value(row[:isc], 0.0, Float64)
                pvanumcells = safe_get_value(row[:pvanumcells], 0, Int)
                
                # Validate parameter values are valid
                if numpanelseries <= 0
                    @warn "Row $i: PV array $name (ID: $index) has invalid number of series panels ($numpanelseries), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if numpanelparallel <= 0
                    @warn "Row $i: PV array $name (ID: $index) has invalid number of parallel panels ($numpanelparallel), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if vmpp <= 0.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid maximum power point voltage ($vmpp V), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if impp <= 0.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid maximum power point current ($impp A), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if voc <= 0.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid open circuit voltage ($voc V), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if isc <= 0.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid short circuit current ($isc A), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if pvanumcells <= 0
                    @warn "Row $i: PV array $name (ID: $index) has invalid number of cells ($pvanumcells), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Temperature parameters
                temperature = haskey(row, :temperature) ? safe_get_value(row[:temperature], 25.0, Float64) : 25.0
                # Validate temperature value is reasonable
                if temperature < -40.0 || temperature > 85.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid temperature value ($temperature °C), setting to default value 25.0 °C"
                    temperature = 25.0
                end

                # Irradiance parameters
                irradiance = haskey(row, :irradiance) ? safe_get_value(row[:irradiance], 1000.0, Float64) : 1000.0
                # Validate irradiance value is reasonable
                if irradiance < 0.0 || irradiance > 2000.0
                    @warn "Row $i: PV array $name (ID: $index) has invalid irradiance value ($irradiance W/m²), setting to default value 1000.0 W/m²"
                    irradiance = 1000.0
                end

                # Additional parameters
                α_isc = haskey(row, :aisctemp) ? safe_get_value(row[:aisctemp], 0.0, Float64) : 0.0

                β_voc = haskey(row, :bvoctemp) ? safe_get_value(row[:bvoctemp], 0.0, Float64) : 0.0
                # Operational status
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Create PVArray object and add to case
                push!(case.pvarray, PVArray(
                    index, name, bus_id, numpanelseries, numpanelparallel,
                    vmpp, impp, voc, isc, pvanumcells,temperature,irradiance, α_isc, β_voc,
                    in_service
                ))
                
                # Update index mapping (if needed)
                # case.pvarray_indices[index] = length(case.pvarrays)
                
                processed_rows += 1
                
            catch e
                @error "Error processing PV array data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "PV array data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading PV array data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_ac_pv_system!(case::JuliaPowerCase, file_path::String, sheet_name::String)
Load AC PV system data from Excel file and add to power system case.
Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing AC PV system data
"""
function load_ac_pv_system!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # Use DataFrame processing
        @info "Reading AC PV system data..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        inverter = DataFrame(XLSX.readtable(file_path, "inverter"))
        inverter = filter(row -> !ismissing(row[:ParentName]), inverter)
        
        # Ensure data is not empty
        if isempty(df)
            @info "AC PV system table is empty"
            return
        end
        
        # Convert column names to lowercase
        rename!(df, lowercase.(names(df)))
        rename!(inverter, lowercase.(names(inverter)))
        
        # Verify necessary columns exist
        required_columns = [:index, :id, :bus, :numpanelseries, :numpanelparallel, :vmpp, :impp, :voc, :isc, :pvanumcells]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "AC PV system table missing required columns: $(join(missing_columns, ", "))"
            return
        end
        
        # Record processed rows and error rows
        processed_rows = 0
        error_rows = 0
        
        # Iterate through each row of data
        for (i, row) in enumerate(eachrow(df))
            try
                # Extract field values from row data
                index = safe_get_value(row[:index], 0, Int)
                
                # Validate index is valid
                if index <= 0
                    @warn "Row $i: Invalid AC PV system index ($index), skipping this row"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                inverterinclude = haskey(row, :inverterincluded) ? safe_get_value(row[:inverterincluded], "", String) : "0"
                if  inverterinclude == "0"
                    # @warn "Row $i: AC PV system $name (ID: $index) PV array does not include inverter, skipping this row"
                    continue
                end

                # Map bus name to bus ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # Use case.bus_name_to_id dictionary to convert bus name to integer ID
                bus_id = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus_id = case.bus_name_to_id[bus_name]
                else
                    @warn "Row $i: AC PV system $name (ID: $index) has bus name '$bus_name' that doesn't exist in bus_name_to_id dictionary, skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Validate bus index is valid
                if bus_id <= 0
                    @warn "Row $i: AC PV system $name (ID: $index) connects to invalid bus ($bus_id), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Check if bus exists
                if !any(b -> b.index == bus_id, case.busesAC)
                    @warn "Row $i: AC PV system $name (ID: $index) connects to non-existent bus ($bus_id), skipping this row"
                    error_rows += 1
                    continue
                end

                # Extract corresponding inverter
                inverter_row = inverter[inverter.parentname .== name, :][1,:]
                
                # Extract power parameters
                p_mw = safe_get_value(inverter_row[:gencat0ackw], 0.0, Float64)/1000
                q_mvar = (haskey(inverter_row, :gencat0kvar) ? safe_get_value(inverter_row[:gencat0kvar], 0.0, Float64) : 0.0)/1000
                vm_ac_pu = haskey(inverter_row, :vm_ac_pu) ? safe_get_value(inverter_row[:vm_ac_pu], 1.0, Float64) : 1.0
                vm_dc_pu = haskey(inverter_row, :vm_dc_pu) ? safe_get_value(inverter_row[:vm_dc_pu], 1.0, Float64) : 1.0

                # Validate voltage values are reasonable
                if vm_ac_pu <= 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid AC side voltage ($vm_ac_pu p.u.), setting to default value 1.0 p.u."
                    vm_ac_pu = 1.0
                end
                
                if vm_dc_pu <= 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid DC side voltage ($vm_dc_pu p.u.), setting to default value 1.0 p.u."
                    vm_dc_pu = 1.0
                end
                
                loss_percent = 1-(haskey(inverter_row, :dcpercenteff) ? safe_get_value(inverter_row[:dcpercenteff], 0.0, Float64) : 0.0)/100
                loss_mw = haskey(inverter_row, :loss_mw) ? safe_get_value(inverter_row[:loss_mw], 0.0, Float64) : 0.0
                
                # Validate losses are reasonable
                if loss_percent < 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid loss percentage ($loss_percent%), setting to 0%"
                    loss_percent = 0.0
                end
                
                if loss_mw < 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid loss power ($loss_mw MW), setting to 0 MW"
                    loss_mw = 0.0
                end
                
                max_p_mw = (haskey(inverter_row, :kwmax) ? safe_get_value(inverter_row[:kwmax], 0.0, Float64) : 0.0)/1000
                min_p_mw = (haskey(inverter_row, :kwmin) ? safe_get_value(inverter_row[:kwmin], 0.0, Float64) : 0.0)/1000
                max_q_mvar = (haskey(inverter_row, :kvarmax) ? safe_get_value(inverter_row[:kvarmax], 0.0, Float64) : 0.0)/1000
                min_q_mvar = (haskey(inverter_row, :min_q_mvar) ? safe_get_value(inverter_row[:min_q_mvar], 0.0, Float64) : 0.0)/1000
                
                # Validate power limits are reasonable
                if min_p_mw > max_p_mw && max_p_mw != 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid active power limits (min: $min_p_mw MW, max: $max_p_mw MW), swapping values"
                    min_p_mw, max_p_mw = max_p_mw, min_p_mw
                end
                
                if min_q_mvar > max_q_mvar && max_q_mvar != 0.0
                    @warn "Row $i: Converter $name (ID: $index) has invalid reactive power limits (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar), swapping values"
                    min_q_mvar, max_q_mvar = max_q_mvar, min_q_mvar
                end
                
                # Extract PV array parameters
                numpanelseries = safe_get_value(row[:numpanelseries], 0, Int)
                numpanelparallel = safe_get_value(row[:numpanelparallel], 0, Int)
                vmpp = safe_get_value(row[:vmpp], 0.0, Float64)
                impp = safe_get_value(row[:impp], 0.0, Float64)
                voc = safe_get_value(row[:voc], 0.0, Float64)
                isc = safe_get_value(row[:isc], 0.0, Float64)
                pvanumcells = safe_get_value(row[:pvanumcells], 0, Int)
                
                # Validate PV parameters are valid
                if numpanelseries <= 0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid number of series panels ($numpanelseries), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if numpanelparallel <= 0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid number of parallel panels ($numpanelparallel), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if vmpp <= 0.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid maximum power point voltage ($vmpp V), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if impp <= 0.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid maximum power point current ($impp A), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if voc <= 0.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid open circuit voltage ($voc V), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if isc <= 0.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid short circuit current ($isc A), skipping this row"
                    error_rows += 1
                    continue
                end
                
                if pvanumcells <= 0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid number of cells ($pvanumcells), skipping this row"
                    error_rows += 1
                    continue
                end
                
                # Environmental parameters
                irradiance = haskey(row, :irradiance) ? safe_get_value(row[:irradiance], 1000.0, Float64) : 1000.0
                # Validate irradiance value is reasonable
                if irradiance < 0.0 || irradiance > 2000.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid irradiance value ($irradiance W/m²), setting to default value 1000.0 W/m²"
                    irradiance = 1000.0
                end
                
                temperature = haskey(row, :temperature) ? safe_get_value(row[:temperature], 25.0, Float64) : 25.0
                # Validate temperature value is reasonable
                if temperature < -40.0 || temperature > 85.0
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid temperature value ($temperature °C), setting to default value 25.0 °C"
                    temperature = 25.0
                end
                
                # Temperature coefficient parameters
                α_isc = haskey(row, :aisctemp) ? safe_get_value(row[:aisctemp], 0.0, Float64) : 0.0
                β_voc = haskey(row, :bvoctemp) ? safe_get_value(row[:bvoctemp], 0.0, Float64) : 0.0
                
                # Control parameters
                control_mode = haskey(inverter_row, :acoperationmode) ? safe_get_value(inverter_row[:acoperationmode], "", String) : "Voltage Control"
                # Validate control mode
                valid_control_modes = ["Voltage Control", "Swing", "Mvar Control", "PF Control"]
                if !(control_mode in valid_control_modes)
                    @warn "Row $i: AC PV system $name (ID: $index) has invalid control mode ($control_mode), setting to default value 'Voltage Control'"
                    control_mode = "Voltage Control"
                end
                
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # Operational status
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # Create ACPVSystem object and add to case
                push!(case.ACPVSystems, ACPVSystem(
                    index, name, bus_id, p_mw, q_mvar, vm_ac_pu, vm_dc_pu,loss_percent, loss_mw,
                    max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                    numpanelseries, numpanelparallel, vmpp, impp, voc, isc, pvanumcells,
                    irradiance, temperature, α_isc, β_voc,
                    control_mode, controllable, in_service
                ))
                
                processed_rows += 1
                
            catch e
                @error "Error processing AC PV system data row $i" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "AC PV system data loading complete: Successfully processed $processed_rows rows, $error_rows errors"
        
    catch e
        @error "Error loading AC PV system data" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_vpps!(case::JuliaPowerCase, file_path::String, sheet_name::String)

Load virtual power plant (VPP) data from Excel file and add to power system case.

Parameters:
- `case::JuliaPowerCase`: Power system case
- `file_path::String`: Excel file path
- `sheet_name::String`: Worksheet name containing virtual power plant data
"""
function load_vpps!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    # Use DataFrame processing
    df = DataFrame(XLSX.readtable(file_path, sheet_name))
    
    # Ensure data is not empty
    if isempty(df)
        @info "Virtual power plant table is empty"
        return
    end
    
    # Convert column names to lowercase
    rename!(df, lowercase.(names(df)))
    
    # Iterate through each row of data
    for row in eachrow(df)
        try
            # Extract basic field values from row data
            index = safe_get_value(row[:index], 0, Int)
            name = safe_get_value(row[:name], "", String)
            description = haskey(row, :description) ? safe_get_value(row[:description], "", String) : ""
            control_area = haskey(row, :control_area) ? safe_get_value(row[:control_area], "", String) : ""
            
            # Capacity and energy parameters
            capacity_mw = haskey(row, :capacity_mw) ? safe_get_value(row[:capacity_mw], 0.0, Float64) : 0.0
            energy_mwh = haskey(row, :energy_mwh) ? safe_get_value(row[:energy_mwh], 0.0, Float64) : 0.0
            
            # Response and ramp parameters
            response_time_s = haskey(row, :response_time_s) ? safe_get_value(row[:response_time_s], 0.0, Float64) : 0.0
            ramp_rate_mw_per_min = haskey(row, :ramp_rate_mw_per_min) ? safe_get_value(row[:ramp_rate_mw_per_min], 0.0, Float64) : 0.0
            availability_percent = haskey(row, :availability_percent) ? safe_get_value(row[:availability_percent], 100.0, Float64) : 100.0
            
            # Operational information
            operator = haskey(row, :operator) ? safe_get_value(row[:operator], "", String) : ""
            in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
            
            # Collect additional kwargs parameters
            kwargs = Dict{Symbol, Any}()
            
            # Resource information
            if haskey(row, :resource_type)
                kwargs[:resource_type] = safe_get_value(row[:resource_type], "", String)
            end
            
            if haskey(row, :resource_id)
                kwargs[:resource_id] = safe_get_value(row[:resource_id], 0, Int)
            end
            
            if haskey(row, :capacity_share_percent)
                kwargs[:capacity_share_percent] = safe_get_value(row[:capacity_share_percent], 0.0, Float64)
            end
            
            if haskey(row, :control_priority)
                kwargs[:control_priority] = safe_get_value(row[:control_priority], 0, Int)
            end
            
            if haskey(row, :resource_response_time_s)
                kwargs[:resource_response_time_s] = safe_get_value(row[:resource_response_time_s], 0.0, Float64)
            end
            
            if haskey(row, :max_duration_h)
                kwargs[:max_duration_h] = safe_get_value(row[:max_duration_h], 0.0, Float64)
            end
            
            # Load information
            if haskey(row, :timestamp) && !ismissing(row[:timestamp])
                # Process timestamp, adjust according to actual format
                if isa(row[:timestamp], String)
                    kwargs[:timestamp] = DateTime(row[:timestamp], dateformat"yyyy-mm-dd HH:MM:SS")
                elseif isa(row[:timestamp], DateTime)
                    kwargs[:timestamp] = row[:timestamp]
                else
                    kwargs[:timestamp] = DateTime(now())
                end
            else
                kwargs[:timestamp] = DateTime(now())
            end
            
            if haskey(row, :p_mw)
                kwargs[:p_mw] = safe_get_value(row[:p_mw], 0.0, Float64)
            end
            
            if haskey(row, :q_mvar)
                kwargs[:q_mvar] = safe_get_value(row[:q_mvar], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_up_mw)
                kwargs[:flexibility_up_mw] = safe_get_value(row[:flexibility_up_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_down_mw)
                kwargs[:flexibility_down_mw] = safe_get_value(row[:flexibility_down_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_duration_h)
                kwargs[:flexibility_duration_h] = safe_get_value(row[:flexibility_duration_h], 0.0, Float64)
            end
            
            # Create VirtualPowerPlant object and add to case
            vpp = VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                                   response_time_s, ramp_rate_mw_per_min, availability_percent,
                                   operator, in_service; kwargs...)
            
            push!(case.vpps, vpp)
            
            # Update index mapping
            case.vpp_indices[index] = length(case.vpps)
            
        catch e
            @warn "Error processing virtual power plant: $e"
            @warn "Problem row: $(row)"
        end
    end
    
    @info "Loaded $(length(case.vpps)) virtual power plants"
end
