function parse_matlab_case_file(filepath)
    # Read file content
    content = read(filepath, String)
    
    # Create empty dictionary to store intermediate results
    mpc = Dict{String, Any}()
    
    # Parse baseMVA
    if occursin(r"mpc\.baseMVA\s*=\s*(\d+)", content)
        basemva_match = match(r"mpc\.baseMVA\s*=\s*(\d+)", content)
        mpc["baseMVA"] = parse(Float64, basemva_match[1])
    end
    
    # Parse version
    if occursin(r"mpc\.version\s*=\s*'(\d+)'", content)
        version_match = match(r"mpc\.version\s*=\s*'(\d+)'", content)
        mpc["version"] = version_match[1]
    end
    
    # Function to parse matrix or string data
    function extract_data(content, key)
        # Split content into lines
        lines = split(content, '\n')
        
        # Find the starting line of the matrix
        start_pattern = "mpc.$key = ["
        end_pattern = "];"
        
        start_idx = 0
        end_idx = 0
        
        # Find the start and end positions of the matrix
        for (i, line) in enumerate(lines)
            if occursin(start_pattern, line)
                start_idx = i
            elseif start_idx > 0 && occursin(end_pattern, line)
                end_idx = i
                break
            end
        end
        
        # If matrix is found
        if start_idx > 0 && end_idx > 0
            # Extract matrix content
            matrix_lines = lines[start_idx+1:end_idx-1]
            
            # Filter out empty lines and comment lines
            matrix_lines = filter(line -> !isempty(strip(line)) && !startswith(strip(line), '%'), matrix_lines)
            
            # Check if it contains string data
            contains_strings = any(line -> occursin("'", line) || occursin("\"", line), matrix_lines)
            
            if contains_strings
                # Process string data
                matrix = []
                for line in matrix_lines
                    # Remove semicolon and comments at the end of line
                    line = replace(line, r";.*$" => "")
                    # Extract content within quotes and numbers
                    parts = String[]
                    current_str = ""
                    in_quotes = false
                    quote_char = nothing
                    
                    for char in line
                        if char in ['\'', '"']
                            if !in_quotes
                                in_quotes = true
                                quote_char = char
                            elseif char == quote_char
                                in_quotes = false
                                if !isempty(current_str)
                                    push!(parts, current_str)
                                    current_str = ""
                                end
                            else
                                current_str *= char
                            end
                        elseif in_quotes
                            current_str *= char
                        elseif !isspace(char) && char != ';'
                            current_str *= char
                        elseif !isempty(current_str)
                            push!(parts, current_str)
                            current_str = ""
                        end
                    end
                    
                    if !isempty(current_str)
                        push!(parts, current_str)
                    end
                    
                    # Filter out empty strings
                    parts = filter(!isempty, parts)
                    
                    if !isempty(parts)
                        push!(matrix, parts)
                    end
                end
                return length(matrix) > 0 ? reduce(vcat, transpose.(matrix)) : nothing
            else
                # Process numerical data
                matrix = []
                for line in matrix_lines
                    # Remove semicolon and comments at the end of line
                    line = replace(line, r";.*$" => "")
                    # Split and convert to numerical values
                    try
                        row = parse.(Float64, split(strip(line)))
                        if !isempty(row)
                            push!(matrix, row)
                        end
                    catch
                        @warn "Unable to parse line: $line"
                        continue
                    end
                end
                return length(matrix) > 0 ? reduce(vcat, transpose.(matrix)) : nothing
            end
        end
        return nothing
    end
    
    # Find all possible matrix names
    matrix_names = String[]
    for line in split(content, '\n')
        m = match(r"mpc\.(\w+)\s*=\s*\[", line)
        if m !== nothing
            push!(matrix_names, m[1])
        end
    end
    
    # Parse each found matrix
    for name in matrix_names
        if name ∉ ["version", "baseMVA"]  # Skip already processed special fields
            matrix = extract_data(content, name)
            if matrix !== nothing
                mpc[name] = matrix
            end
        end
    end
    
    return mpc
end

# Convert dictionary format to JPC structure
function dict_to_jpc(mpc_dict)
    # Create JPC instance
    jpc = JPC()
    
    # Set basic attributes
    if haskey(mpc_dict, "version")
        jpc.version = mpc_dict["version"]
    end
    
    if haskey(mpc_dict, "baseMVA")
        jpc.baseMVA = mpc_dict["baseMVA"]
    end
    
    # Map MATPOWER fields to JPC fields
    field_mapping = Dict(
        "bus" => "busAC",
        "gen" => "genAC",
        "branch" => "branchAC",
        "load" => "loadAC"
        # Add more mappings as needed
    )
    
    # Transfer data from dictionary to JPC structure
    for (matpower_field, jpc_field) in field_mapping
        if haskey(mpc_dict, matpower_field)
            data = mpc_dict[matpower_field]
            # Ensure data is a 2D array
            if ndims(data) == 1
                data = reshape(data, 1, length(data))
            end
            
            # Get reference to corresponding field in JPC structure
            field_ref = getfield(jpc, Symbol(jpc_field))
            
            # Check if data columns match target field
            target_cols = size(field_ref, 2)
            data_cols = size(data, 2)
            
            if data_cols <= target_cols
                # Create new array of appropriate size
                new_data = zeros(size(data, 1), target_cols)
                # Copy data
                new_data[:, 1:data_cols] = data
                # Set field value
                setfield!(jpc, Symbol(jpc_field), new_data)
            else
                # If data columns exceed target field, truncate data
                setfield!(jpc, Symbol(jpc_field), data[:, 1:target_cols])
                @warn "Data $matpower_field has $data_cols columns exceeding target field $jpc_field with $target_cols columns, truncated."
            end
        end
    end
    
    # Handle other unmapped fields
    for field in keys(mpc_dict)
        if field ∉ ["version", "baseMVA"] && !haskey(field_mapping, field)
            # Try to directly match field names
            try
                if hasproperty(jpc, Symbol(field))
                    data = mpc_dict[field]
                    # Ensure data is a 2D array
                    if ndims(data) == 1
                        data = reshape(data, 1, length(data))
                    end
                    
                    # Get reference to corresponding field in JPC structure
                    field_ref = getfield(jpc, Symbol(field))
                    
                    # Check if data columns match target field
                    target_cols = size(field_ref, 2)
                    data_cols = size(data, 2)
                    
                    if data_cols <= target_cols
                        # Create new array of appropriate size
                        new_data = zeros(size(data, 1), target_cols)
                        # Copy data
                        new_data[:, 1:data_cols] = data
                        # Set field value
                        setfield!(jpc, Symbol(field), new_data)
                    else
                        # If data columns exceed target field, truncate data
                        setfield!(jpc, Symbol(field), data[:, 1:target_cols])
                        @warn "Data $field has $data_cols columns exceeding target field with $target_cols columns, truncated."
                    end
                end
            catch e
                @warn "Unable to map field $field to JPC structure: $e"
            end
        end
    end
    
    return jpc
end

function save_to_julia_file(jpc, output_filepath)
    open(output_filepath, "w") do f
        write(f, """function case_data()
    # Create JPC structure instance
    jpc = JPC("$(jpc.version)", $(jpc.baseMVA), $(jpc.success), $(jpc.iterationsAC), $(jpc.iterationsDC))
    
""")
        
        # Get all field names of JPC structure
        field_names = fieldnames(JPC)
        
        # Iterate through all fields
        for field in field_names
            # Skip basic attributes
            if field ∉ [:version, :baseMVA, :success, :iterationsAC, :iterationsDC]
                data = getfield(jpc, field)
                
                # Only process non-empty arrays
                if !isempty(data)
                    write(f, "    # Set $field data\n")
                    write(f, "    jpc.$field = [\n")
                    
                    # Write array data
                    for i in 1:size(data, 1)
                        write(f, "        ")
                        for j in 1:size(data, 2)
                            write(f, "$(data[i,j]) ")
                        end
                        write(f, ";\n")
                    end
                    
                    write(f, "    ]\n\n")
                end
            end
        end
        
        write(f, "    return jpc\nend")
    end
end

function convert_matpower_case(input_filepath, output_filepath)
    try
        println("Parsing MATLAB file...")
        mpc_dict = parse_matlab_case_file(input_filepath)
        
        println("Converting to JPC structure...")
        jpc = dict_to_jpc(mpc_dict)
        
        println("Saving to Julia file...")
        save_to_julia_file(jpc, output_filepath)
        
        println("Conversion completed!")
        println("Input file: $input_filepath")
        println("Output file: $output_filepath")
        
        # Print data statistics
        println("\nData statistics:")
        for field in fieldnames(JPC)
            data = getfield(jpc, field)
            if isa(data, Array)
                if !isempty(data)
                    println("$field matrix size: $(size(data))")
                end
            else
                println("$field: $(data)")
            end
        end
        
        return jpc
    catch e
        println("Error occurred during conversion:")
        println(e)
        return nothing
    end
end
