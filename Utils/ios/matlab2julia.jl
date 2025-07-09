"""
    parse_matlab_case_file(filepath)

Parse a MATPOWER case file in MATLAB format and convert it to a Julia dictionary.
The function extracts baseMVA, version, and all matrix data from the file.

# Arguments
- `filepath::String`: Path to the MATLAB case file

# Returns
- `Dict{String, Any}`: Dictionary containing all parsed data from the case file
"""
function parse_matlab_case_file(filepath)
    # Read file content
    content = read(filepath, String)
    
    # Create empty dictionary to store results
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
    
    """
        extract_data(content, key)
    
    Internal helper function to extract matrix or string data from MATLAB content.
    
    # Arguments
    - `content::String`: Full content of the MATLAB file
    - `key::String`: Name of the matrix to extract
    
    # Returns
    - Matrix data or nothing if extraction fails
    """
    function extract_data(content, key)
        # Split content into lines
        lines = split(content, '\n')
        
        # Find the line where matrix starts
        start_pattern = "mpc.$key = ["
        end_pattern = "];"
        
        start_idx = 0
        end_idx = 0
        
        # Find start and end positions of the matrix
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
            
            # Check if contains string data
            contains_strings = any(line -> occursin("'", line) || occursin("\"", line), matrix_lines)
            
            if contains_strings
                # Process string data
                matrix = []
                for line in matrix_lines
                    # Remove semicolons and comments at end of line
                    line = replace(line, r";.*$" => "")
                    # Extract content in quotes and numbers
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
                # Process numeric data
                matrix = []
                for line in matrix_lines
                    # Remove semicolons and comments at end of line
                    line = replace(line, r";.*$" => "")
                    # Split and convert to numeric values
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

"""
    save_to_julia_file(mpc, output_filepath)

Save a parsed MATPOWER case as a Julia file with a function that returns the data.

# Arguments
- `mpc::Dict{String, Any}`: Dictionary containing the parsed MATPOWER case data
- `output_filepath::String`: Path where the Julia file should be saved

# Returns
- Nothing, but creates a Julia file at the specified location
"""
function save_to_julia_file(mpc, output_filepath)
    open(output_filepath, "w") do f
        write(f, """function case_data()
    mpc = Dict{String, Any}()

""")
        
        # Write version
        if haskey(mpc, "version")
            write(f, "    mpc[\"version\"] = \"$(mpc["version"])\"\n\n")
        end
        
        # Write baseMVA
        if haskey(mpc, "baseMVA")
            write(f, "    mpc[\"baseMVA\"] = $(mpc["baseMVA"])\n\n")
        end
        
        # Write all matrix data
        for key in keys(mpc)
            if key ∉ ["version", "baseMVA"]
                matrix = mpc[key]
                write(f, "    mpc[\"$key\"] = [\n")
                
                # Check matrix dimensions
                if ndims(matrix) == 1
                    # Handle one-dimensional arrays
                    write(f, "        ")
                    for value in matrix
                        if typeof(value) <: AbstractString
                            write(f, "\"$value\" ")  # String values are surrounded by quotes
                        else
                            write(f, "$(value) ")    # Numeric values are written directly
                        end
                    end
                    write(f, "\n")
                else
                    # Handle two-dimensional arrays
                    for i in eachindex(matrix[:,1])
                        write(f, "        ")
                        for j in eachindex(matrix[1,:])
                            value = matrix[i,j]
                            if typeof(value) <: AbstractString
                                write(f, "\"$value\" ")  # String values are surrounded by quotes
                            else
                                write(f, "$(value) ")    # Numeric values are written directly
                            end
                        end
                        write(f, ";\n")
                    end
                end
                
                write(f, "    ]\n\n")
            end
        end
        
        write(f, "    return mpc\nend")
    end
end

"""
    convert_matpower_case(input_filepath, output_filepath)

Convert a MATPOWER case file from MATLAB format to Julia format.
This function handles the entire conversion process, including parsing the MATLAB file,
converting the data structures, and saving the result as a Julia file.

# Arguments
- `input_filepath::String`: Path to the input MATLAB case file
- `output_filepath::String`: Path where the output Julia file should be saved

# Returns
- `Dict{String, Any}` or `nothing`: The parsed case data if successful, nothing if an error occurs
"""
function convert_matpower_case(input_filepath, output_filepath)
    try
        println("Parsing MATLAB file...")
        mpc = parse_matlab_case_file(input_filepath)
        
        println("Saving as Julia file...")
        save_to_julia_file(mpc, output_filepath)
        
        println("Conversion completed!")
        println("Input file: $input_filepath")
        println("Output file: $output_filepath")
        
        # Print data statistics
        println("\nData statistics:")
        for (key, value) in mpc
            if value isa Array
                println("$key matrix size: $(size(value))")
            else
                println("$key: $(value)")
            end
        end
        
        return mpc
    catch e
        println("Error occurred during conversion:")
        println(e)
        return nothing
    end
end
