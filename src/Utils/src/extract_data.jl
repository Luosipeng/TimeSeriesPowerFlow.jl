"""
    extract_data(sheet_name::String, sheets_data::Dict{String, DataFrame})::DataFrame

Extract data from an Excel worksheet.
    
This function processes a worksheet to find and extract the valid data region,
determining the effective number of rows and columns that contain data.

# Arguments
- `sheet_name::String`: Name of the worksheet to extract data from
- `sheets_data::Dict{String, DataFrame}`: Dictionary containing all worksheet data

# Returns
- `DataFrame`: Extracted data subset, excluding headers
"""
function extract_data(sheet_name::String, sheets_data::Dict{String, DataFrame})::DataFrame

    # Check if worksheet exists
    if !haskey(sheets_data, sheet_name)
        # @warn "Worksheet '$sheet_name' does not exist"
        return DataFrame()
    end
    
    # Get current worksheet
    current_sheet::DataFrame = sheets_data[sheet_name]
    
    # Initialize counters
    row_count::Int = 0
    col_count::Int = 0
    
    # Calculate effective number of rows
    sheet_rows::Int = size(current_sheet, 1)
    for row_idx in 1:sheet_rows
        if ismissing(current_sheet[row_idx, 1])
            row_count = row_idx - 1
            break
        end
        row_count = row_idx
    end
    
    # Calculate effective number of columns
    sheet_cols::Int = size(current_sheet, 2)
    for col_idx in 1:sheet_cols
        if ismissing(current_sheet[1, col_idx])
            col_count = col_idx - 1
            break
        end
        col_count = col_idx
    end
    
    # Check if valid data exists
    if row_count == 0 || col_count == 0
        @warn "No valid data found in worksheet '$sheet_name'"
        return DataFrame()
    end
    
    # Return the extracted data (starting from second row, skipping headers)
    return current_sheet[2:row_count, 1:col_count]
end
