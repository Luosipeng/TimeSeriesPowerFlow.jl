"""
    merge_results(results)

Merge power flow calculation results from multiple isolated islands.

# Arguments
- `results`: An array of JPC objects containing power flow results for different islands

# Returns
- `merged_result`: A single JPC object containing the combined results
- `area`: The number of islands that were merged

# Description
This function combines power flow results from multiple isolated islands into a single
comprehensive result. It merges all data matrices (buses, branches, generators, etc.)
from the individual island results and sorts them by their ID numbers.

The function:
1. Creates a new JPC object to hold the merged results
2. Combines basic result fields like success status and iteration counts
3. Merges all data matrices from the input results
4. Sorts the merged matrices by their first column (typically ID numbers)

# Notes
- Assumes all input results use the same base MVA
- Success is determined by the logical AND of all individual results' success flags
- Iteration count is the maximum of all individual results' iteration counts
"""
function merge_results(results)
    # Create a new JPC object
    merged_result = JPC()
    
    # Merge basic result fields
    merged_result.success = all([r.success for r in results])
    merged_result.iterationsAC = maximum([r.iterationsAC for r in results])
    
    # Set baseMVA (assuming all results use the same base power)
    if !isempty(results)
        merged_result.baseMVA = results[1].baseMVA
        merged_result.version = results[1].version
    end
    
    # Merge main data matrices and sort by first column
    for key in [:busAC, :branchAC, :genAC, :loadAC, :loadAC_flex, :loadAC_asymm, 
                :branch3ph, :busDC, :branchDC, :sgenAC, :storage, :sgenDC, 
                :converter, :ext_grid, :hvcb, :microgrid]
        # Collect corresponding matrices from all results
        matrices = []
        for r in results
            if hasproperty(r, key) && !isempty(getproperty(r, key))
                push!(matrices, getproperty(r, key))
            end
        end
        
        if !isempty(matrices)
            # Vertically concatenate all matrices
            combined = vcat(matrices...)
            
            # Sort by first column (if matrix is not empty)
            if !isempty(combined)
                sorted_indices = sortperm(combined[:, 1])
                setproperty!(merged_result, key, combined[sorted_indices, :])
            else
                setproperty!(merged_result, key, combined)
            end
        end
    end
    
    area = length(results)
    return merged_result, area
end
