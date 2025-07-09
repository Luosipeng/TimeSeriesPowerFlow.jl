"""
    dc_preprocess(mpc, opt)

Perform DC power flow preprocessing on multiple islands of a power system network.

# Arguments
- `mpc`: Power system data in MATPOWER case format
- `opt`: Options dictionary containing configuration parameters

# Returns
- `mpc_list`: List of processed power system islands
- `isolated`: Information about isolated components in the network

# Description
This function handles preprocessing for DC power flow analysis by:
1. Extracting islands (electrically separate parts) from the power system
2. Optionally running DC power flow preprocessing on each island in parallel
3. Updating voltage angles in the original data structure

The parallel processing is implemented using Julia's multi-threading capabilities
for both the preprocessing step and the angle update step, which can significantly
improve performance for large systems with multiple islands.
"""
function dc_preprocess(mpc, opt)
    mpc_list, isolated = PowerFlow.extract_islands(mpc)
    n_islands = length(mpc_list)
    if(opt["PF"]["DC_PREPROCESS"]==1)   
        preconditioned_list = Vector{Any}(undef, n_islands)
        
        @threads for i in 1:n_islands
            preconditioned_list[i] = PowerFlow.runprepf(mpc_list[i], opt)
        end
        
        # Update voltage angles using multi-threading
        @threads for i in 1:n_islands
            mpc_list[i]["bus"][:, 9] = preconditioned_list[i]["bus"][:, 9]
        end
    end
    return mpc_list, isolated
end
