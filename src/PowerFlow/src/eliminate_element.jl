"""
Eliminate a specific element from the power system data structure.

Parameters:
- jpc: The power system data structure.
- element_type: The type of element to eliminate (e.g., "bus", "gen", "branch").
- element_id: The identifier of the element to eliminate.

Returns:
- jpc: The updated power system data structure with the specified element removed.
"""
function eliminate_element(jpc)

    # check if the jpc has the element to eliminate
    pv_acsystem = jpc["pv_acsystem"]  
    if !isempty(pv_acsystem)
        bus_ids = pv_acsystem[:, PV_AC_BUS]  # Assuming the first column contains bus IDs
        # Find the indices of the elements to eliminate
        indices_to_eliminate = findall(jpc["genAC"][:, GEN_BUS].== bus_ids)  
        # Eliminate the elements from genAC
        if !isempty(indices_to_eliminate)
            keep_indices = setdiff(1:size(jpc["genAC"], 1), indices_to_eliminate)
            jpc["genAC"] = jpc["genAC"][keep_indices, :]
        end
    end
    
    return jpc
end