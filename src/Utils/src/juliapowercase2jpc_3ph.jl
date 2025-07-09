"""
    JuliaPowerCase2Jpc_3ph(case::JuliaPowerCase)

Convert a JuliaPowerCase object to a three-phase JPC_3ph object.

This function processes a power system case in JuliaPowerCase format and converts it into
a three-phase sequence component model (JPC_3ph) for three-phase power flow analysis.

# Arguments
- `case::JuliaPowerCase`: The power system case data in JuliaPowerCase format

# Returns
- `jpc_3ph::JPC_3ph`: The converted three-phase power system model
- `gs_eg::Vector{Float64}`: Conductance values for external grids
- `bs_eg::Vector{Float64}`: Susceptance values for external grids

# Process
1. Merges virtual nodes
2. Creates a new JPC_3ph object
3. Sets base parameters
4. Processes bus data for all sequence components
5. Processes branch data for all sequence components
6. Processes generator data for all sequence components
7. Processes load data for all sequence components
8. Adds external grid short circuit impedance data
"""
function JuliaPowerCase2Jpc_3ph(case::JuliaPowerCase)
    # 1. Merge virtual nodes
    case = PowerFlow.merge_virtual_nodes(case)

    # 2. Create JPC_3ph object
    jpc_3ph = PowerFlow.JPC_3ph()

    # 3. Set basic parameters
    jpc_3ph.baseMVA = case.baseMVA
    jpc_3ph.basef = case.basef

    # 4. Set bus data
    JPC_3ph_buses_process(case, jpc_3ph)

    # 5. Set branch data
    JPC_3ph_branches_process(case, jpc_3ph)
    
    # 6. Set generator data
    JPC_3ph_gens_process(case, jpc_3ph)
    
    # 7. Set load data
    JPC_3ph_loads_process(case, jpc_3ph)

    # 8. Add external grid short circuit impedance data
    jpc_3ph, _, _ = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 1)
    jpc_3ph, gs_eg, bs_eg = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 2)
    jpc_3ph, _ , _ = JPC_3ph_add_grid_external_sc_impedance(case, jpc_3ph, 0)

    return jpc_3ph, gs_eg, bs_eg
end

"""
    JPC_3ph_buses_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Process bus data from a JuliaPowerCase and populate the JPC_3ph object with sequence component bus data.

This function creates bus matrices for positive, negative, and zero sequence components.
For positive sequence, voltage magnitude is initialized to 1.0 p.u., while for negative and 
zero sequence components, voltage magnitude is initialized to 0.0 p.u.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Returns
- `jpc_3ph::Utils.JPC_3ph`: The updated JPC_3ph object with bus data
"""
function JPC_3ph_buses_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Get bus data and deep copy to prevent unintended modifications
    buses = deepcopy(case.busesAC)
    
    # Create an empty matrix with rows equal to number of buses and 13 columns
    num_buses = length(buses)
    bus_matrix = zeros(num_buses, 13)
    
    for (i, bus) in enumerate(buses)
        # Set initial voltage values (based on sequence)
        vm = 0.0
        va = 0.0
        
        # Fill each row of the matrix
        bus_matrix[i, :] = [
            bus.bus_id,      # Bus ID
            1.0,             # Bus type (all set to PQ buses)
            0.0,             # PD (MW) Active load
            0.0,             # QD (MVAR) Reactive load
            0.0,             # GS (MW) Active generation
            0.0,             # BS (MVAR) Reactive generation
            bus.area_id,     # Area number
            vm,              # Bus voltage magnitude (p.u.)
            va,              # Bus voltage angle (degrees)
            bus.vn_kv,       # Base voltage (kV)
            bus.zone_id,     # Zone number
            bus.max_vm_pu,   # Maximum voltage magnitude (p.u.)
            bus.min_vm_pu,   # Minimum voltage magnitude (p.u.)
        ]
    end
    
    # Store results for all sequences in busAC fields
    jpc_3ph.busAC_0 = bus_matrix
    jpc_3ph.busAC_2 = bus_matrix

    # Create positive sequence bus matrix with voltage magnitude initialized to 1.0
    bus_matrix_1 = zeros(num_buses, 13)
    for (i, bus) in enumerate(buses)
        # Set initial voltage values for positive sequence
        vm = 1.0
        va = 0.0
        
        # Fill each row of the matrix
        bus_matrix_1[i, :] = [
            bus.bus_id,      # Bus ID
            1.0,             # Bus type (all set to PQ buses)
            0.0,             # PD (MW) Active load
            0.0,             # QD (MVAR) Reactive load
            0.0,             # GS (MW) Active generation
            0.0,             # BS (MVAR) Reactive generation
            bus.area_id,     # Area number
            vm,              # Bus voltage magnitude (p.u.)
            va,              # Bus voltage angle (degrees)
            bus.vn_kv,       # Base voltage (kV)
            bus.zone_id,     # Zone number
            bus.max_vm_pu,   # Maximum voltage magnitude (p.u.)
            bus.min_vm_pu,   # Minimum voltage magnitude (p.u.)
        ]
    end

    # Store positive sequence results
    jpc_3ph.busAC_1 = bus_matrix_1
    
    return jpc_3ph
end

"""
    JPC_3ph_branches_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Process branch data from a JuliaPowerCase and populate the JPC_3ph object with sequence component branch data.

This function calculates parameters for lines and transformers for all sequence components
(positive, negative, and zero).

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Returns
- `jpc_3ph::Utils.JPC_3ph`: The updated JPC_3ph object with branch data
"""
function JPC_3ph_branches_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Calculate line parameters
    calculate_3ph_line_parameters(case, jpc_3ph)
    # Calculate transformer parameters
    calculate_3ph_transformer2w_parameters(case, jpc_3ph)
    # Calculate zero sequence branch parameters
    jpc_3ph = calculate_branch_JPC_zero(case, jpc_3ph)

    return jpc_3ph
end

"""
    calculate_3ph_line_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Calculate line parameters for positive and negative sequence components.

This function processes line data from a JuliaPowerCase and converts it to JPC format
for both positive and negative sequence components. For most lines, negative sequence
parameters are equal to positive sequence parameters.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Side effects
- Updates `jpc_3ph.branchAC_1` (positive sequence) and `jpc_3ph.branchAC_2` (negative sequence)
"""
function calculate_3ph_line_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Process line data, convert to JPC format
    nbr = length(case.branchesAC)
    branch_1 = zeros(nbr, 14)  # Positive sequence
    branch_2 = zeros(nbr, 14)  # Negative sequence
    lines = case.branchesAC

    for (i, line) in enumerate(lines)
        # Get from and to bus numbers
        from_bus_idx = line.from_bus
        to_bus_idx = line.to_bus
        
        # Get base voltage (kV) of the from bus
        basekv = jpc_3ph.busAC_1[from_bus_idx, BASE_KV]  # Use three-phase system bus data
        
        # Calculate base impedance
        baseR = (basekv^2) / case.baseMVA
        
        # Consider parallel lines
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # Calculate positive sequence impedance parameters
        r1_pu = line.length_km * line.r_ohm_per_km / baseR / parallel
        x1_pu = line.length_km * line.x_ohm_per_km / baseR / parallel
        
        # Calculate negative sequence impedance parameters (usually same as positive sequence)
        r2_pu = r1_pu  # For most lines, negative sequence resistance equals positive sequence resistance
        x2_pu = x1_pu  # For most lines, negative sequence reactance equals positive sequence reactance
        
        # If line has specific negative sequence parameters, use them:
        if hasfield(typeof(line), :r2_ohm_per_km)
            r2_pu = line.length_km * line.r2_ohm_per_km / baseR / parallel
        end
        if hasfield(typeof(line), :x2_ohm_per_km)
            x2_pu = line.length_km * line.x2_ohm_per_km / baseR / parallel
        end
        
        # Calculate shunt susceptance (p.u.) - usually the same for positive and negative sequence
        b_pu = 2 * π * case.basef * line.length_km * line.c_nf_per_km * 1e-9 * baseR * parallel
        
        # Calculate shunt conductance (p.u.)
        g_pu = 0.0
        if hasfield(typeof(line), :g_us_per_km)
            g_pu = line.g_us_per_km * 1e-6 * baseR * line.length_km * parallel
        end
        
        # Fill positive sequence matrix (branchAC_1)
        branch_1[i, F_BUS] = from_bus_idx
        branch_1[i, T_BUS] = to_bus_idx
        branch_1[i, BR_R] = r1_pu
        branch_1[i, BR_X] = x1_pu
        branch_1[i, BR_B] = b_pu
        
        # Fill negative sequence matrix (branchAC_2)
        branch_2[i, F_BUS] = from_bus_idx
        branch_2[i, T_BUS] = to_bus_idx
        branch_2[i, BR_R] = r2_pu
        branch_2[i, BR_X] = x2_pu
        branch_2[i, BR_B] = b_pu
        
        # Set rated capacity (same for positive and negative sequence)
        rate_a = 100.0  # Default value
        if hasfield(typeof(line), :max_i_ka)
            rate_a = line.max_i_ka * basekv * sqrt(3)  # Rated capacity (MVA)
        end
        branch_1[i, RATE_A] = rate_a
        branch_2[i, RATE_A] = rate_a
        
        # Set branch status (same for positive and negative sequence)
        status = line.in_service ? 1.0 : 0.0
        branch_1[i, BR_STATUS] = status
        branch_2[i, BR_STATUS] = status
        
        # Set angle limits (same for positive and negative sequence)
        branch_1[i, ANGMIN] = -360.0
        branch_1[i, ANGMAX] = 360.0
        branch_2[i, ANGMIN] = -360.0
        branch_2[i, ANGMAX] = 360.0
    end

    # Assign to three-phase system
    jpc_3ph.branchAC_1 = branch_1  # Positive sequence
    jpc_3ph.branchAC_2 = branch_2  # Negative sequence
end

"""
    calculate_3ph_transformer2w_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Calculate two-winding transformer parameters for positive and negative sequence components.

This function processes two-winding transformer data from a JuliaPowerCase and converts
it to JPC format for both positive and negative sequence components. For most transformers,
negative sequence parameters are equal to positive sequence parameters.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Side effects
- Updates `jpc_3ph.branchAC_1` (positive sequence) and `jpc_3ph.branchAC_2` (negative sequence)
"""
function calculate_3ph_transformer2w_parameters(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Process transformer data, convert to JPC format
    transformers = case.transformers_2w_etap
    nbr = length(transformers)
    
    if nbr == 0
        return  # If no transformers, return directly
    end
    
    # Create transformer branch matrices - positive and negative sequence
    branch_1 = zeros(nbr, 14)  # Positive sequence
    branch_2 = zeros(nbr, 14)  # Negative sequence
    
    for (i, transformer) in enumerate(transformers)
        # Get high voltage and low voltage bus numbers
        hv_bus_idx = transformer.hv_bus
        lv_bus_idx = transformer.lv_bus
        
        # Get high voltage side base voltage (kV)
        hv_basekv = jpc_3ph.busAC_1[hv_bus_idx, BASE_KV]  # Use three-phase system bus data
        
        # Calculate impedance parameters
        # Convert transformer impedance percentage to per unit
        z_percent = transformer.z_percent
        x_r_ratio = transformer.x_r
        
        # Calculate resistance and reactance (considering base power conversion)
        s_ratio = transformer.sn_mva / case.baseMVA
        z_pu = z_percent / 100.0 / s_ratio  # Convert to system base (percentage to decimal)
        
        # Calculate positive sequence impedance parameters
        r1_pu = z_pu / sqrt(1 + x_r_ratio^2)
        x1_pu = r1_pu * x_r_ratio
        
        # Calculate negative sequence impedance parameters
        # For transformers, negative sequence impedance usually equals positive sequence impedance
        r2_pu = r1_pu
        x2_pu = x1_pu
        
        # If transformer has specific negative sequence parameters, use them:
        if hasfield(typeof(transformer), :z2_percent)
            z2_pu = transformer.z2_percent / 100.0 / s_ratio
            if hasfield(typeof(transformer), :x2_r)
                x2_r_ratio = transformer.x2_r
            else
                x2_r_ratio = x_r_ratio  # Use positive sequence X/R ratio
            end
            r2_pu = z2_pu / sqrt(1 + x2_r_ratio^2)
            x2_pu = r2_pu * x2_r_ratio
        end
        
        # Consider parallel transformers
        parallel = transformer.parallel
        if parallel > 1
            r1_pu = r1_pu / parallel
            x1_pu = x1_pu / parallel
            r2_pu = r2_pu / parallel
            x2_pu = x2_pu / parallel
        end
        
        # Fill positive sequence branch matrix (branchAC_1)
        branch_1[i, F_BUS] = hv_bus_idx
        branch_1[i, T_BUS] = lv_bus_idx
        branch_1[i, BR_R] = r1_pu
        branch_1[i, BR_X] = x1_pu
        branch_1[i, BR_B] = 0.0  # Transformers usually have no shunt susceptance
        
        # Fill negative sequence branch matrix (branchAC_2)
        branch_2[i, F_BUS] = hv_bus_idx
        branch_2[i, T_BUS] = lv_bus_idx
        branch_2[i, BR_R] = r2_pu
        branch_2[i, BR_X] = x2_pu
        branch_2[i, BR_B] = 0.0  # Transformers usually have no shunt susceptance
        
        # Set tap ratio and phase shift (same for positive and negative sequence)
        tap_ratio = 1.0
        shift_angle = 0.0
        
        # If transformer has tap ratio information
        if hasfield(typeof(transformer), :tap_ratio)
            tap_ratio = transformer.tap_ratio
        end
        if hasfield(typeof(transformer), :shift_angle)
            shift_angle = transformer.shift_angle
        end
        
        branch_1[i, TAP] = tap_ratio
        branch_1[i, SHIFT] = shift_angle
        branch_2[i, TAP] = tap_ratio
        branch_2[i, SHIFT] = shift_angle
        
        # Set rated capacity (same for positive and negative sequence)
        rate_a = transformer.sn_mva
        branch_1[i, RATE_A] = rate_a
        branch_2[i, RATE_A] = rate_a
        
        # Set branch status (same for positive and negative sequence)
        status = transformer.in_service ? 1.0 : 0.0
        branch_1[i, BR_STATUS] = status
        branch_2[i, BR_STATUS] = status
        
        # Set angle limits (same for positive and negative sequence)
        branch_1[i, ANGMIN] = -360.0
        branch_1[i, ANGMAX] = 360.0
        branch_2[i, ANGMIN] = -360.0
        branch_2[i, ANGMAX] = 360.0
    end
    
    # Add transformer branch data to three-phase JPC structure
    if isempty(jpc_3ph.branchAC_1)
        jpc_3ph.branchAC_1 = branch_1
        jpc_3ph.branchAC_2 = branch_2
    else
        jpc_3ph.branchAC_1 = [jpc_3ph.branchAC_1; branch_1]
        jpc_3ph.branchAC_2 = [jpc_3ph.branchAC_2; branch_2]
    end
end

"""
    calculate_branch_JPC_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Calculate zero sequence branch parameters for lines and transformers.

This function creates the zero sequence branch matrix and populates it with parameters
for lines and transformers. Zero sequence parameters often differ significantly from
positive sequence parameters, especially for transformers with different winding connections.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Returns
- `jpc_3ph::Utils.JPC_3ph`: The updated JPC_3ph object with zero sequence branch data
"""
function calculate_branch_JPC_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Initialize branch data matrix
    nbr = length(case.branchesAC)
    branch = zeros(nbr, 14)
    
    # Set default values
    branch[:, RATE_A] .= 250.0
    branch[:, RATE_B] .= 250.0
    branch[:, RATE_C] .= 250.0
    branch[:, TAP] .= 1.0
    branch[:, SHIFT] .= 0.0
    branch[:, BR_STATUS] .= 1.0
    branch[:, ANGMIN] .= -360.0
    branch[:, ANGMAX] .= 360.0
    
    # Add line zero sequence impedance
    add_line_sc_impedance_zero(case, jpc_3ph, branch)
    
    # Add transformer zero sequence impedance
    branch = add_trafo_sc_impedance_zero(case, jpc_3ph, branch)
    
    # Store calculation results in JPC object
    jpc_3ph.branchAC_0 = branch
    
    return jpc_3ph
end

"""
    add_line_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)

Add zero sequence impedance parameters for lines to the branch matrix.

This function calculates and adds zero sequence parameters for transmission lines.
If specific zero sequence data is not available, it approximates using positive sequence data.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object
- `branch`: The branch matrix to be updated with zero sequence parameters

# Returns
- Updated branch matrix with zero sequence line parameters
"""
function add_line_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)
    # Check if there is line data
    if isempty(case.branchesAC)
        return
    end
    
    # Process all lines (excluding transformers)
    line_indices = findall(l -> !hasfield(typeof(l), :is_transformer) || !l.is_transformer, case.branchesAC)
    
    if isempty(line_indices)
        return
    end
    
    for i in line_indices
        line = case.branchesAC[i]
        
        # Get from and to bus numbers
        fb = line.from_bus
        tb = line.to_bus
        
        # Get length and parallel data
        length_km = line.length_km
        parallel = hasfield(typeof(line), :parallel) ? line.parallel : 1.0
        
        # Calculate base impedance (note division by 3, consistent with original code)
        base_kv = jpc_3ph.busAC_0[fb, BASE_KV]
        baseR = (base_kv^2) / (3 * case.baseMVA)
        
        # Check if zero sequence parameters exist
        if hasfield(typeof(line), :r0_ohm_per_km) && hasfield(typeof(line), :x0_ohm_per_km) && hasfield(typeof(line), :c0_nf_per_km)
            # Calculate zero sequence impedance
            r0_pu = length_km * line.r0_ohm_per_km / baseR / parallel
            x0_pu = length_km * line.x0_ohm_per_km / baseR / parallel
            
            # Calculate zero sequence shunt susceptance
            b0_pu = 2 * π * case.basef * length_km * line.c0_nf_per_km * 1e-9 * baseR * parallel
        else
            # If no zero sequence data, use approximation from positive sequence data
            r0_pu = length_km * line.r_ohm_per_km * 3.0 / baseR / parallel
            x0_pu = length_km * line.x_ohm_per_km * 3.0 / baseR / parallel
            b0_pu = 0.0
        end
        
        # Fill branch matrix
        branch[i, F_BUS] = fb
        branch[i, T_BUS] = tb
        branch[i, BR_R] = r0_pu
        branch[i, BR_X] = x0_pu
        branch[i, BR_B] = b0_pu
        branch[i, BR_STATUS] = line.in_service ? 1.0 : 0.0
    end
    return branch  # Return updated branch matrix
end

"""
    add_trafo_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)

Add zero sequence impedance parameters for transformers to the branch matrix.

This function calculates and adds zero sequence parameters for two-winding transformers,
considering different transformer vector groups (winding connections). For some connections,
transformers are represented as shunt admittances rather than branches in the zero sequence network.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object
- `branch`: The branch matrix to be updated with zero sequence parameters

# Returns
- Updated branch matrix with zero sequence transformer parameters
"""
function add_trafo_sc_impedance_zero(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, branch)
    # Process zero sequence impedance for two-winding ETAP transformers
    transformers = case.transformers_2w_etap
    
    if isempty(transformers)
        return branch  # If no transformers, return original branch
    end
    
    # Get number of lines to determine starting index for transformers in branch matrix
    n_lines = count(l -> !hasfield(typeof(l), :is_transformer) || !l.is_transformer, case.branchesAC)
    
    # Check if branch matrix has enough rows
    n_transformers = length(transformers)
    current_rows = size(branch, 1)
    
    if current_rows < n_lines + n_transformers
        # If branch matrix doesn't have enough rows, extend it
        additional_rows = n_lines + n_transformers - current_rows
        branch = vcat(branch, zeros(additional_rows, size(branch, 2)))
    end
    
    # Create arrays to store bus, conductance and susceptance
    buses_all = Int[]
    gs_all = Float64[]
    bs_all = Float64[]
    
    # Process each transformer
    for (i, transformer) in enumerate(transformers)
        branch_idx = n_lines + i  # Index in branch matrix
        
        # Get high voltage and low voltage bus numbers
        hv_bus = transformer.hv_bus
        lv_bus = transformer.lv_bus
        
        # Get vector group
        vector_group = transformer.vector_group
        
        # Skip unsupported vector groups
        if lowercase(vector_group) in ["yy", "yd", "dy", "dd"]
            continue
        end
        
        # Set basic parameters in branch matrix
        branch[branch_idx, F_BUS] = hv_bus
        branch[branch_idx, T_BUS] = lv_bus
        
        # By default, transformer doesn't participate in zero sequence network
        branch[branch_idx, BR_STATUS] = 0.0
        
        # Get zero sequence impedance parameters (if no zero sequence parameters, use positive sequence)
        z0_percent = hasfield(typeof(transformer), :z0_percent) ? transformer.z0_percent : transformer.z_percent
        x0_r0 = hasfield(typeof(transformer), :x0_r0) ? transformer.x0_r0 : transformer.x_r
        
        # Get rated power (MVA) and parallel count
        sn_trafo_mva = transformer.sn_mva
        parallel = transformer.parallel
        
        # Calculate zero sequence impedance (similar to positive sequence, but using zero sequence parameters)
        s_ratio = sn_trafo_mva / case.baseMVA
        z0_pu = z0_percent / 100.0 / s_ratio  # Convert to system base
        
        r0_pu = z0_pu / sqrt(1 + x0_r0^2)
        x0_pu = r0_pu * x0_r0
        
        # Consider parallel transformers
        if parallel > 1
            r0_pu = r0_pu / parallel
            x0_pu = x0_pu / parallel
        end
        
        # Process zero sequence network based on transformer connection type
        if lowercase(vector_group) == "ynyn"
            # YNyn type transformer appears as a normal branch in zero sequence network
            branch[branch_idx, BR_R] = r0_pu
            branch[branch_idx, BR_X] = x0_pu
            branch[branch_idx, BR_B] = 0.0
            branch[branch_idx, BR_STATUS] = transformer.in_service ? 1.0 : 0.0
            
        elseif lowercase(vector_group) == "dyn"
            # Dyn type transformer in zero sequence network: HV side doesn't participate, LV side grounded
            # Add equivalent impedance as shunt admittance to LV bus
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, lv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "ynd"
            # YNd type transformer in zero sequence network: LV side doesn't participate, HV side grounded
            # Add equivalent impedance as shunt admittance to HV bus
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, hv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
                elseif lowercase(vector_group) == "yyn"
            # YYn type transformer in zero sequence network: add equivalent impedance as shunt admittance to LV bus
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, lv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "yny"
            # YNy type transformer in zero sequence network: add equivalent impedance as shunt admittance to HV bus
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                push!(buses_all, hv_bus)
                push!(gs_all, real(y0))
                push!(bs_all, imag(y0))
            end
            
        elseif lowercase(vector_group) == "yzn"
            # YZn type transformer in zero sequence network: requires special impedance relationship
            # Usually add equivalent impedance multiplied by a coefficient as shunt admittance to LV bus
            if transformer.in_service
                y0 = 1.0 / Complex(r0_pu, x0_pu)
                # Coefficient 1.1547 = sqrt(3)/sqrt(2), related to transformer connection
                push!(buses_all, lv_bus)
                push!(gs_all, 1.1547 * real(y0))
                push!(bs_all, 1.1547 * imag(y0))
            end
            
        else
            # Unsupported vector group
            @warn "Transformer vector group $(vector_group) not supported or not implemented, transformer index: $i"
        end
        
        # Set transformer tap ratio and phase shift (same as positive sequence)
        # Get transformer rated voltages
        vn_trafo_hv = transformer.vn_hv_kv
        vn_trafo_lv = transformer.vn_lv_kv
        
        # Get bus rated voltages
        vn_bus_hv = jpc_3ph.busAC_0[hv_bus, BASE_KV]
        vn_bus_lv = jpc_3ph.busAC_0[lv_bus, BASE_KV]
        
        # Calculate tap ratio
        ratio = (vn_trafo_hv / vn_bus_hv) / (vn_trafo_lv / vn_bus_lv)
        
        # Get phase shift angle
        shift = hasfield(typeof(transformer), :shift_degree) ? transformer.shift_degree : 0.0
        
        # Set transformer tap ratio and phase shift
        branch[branch_idx, TAP] = ratio
        branch[branch_idx, SHIFT] = shift
        
        # Set rated capacity and angle limits (same as positive sequence)
        branch[branch_idx, RATE_A] = case.baseMVA 
        branch[branch_idx, ANGMIN] = -360.0
        branch[branch_idx, ANGMAX] = 360.0
    end
    
    # Add conductance and susceptance to bus matrix
    # Combine values for the same bus
    if !isempty(buses_all)
        buses_unique = unique(buses_all)
        for bus in buses_unique
            indices = findall(x -> x == bus, buses_all)
            gs_sum = sum(gs_all[indices])
            bs_sum = sum(bs_all[indices])
            
            jpc_3ph.busAC_0[bus, GS] += gs_sum
            jpc_3ph.busAC_0[bus, BS] += bs_sum
        end
    end
    
    return branch  # Return updated branch matrix
end

"""
    JPC_3ph_gens_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Process generator data from a JuliaPowerCase and populate the JPC_3ph object with sequence component generator data.

This function handles different types of generators including conventional generators, static generators,
and external grids. It also sets appropriate bus types (PQ, PV, or reference) based on the generators.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Side effects
- Updates `jpc_3ph.genAC_1`, `jpc_3ph.genAC_2`, and `jpc_3ph.genAC_0` with generator data
- Updates bus types in `jpc_3ph.busAC_1`, `jpc_3ph.busAC_2`, and `jpc_3ph.busAC_0`
"""
function JPC_3ph_gens_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Count different types of generation devices
    n_gen = length(case.gensAC)
    n_sgen = length(case.sgensAC)
    n_ext = length(case.ext_grids)
    
    # Calculate total number of generation devices
    total_gens = n_gen + n_sgen + n_ext
    
    if total_gens == 0
        return  # If no generation devices, return directly
    end
    
    # Create generator matrix with rows equal to number of generation devices and 26 columns
    gen_data = zeros(total_gens, 26)
    
    # Process external grids (usually as slack/reference nodes)
    for (i, ext) in enumerate(case.ext_grids)
        if !ext.in_service
            continue
        end
        
        bus_idx = ext.bus
        
        # Fill generator data
        gen_data[i, :] = [
            bus_idx,        # Bus number generator is connected to
            0.0,            # Active power output (MW)
            0.0,            # Reactive power output (MVAr)
            9999.0,         # Maximum reactive power output (MVAr)
            -9999.0,        # Minimum reactive power output (MVAr)
            ext.vm_pu,      # Voltage magnitude setpoint (p.u.)
            case.baseMVA,   # Generator base MVA
            1.0,            # Generator status (1=on, 0=off)
            9999.0,         # Maximum active power output (MW)
            -9999.0,        # Minimum active power output (MW)
            0.0,            # Active power output at lower end of PQ capability curve (MW)
            0.0,            # Active power output at upper end of PQ capability curve (MW)
            0.0,            # Minimum reactive power output at PC1 (MVAr)
            0.0,            # Maximum reactive power output at PC1 (MVAr)
            0.0,            # Minimum reactive power output at PC2 (MVAr)
            0.0,            # Maximum reactive power output at PC2 (MVAr)
            0.0,            # AGC ramp rate (MW/min)
            0.0,            # 10-minute reserve ramp rate (MW)
            0.0,            # 30-minute reserve ramp rate (MW)
            0.0,            # Reactive power ramp rate (MVAr/min)
            1.0,            # Area participation factor
            2.0,            # Generator model (2=polynomial cost model)
            0.0,            # Startup cost ($)
            0.0,            # Shutdown cost ($)
            3.0,            # Number of polynomial cost function coefficients
            0.0             # Cost function parameters (to be extended later)
        ]
        
        # Update bus type to reference node (REF/slack) - set for all sequence components
        jpc_3ph.busAC_1[bus_idx, 2] = 3  # 3 means REF node
        jpc_3ph.busAC_2[bus_idx, 2] = 3  # 3 means REF node
        jpc_3ph.busAC_0[bus_idx, 2] = 3  # 3 means REF node
    end
    
    # Process conventional generators (usually as PV nodes)
    offset = n_ext
    for (i, gen) in enumerate(case.gensAC)
        if !gen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = gen.bus
        
        # Calculate reactive power (if not directly given)
        q_mvar = 0.0
        if hasfield(typeof(gen), :q_mvar)
            q_mvar = gen.q_mvar
        else
            # Calculate reactive power from power factor
            p_mw = gen.p_mw * gen.scaling
            if gen.cos_phi > 0 && p_mw > 0
                q_mvar = p_mw * tan(acos(gen.cos_phi))
            end
        end
        
        # Base MVA
        mbase = gen.sn_mva > 0 ? gen.sn_mva : case.baseMVA
        
        # Ramp rate parameters
        ramp_agc = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                   gen.ramp_up_rate_mw_per_min : 
                   (gen.max_p_mw - gen.min_p_mw) / 10
        ramp_10 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 10 : 
                  gen.max_p_mw - gen.min_p_mw
        ramp_30 = hasfield(typeof(gen), :ramp_up_rate_mw_per_min) ? 
                  gen.ramp_up_rate_mw_per_min * 30 : 
                  gen.max_p_mw - gen.min_p_mw
        
        # Fill generator data
        gen_data[idx, :] = [
            bus_idx,                               # Bus number generator is connected to
            gen.p_mw * gen.scaling,                # Active power output (MW)
            q_mvar,                                # Reactive power output (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power output (MVAr)
            gen.min_q_mvar,                        # Minimum reactive power output (MVAr)
            gen.vm_pu,                             # Voltage magnitude setpoint (p.u.)
            mbase,                                 # Generator base MVA
            1.0,                                   # Generator status (1=on, 0=off)
            gen.max_p_mw,                          # Maximum active power output (MW)
            gen.min_p_mw,                          # Minimum active power output (MW)
            gen.min_p_mw,                          # Active power at lower end of PQ capability curve (MW)
            gen.max_p_mw,                          # Active power at upper end of PQ capability curve (MW)
            gen.min_q_mvar,                        # Minimum reactive power at PC1 (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power at PC1 (MVAr)
            gen.min_q_mvar,                        # Minimum reactive power at PC2 (MVAr)
            gen.max_q_mvar,                        # Maximum reactive power at PC2 (MVAr)
            ramp_agc,                              # AGC ramp rate (MW/min)
            ramp_10,                               # 10-minute reserve ramp rate (MW)
            ramp_30,                               # 30-minute reserve ramp rate (MW)
            (gen.max_q_mvar - gen.min_q_mvar)/10,  # Reactive power ramp rate (MVAr/min)
            1.0,                                   # Area participation factor
            2.0,                                   # Generator model (2=polynomial cost model)
            0.0,                                   # Startup cost ($)
            0.0,                                   # Shutdown cost ($)
            3.0,                                   # Number of polynomial cost function coefficients
            0.0                                    # Cost function parameters (to be extended later)
        ]
        
        # If bus is not already set as reference node, set it as PV node - for all sequence components
        if jpc_3ph.busAC_1[bus_idx, 2] != 3  # 3 means REF node
            jpc_3ph.busAC_1[bus_idx, 2] = 2  # 2 means PV node
        end
        if jpc_3ph.busAC_2[bus_idx, 2] != 3
            jpc_3ph.busAC_2[bus_idx, 2] = 2
        end
        if jpc_3ph.busAC_0[bus_idx, 2] != 3
            jpc_3ph.busAC_0[bus_idx, 2] = 2
        end
    end
    
    # Process static generators (usually as PQ nodes, but could be PV nodes if they have voltage control capability)
    offset = n_ext + n_gen
    for (i, sgen) in enumerate(case.sgensAC)
        if !sgen.in_service
            continue
        end
        
        idx = i + offset
        bus_idx = sgen.bus
        
        # Fill generator data
        gen_data[idx, :] = [
            bus_idx,                                # Bus number generator is connected to
            sgen.p_mw * sgen.scaling,               # Active power output (MW)
            sgen.q_mvar * sgen.scaling,             # Reactive power output (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power output (MVAr)
            sgen.min_q_mvar,                        # Minimum reactive power output (MVAr)
            1.0,                                    # Voltage magnitude setpoint (p.u.)
            case.baseMVA,                           # Generator base MVA
            1.0,                                    # Generator status (1=on, 0=off)
            sgen.max_p_mw,                          # Maximum active power output (MW)
            sgen.min_p_mw,                          # Minimum active power output (MW)
            sgen.min_p_mw,                          # Active power at lower end of PQ capability curve (MW)
            sgen.max_p_mw,                          # Active power at upper end of PQ capability curve (MW)
            sgen.min_q_mvar,                        # Minimum reactive power at PC1 (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power at PC1 (MVAr)
            sgen.min_q_mvar,                        # Minimum reactive power at PC2 (MVAr)
            sgen.max_q_mvar,                        # Maximum reactive power at PC2 (MVAr)
            (sgen.max_p_mw - sgen.min_p_mw) / 10,   # AGC ramp rate (MW/min)
            sgen.max_p_mw - sgen.min_p_mw,          # 10-minute reserve ramp rate (MW)
            sgen.max_p_mw - sgen.min_p_mw,          # 30-minute reserve ramp rate (MW)
            (sgen.max_q_mvar - sgen.min_q_mvar)/10, # Reactive power ramp rate (MVAr/min)
            1.0,                                    # Area participation factor
            2.0,                                    # Generator model (2=polynomial cost model)
            0.0,                                    # Startup cost ($)
            0.0,                                    # Shutdown cost ($)
            3.0,                                    # Number of polynomial cost function coefficients
            0.0                                     # Cost function parameters (to be extended later)
        ]
        
        # If static generator is controllable and bus is not already set as REF or PV node, set it as PV node
        if sgen.controllable
            if jpc_3ph.busAC_1[bus_idx, 2] == 1  # 1 means PQ node
                jpc_3ph.busAC_1[bus_idx, 2] = 2  # 2 means PV node
            end
            if jpc_3ph.busAC_2[bus_idx, 2] == 1
                jpc_3ph.busAC_2[bus_idx, 2] = 2
            end
            if jpc_3ph.busAC_0[bus_idx, 2] == 1
                jpc_3ph.busAC_0[bus_idx, 2] = 2
            end
        end
    end
    
    # Remove unused rows (corresponding to out-of-service generation devices)
    active_rows = findall(x -> x > 0, gen_data[:, 8])  # Column 8 is GEN_STATUS
    gen_data = gen_data[active_rows, :]
    
    # Store generator data in three-phase JPC structure for all sequence components
    jpc_3ph.genAC_1 = copy(gen_data)  # Positive sequence
    jpc_3ph.genAC_2 = copy(gen_data)  # Negative sequence
    jpc_3ph.genAC_0 = copy(gen_data)  # Zero sequence
    
    # Ensure at least one slack bus - set for all sequence components
    if !any(jpc_3ph.busAC_1[:, 2] .== 3) && size(gen_data, 1) > 0  # 3 means REF node
        # If no slack bus, select the first generator bus as slack
        first_gen_bus = Int(gen_data[1, 1])
        jpc_3ph.busAC_1[first_gen_bus, 2] = 3  # 3 means REF node
        jpc_3ph.busAC_2[first_gen_bus, 2] = 3
        jpc_3ph.busAC_0[first_gen_bus, 2] = 3
    end
end

"""
    JPC_3ph_loads_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)

Process load data from a JuliaPowerCase and populate the JPC_3ph object with sequence component load data.

This function handles different types of loads (wye and delta connected) and distributes
their power appropriately across sequence components. For balanced loads, only positive
sequence components have non-zero values.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object to be populated

# Side effects
- Updates `jpc_3ph.loadAC_1`, `jpc_3ph.loadAC_2`, and `jpc_3ph.loadAC_0` with load data
- Updates PD and QD fields in bus matrices for all sequence components
"""
function JPC_3ph_loads_process(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph)
    # Process load data, convert to three-phase JPC format and update busAC's PD and QD
    
    # Filter in-service loads
    in_service_loads = filter(load -> load.in_service == true, case.loadsAC)
    
    # If no in-service loads, return directly
    if isempty(in_service_loads)
        return
    end
    
    # Create an empty matrix with rows equal to number of loads and 8 columns
    num_loads = length(in_service_loads)
    load_matrix = zeros(num_loads, 8)
    
    # Create three dictionaries to accumulate loads connected to the same bus (three sequence components)
    bus_load_sum_1 = Dict{Int, Vector{Float64}}()  # Positive sequence
    bus_load_sum_2 = Dict{Int, Vector{Float64}}()  # Negative sequence
    bus_load_sum_0 = Dict{Int, Vector{Float64}}()  # Zero sequence
    
    for (i, load) in enumerate(in_service_loads)
        # Calculate actual active and reactive load (considering scaling factor)
        actual_p_mw = load.p_mw * load.scaling
        actual_q_mvar = load.q_mvar * load.scaling
        
        # Distribute load to sequence components based on load type
        if load.type == "wye"
            # Y-connected: balanced three-phase load, mainly positive sequence
            p_1 = actual_p_mw    # Positive sequence component (per phase power)
            q_1 = actual_q_mvar
            p_2 = 0.0                  # Negative sequence component (zero for balanced load)
            q_2 = 0.0
            p_0 = 0.0                  # Zero sequence component (zero for balanced load)
            q_0 = 0.0
        elseif load.type == "delta"
            # Delta-connected: line voltage, zero sequence component is zero
            p_1 = actual_p_mw    # Positive sequence component
            q_1 = actual_q_mvar
            p_2 = 0.0                  # Negative sequence component (zero for balanced load)
            q_2 = 0.0
            p_0 = 0.0                  # Zero sequence component (no zero sequence current in delta connection)
            q_0 = 0.0
        else
            # Default handling as Y-connected
            p_1 = actual_p_mw
            q_1 = actual_q_mvar
            p_2 = 0.0
            q_2 = 0.0
            p_0 = 0.0
            q_0 = 0.0
        end
        
        # Fill basic load matrix (using positive sequence component)
        load_matrix[i, :] = [
            i,                         # Load number
            load.bus,                  # Bus number load is connected to
            1.0,                       # Load status (1=in service)
            p_1,                       # Positive sequence active load (MW)
            q_1,                       # Positive sequence reactive load (MVAr)
            load.const_z_percent/100,  # Constant impedance load percentage
            load.const_i_percent/100,  # Constant current load percentage
            load.const_p_percent/100   # Constant power load percentage
        ]
        
        # Accumulate loads connected to the same bus (positive sequence)
        bus_idx = load.bus
        if haskey(bus_load_sum_1, bus_idx)
            bus_load_sum_1[bus_idx][1] += p_1
            bus_load_sum_1[bus_idx][2] += q_1
        else
            bus_load_sum_1[bus_idx] = [p_1, q_1]
        end
        
        # Accumulate loads connected to the same bus (negative sequence)
        if haskey(bus_load_sum_2, bus_idx)
            bus_load_sum_2[bus_idx][1] += p_2
            bus_load_sum_2[bus_idx][2] += q_2
        else
            bus_load_sum_2[bus_idx] = [p_2, q_2]
        end
        
        # Accumulate loads connected to the same bus (zero sequence)
        if haskey(bus_load_sum_0, bus_idx)
            bus_load_sum_0[bus_idx][1] += p_0
            bus_load_sum_0[bus_idx][2] += q_0
        else
            bus_load_sum_0[bus_idx] = [p_0, q_0]
        end
    end
    
    # Create load matrices for three sequence components
    load_matrix_1 = copy(load_matrix)  # Positive sequence load matrix
    load_matrix_2 = copy(load_matrix)  # Negative sequence load matrix
    load_matrix_0 = copy(load_matrix)  # Zero sequence load matrix
    
    # Update power values in negative and zero sequence load matrices
    for (i, load) in enumerate(in_service_loads)
        # Update negative sequence matrix
        load_matrix_2[i, 4] = 0.0  # Balanced load's negative sequence active power is 0
        load_matrix_2[i, 5] = 0.0  # Balanced load's negative sequence reactive power is 0
        
        # Update zero sequence matrix
        load_matrix_0[i, 4] = 0.0  # Balanced load's zero sequence active power is 0
        load_matrix_0[i, 5] = 0.0  # Balanced load's zero sequence reactive power is 0
    end
    
    # Store load data in three-phase JPC structure for three sequence components
    jpc_3ph.loadAC_1 = load_matrix_1  # Positive sequence
    jpc_3ph.loadAC_2 = load_matrix_2  # Negative sequence
    jpc_3ph.loadAC_0 = load_matrix_0  # Zero sequence
    
    # Update PD and QD fields in busAC matrices for three sequence components
    
    # Update positive sequence busAC_1
    for (bus_idx, load_values) in bus_load_sum_1
        # Find corresponding bus row
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_1[:, 1])
        
        if !isnothing(bus_row)
            # Update PD (column 3) and QD (column 4)
            jpc_3ph.busAC_1[bus_row, PD] = load_values[1]  # PD - Active load (MW)
            jpc_3ph.busAC_1[bus_row, QD] = load_values[2]  # QD - Reactive load (MVAr)
        end
    end
    
    # Update negative sequence busAC_2
    for (bus_idx, load_values) in bus_load_sum_2
        # Find corresponding bus row
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_2[:, 1])
        
        if !isnothing(bus_row)
            # Update PD (column 3) and QD (column 4)
            jpc_3ph.busAC_2[bus_row, PD] = load_values[1]  # PD - Active load (MW)
            jpc_3ph.busAC_2[bus_row, QD] = load_values[2]  # QD - Reactive load (MVAr)
        end
    end
    
    # Update zero sequence busAC_0
    for (bus_idx, load_values) in bus_load_sum_0
        # Find corresponding bus row
        bus_row = findfirst(x -> x == bus_idx, jpc_3ph.busAC_0[:, 1])
        
        if !isnothing(bus_row)
            # Update PD (column 3) and QD (column 4)
            jpc_3ph.busAC_0[bus_row, PD] = load_values[1]  # PD - Active load (MW)
            jpc_3ph.busAC_0[bus_row, QD] = load_values[2]  # QD - Reactive load (MVAr)
        end
    end
end

"""
    JPC_3ph_add_grid_external_sc_impedance(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, sequence::Int)

Add external grid short circuit impedance to the appropriate sequence component.

This function calculates and adds short circuit impedance for external grids to the
specified sequence component. For positive sequence (sequence=1), no additional impedance
is added. For negative (sequence=2) and zero (sequence=0) sequences, short circuit impedance
is added as shunt admittance to the connected buses.

# Arguments
- `case::JuliaPowerCase`: The power system case data
- `jpc_3ph::Utils.JPC_3ph`: The three-phase JPC object
- `sequence::Int`: Sequence component identifier (0=zero, 1=positive, 2=negative)

# Returns
- `jpc_3ph::Utils.JPC_3ph`: The updated JPC_3ph object
- `gs_eg::Vector{Float64}`: Conductance values for external grids (only for sequence=2)
- `bs_eg::Vector{Float64}`: Susceptance values for external grids (only for sequence=2)
"""
function JPC_3ph_add_grid_external_sc_impedance(case::JuliaPowerCase, jpc_3ph::Utils.JPC_3ph, sequence::Int)    
    # If no external grids, return empty arrays and original jpc_3ph
    if isempty(case.ext_grids)
        return jpc_3ph, Float64[], Float64[]
    end
    
    # Only process external grid impedance for negative (sequence=2) or zero (sequence=0) sequence
    if sequence == 2 || sequence == 0
        # Collect data for all external grids
        external_buses = Int[]
        Y_grid_real = Float64[]
        Y_grid_imag = Float64[]
        
        # Iterate through all external grids
        for ext_grid in case.ext_grids
            # Get bus ID the external grid is connected to
            external_bus = ext_grid.bus
            
            # Calculate external grid short circuit impedance
            c = 1.1
            s_sc = ext_grid.s_sc_max_mva / jpc_3ph.baseMVA
            rx = ext_grid.rx_max
            z_grid = c / (s_sc / 3)
            x_grid = z_grid / sqrt(1 + rx^2)
            r_grid = x_grid * rx
            
            # Calculate admittance
            Y_grid = 1 / (r_grid + 1im * x_grid)
            
            push!(external_buses, external_bus)
            push!(Y_grid_real, real(Y_grid))
            push!(Y_grid_imag, imag(Y_grid))
        end
        
        # Sum admittance values for the same bus
        buses, gs, bs = PowerFlow.sum_by_group(external_buses, Y_grid_real, Y_grid_imag)
        
        # Update corresponding bus data based on sequence component
        if sequence == 2
            # Update negative sequence busAC_2 bus data
            for (i, bus_id) in enumerate(buses)
                # Find corresponding bus index in jpc_3ph.busAC_2
                bus_idx = findfirst(x -> x == bus_id, jpc_3ph.busAC_2[:, 1])
                
                if !isnothing(bus_idx)
                                        # Update bus conductance and susceptance values
                    jpc_3ph.busAC_2[bus_idx, GS] = gs[i] * jpc_3ph.baseMVA
                    jpc_3ph.busAC_2[bus_idx, BS] = bs[i] * jpc_3ph.baseMVA
                end
            end
        elseif sequence == 0
            # Update zero sequence busAC_0 bus data
            for (i, bus_id) in enumerate(buses)
                # Find corresponding bus index in jpc_3ph.busAC_0
                bus_idx = findfirst(x -> x == bus_id, jpc_3ph.busAC_0[:, 1])
                
                if !isnothing(bus_idx)
                    # Update bus conductance and susceptance values
                    jpc_3ph.busAC_0[bus_idx, GS] = gs[i] * jpc_3ph.baseMVA
                    jpc_3ph.busAC_0[bus_idx, BS] = bs[i] * jpc_3ph.baseMVA
                end
            end
        end
        
        # Return updated jpc_3ph and calculated conductance, susceptance values
        gs_eg = gs .* jpc_3ph.baseMVA
        bs_eg = bs .* jpc_3ph.baseMVA
        return jpc_3ph, gs_eg, bs_eg
    else
        # For positive sequence (sequence=1), don't process external grid impedance, return directly
        return jpc_3ph, Float64[], Float64[]
    end
end


                
