"""
    runupf(case, jpc_3ph, gs_eg, bs_eg, opt)

This function is used to run a local unbalanced power flow analysis on unsymmetrical load nodes

#Step1: find the unbalanced nodes in the system
#Step2: find the interface branches of the unbalanced nodes and balanced nodes
"""
function runupf(case, jpc_3ph, gs_eg, bs_eg, opt)
    jpc0, jpc1, jpc2 = extratct_jpc_3ph(jpc_3ph)
    jpc0, i2e0 = PowerFlow.ext2int(jpc0)
    jpc1, i2e1 = PowerFlow.ext2int(jpc1)
    jpc2, i2e2 = PowerFlow.ext2int(jpc2)
    _, bus0, branch0, gen0, _, _, _, _, _, _ = get_pf_variables_from_JPC(jpc0)
    baseMVA, bus1, branch1, gen1, ref, pv, pq, _, _, _ = get_pf_variables_from_JPC(jpc1)
    _, bus2, branch2, gen2, _, _, _, _, _, _ = get_pf_variables_from_JPC(jpc2)

    s_abc_delta, s_abc = load_mapping(case, jpc1)
    jpc0, jpc1, jpc2, y_0_pu, y_1_pu, y_2_pu, y_0_f, y_1_f, _, y_0_t, y_1_t, _ = get_y_bus(jpc0, jpc1, jpc2)

    nb = size(bus1, 1) # number of buses in the system
    v_012_it = vcat(
    [reshape(Vector{ComplexF64}(ppc["busAC"][:, VM] .* exp.(im .* deg2rad.(ppc["busAC"][:, VA]))), 1, nb)
     for ppc in (jpc0, jpc1, jpc2)]...)
    #Initialize
    # For Delta transformation:
    # Voltage changed from line-earth to line-line using V_T
    # s_abc/v_abc will now give line-line currents. This is converted to line-earth
    # current using I-T
    v_del_xfmn = [1 -1 0;
                 0 1 -1;
                 -1 0 1]

    i_del_xfmn = [1 0 -1;
                 -1 1 0;
                 0 -1 1]
    v_abc_it = sequence_to_phase(v_012_it)
    
    ## Iteration using Power mismatch criterion
    outer_tolerance_mva = 3e-8
    Count = 0

    # Initialize mismatch as boolean array
    s_mismatch = fill(true, (2,1))

    # Record start time
    t0 = time()

    # Iteration loop
    max_iteration = 30
    while any(s_mismatch .> outer_tolerance_mva) && Count < max_iteration
        # Calculate per unit power
        s_abc_pu = -s_abc ./ jpc1["baseMVA"]
        s_abc_delta_pu = -s_abc_delta ./ jpc1["baseMVA"]
        # Calculate star connection current
        i_abc_it_wye = conj.(s_abc_pu ./ v_abc_it)

        # Calculate delta connection current
        i_abc_it_delta = i_del_xfmn * conj.(s_abc_delta_pu ./ (v_del_xfmn * v_abc_it))

        # For buses with both delta and wye loads we need to sum of their currents
        # to sum up the currents
        i_abc_it = i_abc_it_wye + i_abc_it_delta
        i012_it = phase_to_sequence(i_abc_it)
        v1_for_s1 = v_012_it[2, :]
        i1_for_s1 = -i012_it[2, :]
        v0_pu_it = transpose(v_012_it[1,:])
        v2_pu_it = transpose(v_012_it[3,:])
        i0_pu_it = transpose(i012_it[1,:])
        i2_pu_it = transpose(i012_it[3,:])
        s1 = v1_for_s1 .* conj.(i1_for_s1)

        # Current used to find S1 Positive sequence power
        jpc1["busAC"][pq, PD] = real(s1[pq]) * jpc1["baseMVA"]
        jpc1["busAC"][pq, QD] = imag(s1[pq]) * jpc1["baseMVA"]

        #run newton raphson power flow
        jpc1, success, iterations = run_newton_raphson_pf(jpc1,opt)
        internal_vars_jpc1 = get_internal_variables(jpc1)
        Ybus,Yf,Yt = internal_vars_jpc1[:Ybus],internal_vars_jpc1[:Yf],internal_vars_jpc1[:Yt]
        V = internal_vars_jpc1[:V]
        bus,branch,gen = internal_vars_jpc1[:bus],internal_vars_jpc1[:branch],internal_vars_jpc1[:gen]
        bus, gen, branch = PowerFlow.pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq)
        #return result
        jpc1["busAC"] = bus
        jpc1["genAC"] = gen
        jpc1["branchAC"] = branch

        # Conduct Negative and Zero sequence power flow
        v0_pu_it = V_from_I(y_0_pu, i0_pu_it)
        v2_pu_it = V_from_I(y_2_pu, i2_pu_it)

        # Evaluate Positive Sequence Power Mismatch
        i1_from_v_it = vec(I1_from_V012(v_012_it, y_1_pu))
        s_from_voltage = S_from_VI_elementwise(v1_for_s1, i1_from_v_it)
        v1_pu_it = V1_from_jpc(jpc1)
        v_012_new = combine_X012(v0_pu_it, v1_pu_it, v2_pu_it)
        s_mismatch = abs.(abs.(s1[pq]) .- abs.(s_from_voltage[pq]))
        v_012_it = v_012_new
        v_abc_it = sequence_to_phase(v_012_it)
        Count += 1
    end

    success = (Count < max_iteration)
    jpc_3ph["success"] = success
    jpc_3ph["iterations"] = Count
    
    for ppc in [jpc0, jpc1, jpc2]
        ppc["success"] = success
    end
    # TODO: Add reference to paper to explain the following steps
    # This is required since the ext_grid power results are not correct if its
    # not done
    ref, pv, pq = PowerFlow.bustypes(jpc0["busAC"], jpc0["genAC"])
    jpc0["busAC"][ref, GS] -= gs_eg
    jpc0["busAC"][ref, BS] -= bs_eg
    y_0_pu, y_0_f, y_0_t = PowerFlow.makeYbus(jpc0["baseMVA"], jpc0["busAC"], jpc0["branchAC"])
    # revert the change, otherwise repeated calculation with recycled elements will fail
    jpc0["busAC"][ref, GS] += gs_eg
    jpc0["busAC"][ref, BS] += bs_eg
    # Bus, Branch, and Gen  power values
    bus0, gen0, branch0 = PowerFlow.pfsoln(baseMVA, bus0, gen0, branch0, y_0_pu, y_0_f, y_0_t,
                                vec(v_012_it[1, :]), ref, pv, pq)
    bus1, gen1, branch1 = PowerFlow.pfsoln(baseMVA, bus1, gen1, branch1, y_1_pu, y_1_f, y_1_t,
                                vec(v_012_it[2, :]), ref, pv, pq)
    bus2, gen2, branch2 = PowerFlow.pfsoln(baseMVA, bus2, gen2, branch2, y_1_pu, y_1_f, y_1_t,
                                vec(v_012_it[3, :]) , ref, pv, pq)

    #store the results in the jpc structure
    jpc0["busAC"] = bus0
    jpc0["genAC"] = gen0
    jpc0["branchAC"] = branch0

    jpc1["busAC"] = bus1
    jpc1["genAC"] = gen1
    jpc1["branchAC"] = branch1

    jpc2["busAC"] = bus2
    jpc2["genAC"] = gen2
    jpc2["branchAC"] = branch2

    i_012_res = current_from_voltage_results(y_0_pu, y_1_pu, v_012_it)
    s_012_res = S_from_VI_elementwise(v_012_it, i_012_res) * jpc1["baseMVA"]
    eg_idx = case.ext_grids[1].bus
    eg_bus_idx_ppc = Int64.(real(jpc1["genAC"][eg_idx, GEN_BUS]))

    jpc0["genAC"][eg_idx, PG] = real.(s_012_res[1, eg_bus_idx_ppc])
    jpc1["genAC"][eg_idx, PG] = real.(s_012_res[2, eg_bus_idx_ppc])
    jpc2["genAC"][eg_idx, PG] = real.(s_012_res[3, eg_bus_idx_ppc])
    jpc0["genAC"][eg_idx, QG] = imag.(s_012_res[1, eg_bus_idx_ppc])
    jpc1["genAC"][eg_idx, QG] = imag.(s_012_res[2, eg_bus_idx_ppc])
    jpc2["genAC"][eg_idx, QG] = imag.(s_012_res[3, eg_bus_idx_ppc])

    # net["jpc0"] = jpc0
    # net["jpc1"] = jpc1
    # net["jpc2"] = jpc2

    get_bus_v_results_3ph(jpc_3ph, jpc0, jpc1, jpc2)
    bus_pq = get_p_q_results_3ph(case, jpc_3ph)

    get_full_branch_zero(jpc_3ph)
    get_branch_results_3ph(jpc_3ph, jpc0, jpc1, jpc2, bus_pq)
    # get_gen_results_3ph(case, c_3ph, jpc0, jpc1, jpc2, bus_pq)
    get_bus_results_3ph(case, jpc_3ph, bus_pq)

    return jpc_3ph
end

"""
    load_mapping(case::JuliaPowerCase, jpc1::JPC)

Maps loads from the case to the power flow model, handling both wye and delta connections
for all three phases.
"""
function load_mapping(case::JuliaPowerCase, jpc1::JPC)

    params = Dict{String, Any}()    
    phases = ["a", "b", "c"]
    load_types = ["wye", "delta"]
    load_elements = ["loadsAC", "loadsAC_asymm", "gensAC", "asymmetric_sgen"]

    # Get number of buses
    nb = size(jpc1.busAC, 1)

    for phase in phases
        for typ in load_types
            # Initialize power arrays
            params["S$(phase)$(typ)"] = (jpc1.busAC[:, PD] .+
                                        jpc1.busAC[:, QD] .* 1im) .* 0
            params["p$(phase)$(typ)"] = Float64[]  # p values from loads/sgens
            params["q$(phase)$(typ)"] = Float64[]  # q values from loads/sgens
            params["P$(phase)$(typ)"] = Float64[]  # Aggregated Active Power
            params["Q$(phase)$(typ)"] = Float64[]  # Aggregated reactive Power
            params["b$(phase)$(typ)"] = Int64[]    # bus map for phases
            params["b$(typ)"] = Int64[]            # aggregated bus map(s_abc)
            
            for element in load_elements
                get_elements(params, case, element, phase, typ)
            end
            
            # Mapping constant power loads to buses
            if !isempty(params["b$(typ)"])
                params["b$(phase)$(typ)"] = params["b$(typ)"]
                
                params["b$(phase)$(typ)"], params["P$(phase)$(typ)"],
                temp_q = PowerFlow.sum_by_group(params["b$(phase)$(typ)"],
                                    params["p$(phase)$(typ)"],
                                    params["q$(phase)$(typ)"] .* 1im)
                
                params["Q$(phase)$(typ)"] = imag.(temp_q)
                
                params["S$(phase)$(typ)"][params["b$(phase)$(typ)"]] .= 
                    params["P$(phase)$(typ)"] .+ params["Q$(phase)$(typ)"] .* 1im
            end
        end
    end
    Sabc_del = vcat(transpose(params["Sadelta"]), transpose(params["Sbdelta"]), transpose(params["Scdelta"]))
    Sabc_wye = vcat(transpose(params["Sawye"]), transpose(params["Sbwye"]), transpose(params["Scwye"]))
    return Sabc_del, Sabc_wye
end

"""
    get_elements(params, case::JuliaPowerCase, element, phase, typ)

This function is used to get the elements of the load mapping.
Automatically skips elements that don't exist in case.
"""
function get_elements(params, case::JuliaPowerCase, element, phase, typ)
    # First check if the element exists in the case
    if !hasfield(JuliaPowerCase, Symbol(element))
        # If the element doesn't exist, return the original parameters without modification
        return params
    end

    sign = endswith(element, "sgen") ? -1 : 1
    elm = getfield(case, Symbol(element))
    
    # Check if elm is empty
    if isempty(elm)
        return params
    end
    
    # Process data based on element type
    if element == "loadsAC" || element == "gensAC"
        # Use accurate field names for Load structure
        for item in elm
            if item.in_service == true && item.type == typ
                # Get relevant data
                p_value = item.p_mw / 3 * item.scaling * sign
                q_value = item.q_mvar / 3 * item.scaling * sign
                bus_value = convert(Int64, item.bus)
                
                push!(params["p$(phase)$(typ)"], p_value)
                push!(params["q$(phase)$(typ)"], q_value)
                push!(params["b$(typ)"], bus_value)
            end
        end
    elseif startswith(element, "asymmetric")
        # Process asymmetric loads/generators
        for item in elm
            if item.in_service == true && item.type == typ
                # Select correct power field based on phase
                p_field = Symbol("p_$(phase)_mw")
                q_field = Symbol("q_$(phase)_mvar")
                
                if hasproperty(item, p_field) && hasproperty(item, q_field)
                    p_value = getproperty(item, p_field) * item.scaling * sign
                    q_value = getproperty(item, q_field) * item.scaling * sign
                    bus_value = convert(Int64, item.bus)
                    
                    push!(params["p$(phase)$(typ)"], p_value)
                    push!(params["q$(phase)$(typ)"], q_value)
                    push!(params["b$(typ)"], bus_value)
                end
            end
        end
    end
    
    return params
end

"""
    get_pf_variables_from_JPC(jpc::JPC)

Extract power flow variables from JPC structure.
"""
function get_pf_variables_from_JPC(jpc::JPC)
    # No need for deep copy since we're only reading data
    if isnothing(jpc)
        throw(ArgumentError("jpc cannot be nothing"))
    end

    # Get data from JPC structure
    baseMVA = jpc.baseMVA
    bus = jpc.busAC
    branch = jpc.branchAC
    gen = jpc.genAC

    # Get indices list for each bus type
    ref, pv, pq = PowerFlow.bustypes(bus, gen)

    # Generator information
    # Assuming GEN_STATUS is a constant representing the generator status column index
    GEN_STATUS = 8  # Adjust according to actual situation
    GEN_BUS = 1     # Adjust according to actual situation
    VG = 6          # Adjust according to actual situation
    
    on = Int.(findall(gen[:, GEN_STATUS] .> 0))  # Which generators are running
    gbus = Int.(gen[on, GEN_BUS])                # Which buses these generators are connected to

    # Initial state
    # Assuming VM and VA are constants representing voltage magnitude and angle column indices
    VM = 8  # Adjust according to actual situation
    VA = 9  # Adjust according to actual situation
    
    V0 = bus[:, VM] .* exp.(1im * pi/180 * bus[:, VA])
    V0[gbus] = gen[on, VG] ./ abs.(V0[gbus]) .* V0[gbus]

    return baseMVA, bus, branch, gen, ref, pv, pq, on, gbus, V0
end

"""
    sum_by_group(bus, first_val, second_val)

Group values by bus and sum them within each group.
"""
function sum_by_group(bus, first_val, second_val)
    # Create sorting index
    order = sortperm(bus)
    
    # Rearrange arrays according to sorting index
    sorted_bus = bus[order]
    sorted_first_val = first_val[order]
    sorted_second_val = second_val[order]
    
    # Find positions of unique elements
    n = length(sorted_bus)
    index = trues(n)
    for i in 1:(n-1)
        index[i] = sorted_bus[i+1] != sorted_bus[i]
    end
    
    # Calculate cumulative sum for first group of values
    cumsum_first = cumsum(sorted_first_val)
    
    # Calculate cumulative sum for second group of values
    cumsum_second = cumsum(sorted_second_val)
    
    # Extract unique buses and corresponding cumulative sums
    unique_buses = sorted_bus[index]
    cumsum_first_unique = cumsum_first[index]
    cumsum_second_unique = cumsum_second[index]
    
    # Calculate sum for each group
    result_first = similar(cumsum_first_unique)
    result_second = similar(cumsum_second_unique)
    
    # First element remains unchanged
    result_first[1] = cumsum_first_unique[1]
    result_second[1] = cumsum_second_unique[1]
    
    # Calculate differences to get sum for each group
    for i in 2:length(unique_buses)
        result_first[i] = cumsum_first_unique[i] - cumsum_first_unique[i-1]
        result_second[i] = cumsum_second_unique[i] - cumsum_second_unique[i-1]
    end
    
    return unique_buses, result_first, result_second
end

"""
    extratct_jpc_3ph(jpc_3ph::PowerFlow.JPC_3ph)

Extract individual sequence components from a three-phase JPC object.
"""
function extratct_jpc_3ph(jpc_3ph::PowerFlow.JPC_3ph)
    # Create copies of JPC0, JPC1, JPC2
    jpc0 = PowerFlow.JPC()
    jpc1 = PowerFlow.JPC()
    jpc2 = PowerFlow.JPC()

    # Copy attributes from jpc_3ph to jpc0, jpc1, jpc2
    # Assign to jpc0
    jpc0.busAC = jpc_3ph.busAC_0
    jpc0.branchAC = jpc_3ph.branchAC_0
    jpc0.genAC = jpc_3ph.genAC_0
    jpc0.loadAC = jpc_3ph.loadAC_0

    # Assign to jpc1
    jpc1.busAC = jpc_3ph.busAC_1
    jpc1.branchAC = jpc_3ph.branchAC_1
    jpc1.genAC = jpc_3ph.genAC_1
    jpc1.loadAC = jpc_3ph.loadAC_1

    # Assign to jpc2
    jpc2.busAC = jpc_3ph.busAC_2
    jpc2.branchAC = jpc_3ph.branchAC_2
    jpc2.genAC = jpc_3ph.genAC_2
    jpc2.loadAC = jpc_3ph.loadAC_2

    return jpc0, jpc1, jpc2
end

"""
    get_y_bus(jpc0, jpc1, jpc2)

Build admittance matrices for zero, positive, and negative sequence networks.
"""
function get_y_bus(jpc0, jpc1, jpc2)
        # build admittance matrices
        y_0_bus, y_0_f, y_0_t = PowerFlow.makeYbus(jpc0["baseMVA"], jpc0["busAC"], jpc0["branchAC"])
        y_1_bus, y_1_f, y_1_t = PowerFlow.makeYbus(jpc1["baseMVA"], jpc1["busAC"], jpc1["branchAC"])
        y_2_bus, y_2_f, y_2_t = PowerFlow.makeYbus(jpc2["baseMVA"], jpc2["busAC"], jpc2["branchAC"])

    return jpc0, jpc1, jpc2, y_0_bus, y_1_bus, y_2_bus, y_0_f, y_1_f, y_2_f, y_0_t, y_1_t, y_2_t
end

"""
    phase_shift_unit_operator(angle_deg)

Create a complex number representing a phase shift of the given angle in degrees.
"""
function phase_shift_unit_operator(angle_deg)
    return 1 * exp(im * deg2rad(angle_deg))
end

"""
    sequence_to_phase(X012)

Transform sequence components (zero, positive, negative) to phase components (a, b, c).
"""
function sequence_to_phase(X012)
    a = phase_shift_unit_operator(120)
    asq = phase_shift_unit_operator(-120)

    Tabc = [1 1 1;
            1 asq a;
            1 a asq]
    return Tabc * X012
end

"""
    phase_to_sequence(Xabc)

Transform phase components (a, b, c) to sequence components (zero, positive, negative).
"""
function phase_to_sequence(Xabc)
    a = phase_shift_unit_operator(120)
    asq = phase_shift_unit_operator(-120)  
    T012 = [1 1 1;
        1 a asq;
        1 asq a] ./ 3 
    return T012 * Xabc
end

"""
    run_newton_raphson_pf(jpc::JPC, opt)

Run Newton-Raphson power flow calculation.
"""
function run_newton_raphson_pf(jpc::JPC, opt)
    baseMVA, bus, branch, gen, ref, pv, pq, _, _, V0 = get_pf_variables_from_JPC(jpc)
    
    # Build grid model
    Ybus, Yf, Yt = PowerFlow.makeYbus(baseMVA, bus, branch)
    Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V0)
    
    # Execute Newton-Raphson power flow calculation
    V, success, iterations, norm_history = PowerFlow.newtonpf(
        baseMVA, bus, gen, Ybus, V0, ref, pv, pq, 
        opt["PF"]["PF_TOL"], opt["PF"]["PF_MAX_IT"], opt["PF"]["NR_ALG"]
    )
    
    # Update success status and iteration count in JPC object
    jpc.success = success
    jpc.iterationsAC = iterations
    
    # Store internal variables in global dictionary using JPC object's memory address as key
    jpc_id = objectid(jpc)
    _JPC_INTERNAL_STORAGE[jpc_id] = Dict(
        :bus => bus,
        :gen => gen,
        :branch => branch,
        :baseMVA => baseMVA,
        :V => V,
        :pv => pv,
        :pq => pq,
        :ref => ref,
        :Sbus => Sbus,
        :Ybus => Ybus,
        :Yf => Yf,
        :Yt => Yt,
    )
    
    return jpc, success, iterations
end

"""
    get_internal_variables(jpc::JPC)

Retrieve stored internal variables for a JPC object.
"""
function get_internal_variables(jpc::JPC)
    jpc_id = objectid(jpc)
    if haskey(_JPC_INTERNAL_STORAGE, jpc_id)
        return _JPC_INTERNAL_STORAGE[jpc_id]
    else
        return nothing
    end
end

"""
    V_from_I(Y, I)

Calculate voltage from current using V = Y\\I.
"""
function V_from_I(Y, I)
    # Ensure I is in column vector form
    I_vec = vec(collect(I))  # Convert any form of I to column vector
    
    # Solve linear equation system
    V = Y \ I_vec
    
    # Return result
    return transpose(V)
end

"""
    I1_from_V012(V012, Y)

Calculate positive sequence current from sequence voltages and admittance matrix.
"""
function I1_from_V012(V012, Y)
    # Extract positive sequence component from symmetrical components
    V1 = reshape(transpose(V012[2, :]), :, 1)  # Convert to column vector
    
    # Check if Y is a sparse matrix
    if issparse(Y)
        # If sparse, convert to dense matrix before calculation
        i1 = Array(Matrix(Y) * V1)
        return transpose(i1)
    else
        # If dense, calculate directly
        i1 = Array(Y * V1)
        return transpose(i1)
    end
end

"""
    I0_from_V012(V012, Y)

Calculate zero sequence current from sequence voltages and admittance matrix.
"""
function I0_from_V012(V012, Y)
    # Extract zero sequence component from symmetrical components
    V0 = reshape(transpose(V012[1, :]), :, 1)  # Convert to column vector
    
    # Check if Y is a sparse matrix
    if issparse(Y)
        # If sparse, convert to dense matrix before calculation
        i0 = Array(Matrix(Y) * V0)
        return transpose(i0)
    else
        # If dense, calculate directly
        i0 = Array(Y * V0)
        return transpose(i0)
    end
end

"""
    I2_from_V012(V012, Y)

Calculate negative sequence current from sequence voltages and admittance matrix.
"""
function I2_from_V012(V012, Y)
    # Extract negative sequence component from symmetrical components
    V2 = reshape(transpose(V012[3, :]), :, 1)  # Convert to column vector
    
    # Check if Y is a sparse matrix
    if issparse(Y)
        # If sparse, convert to dense matrix before calculation
        i2 = Array(Matrix(Y) * V2)
        return transpose(i2)
    else
        # If dense, calculate directly
        i2 = Array(Y * V2)
        return transpose(i2)
    end
end

"""
    S_from_VI_elementwise(v, i)

Calculate complex power from voltage and current: S = V * conj(I).
"""
function S_from_VI_elementwise(v, i)
    return v .* conj(i)
end

"""
    V1_from_jpc(jpc)

Extract positive sequence voltage from JPC structure.
"""
function V1_from_jpc(jpc)
    return transpose(
        jpc["busAC"][:, VM] .* exp.(1im .* deg2rad.(jpc["busAC"][:, VA]))
    )
end

"""
    combine_X012(X0, X1, X2)

Combine zero, positive, and negative sequence components into one matrix.
"""
function combine_X012(X0, X1, X2)
    return vcat(X0, X1, X2)
end

"""
    current_from_voltage_results(y_0_pu, y_1_pu, v_012_pu)

Calculate sequence currents from sequence voltages and admittance matrices.
"""
function current_from_voltage_results(y_0_pu, y_1_pu, v_012_pu)
    I012_pu = combine_X012(I0_from_V012(v_012_pu, y_0_pu),
                          I1_from_V012(v_012_pu, y_1_pu),
                          I2_from_V012(v_012_pu, y_1_pu))
    return I012_pu
end

"""
    get_bus_v_results_3ph(jpc_3ph, jpc0, jpc1, jpc2)

Calculate and store three-phase voltage results for all buses.
"""
function get_bus_v_results_3ph(jpc_3ph, jpc0, jpc1, jpc2)

    V012_pu = zeros(ComplexF64, 3, length(jpc1["busAC"][:, BUS_I]))

    V012_pu[1, :] = jpc0["busAC"][:, VM] .* exp.(1im * pi/180 * jpc0["busAC"][:, VA])
    V012_pu[2, :] = jpc1["busAC"][:, VM] .* exp.(1im * pi/180 * jpc1["busAC"][:, VA])
    V012_pu[3, :] = jpc2["busAC"][:, VM] .* exp.(1im * pi/180 * jpc2["busAC"][:, VA])
    # Uncomment to get results in kV units rather than per unit
    # bus_base_kv = ppc0["bus"][:, BASE_KV] ./ sqrt(3)
    # V012_pu = V012_pu .* bus_base_kv
    
    Vabc_pu = sequence_to_phase(V012_pu)
    
    jpc_3ph["res_bus_3ph"] = Array{Float64}(undef, size(jpc1["busAC"], 1), 14)
    jpc_3ph["res_bus_3ph"][:,RES_3PH_BUS] = jpc_3ph["busAC_1"][:, BUS_I]
    # Calculate voltage magnitude
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VM_A] = abs.(Vabc_pu[1, :])
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VM_B] = abs.(Vabc_pu[2, :])
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VM_C] = abs.(Vabc_pu[3, :])
    
    # Voltage angle
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VA_A] = angle.(Vabc_pu[1, :]) .* 180 ./ π
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VA_B] = angle.(Vabc_pu[2, :]) .* 180 ./ π
    jpc_3ph["res_bus_3ph"][:,RES_3PH_VA_C] = angle.(Vabc_pu[3, :]) .* 180 ./ π
    
    jpc_3ph["res_bus_3ph"][:,RES_3PH_UNBALANCED] = abs.(V012_pu[3, :] ./ V012_pu[2, :]) .* 100
    
    return jpc_3ph
end

"""
    get_p_q_results_3ph(case, jpc_3ph::PowerFlow.JPC_3ph)

Calculate and aggregate active and reactive power results for three-phase system.
"""
function get_p_q_results_3ph(case, jpc_3ph::PowerFlow.JPC_3ph)

    # Create result matrix (bus, p in kw, q in kvar)
    bus_pq = zeros(Float64, length(jpc_3ph["busAC_1"][:,BUS_I]), 6)
    b, pA, pB, pC, qA, qB, qC = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

    # Todo: Voltage dependent loads
    elements = ["storages", "sgensAC", "loadsAC"]
    elements_3ph = ["loadsAC_asymm"]
    # elements_3ph = ["loadsAC_asymm", "asymmetric_sgen"]
    
    for element in elements
        # Check if element exists in jpc
        if  length(getproperty(case, Symbol(element)))==0
            continue
        end
        
        sign = element in ["sgensAC", "asymmetric_sgen"] ? -1 : 1
        if length(getproperty(case, Symbol(element))) > 0
            write_pq_results_to_element(case, jpc_3ph, element, "3ph")
            p_el, q_el, bus_el = get_p_q_b(case, jpc_3ph, element,"3ph")
            pA = vcat(pA, sign .* p_el ./ 3)
            pB = vcat(pB, sign .* p_el ./ 3)
            pC = vcat(pC, sign .* p_el ./ 3)
            qA = vcat(qA, sign .* q_el ./ 3 )
            qB = vcat(qB, sign .* q_el ./ 3 )
            qC = vcat(qC,  sign .* q_el ./ 3 )
            b = vcat(b, bus_el)
        end
    end
    
    for element in elements_3ph
        if length(getproperty(case, Symbol(element)))==0
            continue  # Skip this iteration if it doesn't exist
        end
        sign = element in ["sgensAC", "asymmetric_sgen"] ? -1 : 1
        if length(getproperty(case, Symbol(element))) > 0
            write_pq_results_to_element_3ph(net, element)
            p_el_A, q_el_A, p_el_B, q_el_B, p_el_C, q_el_C, bus_el = get_p_q_b_3ph(net, element)
            pA = vcat(pA, sign .* p_el_A)
            pB = vcat(pB, sign .* p_el_B)
            pC = vcat(pC, sign .* p_el_C)
            qA = vcat(qA, sign .* q_el_A )
            qB = vcat(qB, sign .* q_el_B )
            qC = vcat(qC, sign .* q_el_C )
            b = vcat(b, bus_el)
        end
    end

    # Summarize pq results for each element by group, to be written to jpc['bus'] later
    b_pp, vp_A, vq_A, vp_B, vq_B, vp_C, vq_C = _sum_by_group_nvals(Int64.(b), pA, qA, pB, qB, pC, qC)
    bus_pq[b_pp, 1] = vp_A
    bus_pq[b_pp, 2] = vq_A
    bus_pq[b_pp, 3] = vp_B
    bus_pq[b_pp, 4] = vq_B
    bus_pq[b_pp, 5] = vp_C
    bus_pq[b_pp, 6] = vq_C
    
    return bus_pq
end

"""
    write_pq_results_to_element(case, jpc_3ph, element, suffix=nothing)

Write power results to the specified element in jpc_3ph.
"""
function write_pq_results_to_element(case, jpc_3ph, element, suffix=nothing)

    # info element
    el_data = getproperty(case, Symbol(element))
    res_ = "res_" * element
    if suffix !== nothing
        res_ *= "_" * suffix
    end

    # Wards and xwards have different names in their element table, but not in res table. Also no scaling -> Fix...
    scaling = element in ["ward", "xward"] ? 1.0 : [load.scaling for load in el_data]

    element_in_service = [el_data.in_service for el_data in el_data]
    # P result in mw to element
    jpc_3ph["res_loadsAC_3ph"] = hcat([load.p_mw for load in el_data].* scaling .* element_in_service,[load.q_mvar for load in el_data].* scaling .* element_in_service)

    return jpc_3ph
end

"""
    get_p_q_b(case, jpc_3ph, element, suffix=nothing)

Get active and reactive power values and bus indices for a specific element.
"""
function get_p_q_b(case, jpc_3ph, element, suffix=nothing)
    res_ = "res_" * element
    if suffix !== nothing
        res_ *= "_$(suffix)"
    end

    el_data = getproperty(case, Symbol(element))
    # bus values are needed for stacking
    b = [el_data.bus for el_data in el_data]
    p = jpc_3ph["res_loadsAC_3ph"][:, 1]  # Assuming column 1 is active power
    q = jpc_3ph["res_loadsAC_3ph"][:, 2]  # Assuming column 2 is reactive power
    
    return p, q, b
end

"""
    write_pq_results_to_element_3ph(net, element)

Get p_mw and q_mvar for a specific pq element ("load", "sgen"...).
This function basically writes values from element table to res_element table.

Parameters:
    net: pandapower network
    element: element name (str)
Returns:
    net with updated results
"""
function write_pq_results_to_element_3ph(net, element)
    # info element
    el_data = net[element]
    res_ = "res_" * element * "_3ph"
    
    # Create results DataFrame if it doesn't exist
    if !haskey(net, res_)
        net[res_] = DataFrame()
    end

    scaling = el_data.scaling
    element_in_service = el_data.in_service .== true

    # Assumed column indices for p_mw and q_mvar
    p_mw_idx = columnindex(el_data, :p_a_mw) # Assumed column index for total active power
    q_mvar_idx = columnindex(el_data, :q_a_mvar) # Assumed column index for total reactive power

    # Create result arrays
    if element in ["load", "sgen"]
        # For load and sgen, distribute total power equally to three phases
        p_a_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
        p_b_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
        p_c_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
    else
        # For other elements, use specific power for each phase
        p_a_mw = el_data.p_a_mw .* scaling .* element_in_service
        p_b_mw = el_data.p_b_mw.* scaling .* element_in_service
        p_c_mw = el_data.p_c_mw .* scaling .* element_in_service
    end

    # Store results in net
    net[res_].p_a_mw = p_a_mw
    net[res_].p_b_mw = p_b_mw
    net[res_].p_c_mw = p_c_mw

    
    # Calculate reactive power results
    if element in ["load", "sgen"]
        q_a_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
        q_b_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
        q_c_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
    else
        q_a_mvar = el_data.q_a_mvar.* scaling .* element_in_service
        q_b_mvar = el_data.q_b_mvar .* scaling .* element_in_service
        q_c_mvar = el_data.q_c_mvar.* scaling .* element_in_service
    end
    
    # Store reactive power results in net
    net[res_].q_a_mvar = q_a_mvar
    net[res_].q_b_mvar = q_b_mvar
    net[res_].q_c_mvar = q_c_mvar

    return net
end

"""
    get_p_q_b_3ph(net, element)

Get three-phase active and reactive power values and bus indices for a specific element.
"""
function get_p_q_b_3ph(net, element)

    res_ = "res_" * element * "_3ph"

    # bus values are needed for stacking
    b = net[element].bus # In Julia, we don't need .values
    pA = net[res_].p_a_mw
    pB = net[res_].p_b_mw
    pC = net[res_].p_c_mw
    
    # Using ternary operator ? : in Julia
    qA = net[res_].q_a_mvar 
    qB = net[res_].q_b_mvar 
    qC = net[res_].q_c_mvar
    
    return pA, qA, pB, qB, pC, qC, b
end

"""
    _sum_by_group_nvals(bus, vals...)

Group values by bus and sum them within each group for multiple value arrays.
"""
function _sum_by_group_nvals(bus, vals...)
    # Get sorting order
    order = sortperm(bus)
    sorted_bus = bus[order]
    
    # Create index array marking the start position of each unique bus
    index = ones(Bool, length(sorted_bus))
    if length(sorted_bus) > 1
        index[2:end] = sorted_bus[2:end] .!= sorted_bus[1:end-1]
    end
    
    # Extract unique bus values
    unique_bus = sorted_bus[index]
    
    # Create result arrays
    newvals = ntuple(i -> zeros(eltype(vals[i]), length(unique_bus)), length(vals))
    
    # Group and sum each value
    for (i, val) in enumerate(vals)
        sorted_val = val[order]
        # Calculate cumulative sum
        cumulative_sum = cumsum(sorted_val)
        
        # Extract cumulative sum for each group
        group_sums = cumulative_sum[index]
        
        # Calculate sum for each group (by differencing)
        if length(group_sums) > 1
            group_sums[2:end] = group_sums[2:end] .- group_sums[1:end-1]
        end
        
        # Store result
        newvals[i] .= group_sums
    end
    
    # Return result tuple, first element is unique bus values
    return (unique_bus, newvals...)
end

"""
    get_branch_results_3ph(jpc_3ph, jpc0, jpc1, jpc2, pq_buses)

Extract branch results and write them to the appropriate dataframes.

Parameters:
    results: result of runpf loadflow calculation
    p: dict to store "res_line" and "res_trafo" Dataframes
"""
function get_branch_results_3ph(jpc_3ph, jpc0, jpc1, jpc2, pq_buses)
    I012_f, S012_f, V012_f, I012_t, S012_t, V012_t = get_branch_flows_3ph(jpc0, jpc1, jpc2)
    # get_line_results_3ph(net, jpc0, jpc1, jpc2, I012_f, V012_f, I012_t, V012_t)
    # get_trafo_results_3ph(net, jpc0, jpc1, jpc2, I012_f, V012_f, I012_t, V012_t)
    # _get_trafo3w_results(net, jpc, s_ft, i_ft)
    # _get_impedance_results(net, jpc, i_ft)
    # _get_xward_branch_results(net, jpc, bus_lookup_aranged, pq_buses)
    # _get_switch_results(net, i_ft)
end

"""
    get_branch_flows_3ph(jpc0, jpc1, jpc2)

Calculate branch flows for three-phase system.
"""
function get_branch_flows_3ph(jpc0, jpc1, jpc2)

    br_from_idx = Int64.(real(jpc1["branchAC"][:, F_BUS]))
    br_to_idx = Int64.(real(jpc1["branchAC"][:, T_BUS]))
    
    # Calculate from-end voltage - use vec instead of flatten
    V012_f = [vec(jpc["busAC"][br_from_idx, VM] .* jpc["busAC"][br_from_idx, BASE_KV] .*
                 exp.(1im .* deg2rad.(jpc["busAC"][br_from_idx, VA]))) for jpc in [jpc0, jpc1, jpc2]]
    V012_f = transpose(hcat(V012_f...))  # Convert to array shape equivalent to Python
    
    # Calculate to-end voltage
    V012_t = [vec(jpc["busAC"][br_to_idx, VM] .* jpc["busAC"][br_to_idx, BASE_KV] .*
                 exp.(1im .* deg2rad.(jpc["busAC"][br_to_idx, VA]))) for jpc in [jpc0, jpc1, jpc2]]
    V012_t = transpose(hcat(V012_t...))  # Convert to array shape equivalent to Python
    
    # Calculate from-end complex power
    S012_f = [real(jpc["branchAC"][:, PF]) .+ 1im .* real(jpc["branchAC"][:, QF]) for jpc in [jpc0, jpc1, jpc2]]
    S012_f = transpose(hcat(S012_f...))  # Convert to array shape equivalent to Python
    
    # Calculate to-end complex power
    S012_t = [real(jpc["branchAC"][:, PT]) .+ 1im .* real(jpc["branchAC"][:, QT]) for jpc in [jpc0, jpc1, jpc2]]
    S012_t = transpose(hcat(S012_t...))  # Convert to array shape equivalent to Python
    
    # Calculate current
    I012_f = I_from_SV_elementwise(S012_f, V012_f ./ sqrt(3))
    I012_t = I_from_SV_elementwise(S012_t, V012_t ./ sqrt(3))

    return I012_f, S012_f, V012_f, I012_t, S012_t, V012_t
end

"""
    I_from_SV_elementwise(S, V)

Calculate current from complex power and voltage: I = conj(S/V).
"""
function I_from_SV_elementwise(S, V)
    # Create zero array of same size as S
    result = zeros(Complex{Float64}, size(S))
    
    # For elements where V is non-zero, calculate conj(S/V)
    for i in eachindex(S)
        if V[i] != 0
            result[i] = conj(S[i] / V[i])
        end
    end
    
    return result
end

"""
    get_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2, pq_bus)

Calculate and store three-phase generator results.
"""
function get_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2, pq_bus)

    eg_end = size(case.ext_grids,1)
    if size(case.gensAC,1)>0
        gen_end = eg_end + length( [el_data.in_service for el_data in case.gensAC] )
    else
        gen_end = eg_end
    end
    
    # Get external grid results
    b, pA, qA, pB, qB, pC, qC = get_ext_grid_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2)

    # Get generator results
    if gen_end > eg_end
        b, pA, qA, pB, qB, pC, qC = get_pp_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2, b, pA, qA, pB, qB, pC, qC)
    end

    # Sum values by group
    b_pp, pA_sum, qA_sum, pB_sum, qB_sum, pC_sum, qC_sum = _sum_by_group_nvals(
        convert.(Int64, b), pA, qA, pB, qB, pC, qC)
    
    # Update bus power
    b_jpc = b_pp
    pq_bus[b_jpc, 1] .-= pA_sum
    pq_bus[b_jpc, 2] .-= qA_sum
    pq_bus[b_jpc, 3] .-= pB_sum
    pq_bus[b_jpc, 4] .-= qB_sum
    pq_bus[b_jpc, 5] .-= pC_sum
    pq_bus[b_jpc, 6] .-= qC_sum
end

"""
    get_ext_grid_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2)

Calculate and store three-phase external grid results.
"""
function get_ext_grid_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2)

    # Get external grid results
    n_res_eg = length(case.ext_grids)
    
    eg_idx_net = [el_data.index for el_data in case.ext_grids]
    # Fix: Use bus index of generator instead of generator index
    eg_bus_idx_net = [el_data.bus for el_data in case.ext_grids]
    
    # Read results for these buses from jpc
    V012 = zeros(ComplexF64, (3, n_res_eg))
    V012_temp= [net["busAC"][eg_bus_idx_net, VM] .* net["busAC"][eg_bus_idx_net, BASE_KV] .*
                          exp.(1im .* deg2rad.(net["busAC"][eg_bus_idx_net, VA]))
                          for net in [jpc0, jpc1, jpc2]]
    V012[:, eg_idx_net] = transpose(hcat(V012_temp...))  # Convert to array shape equivalent to Python
    S012 = zeros(ComplexF64, (3, n_res_eg))
    S012_temp = [(net["genAC"][eg_idx_net, PG] .+ 1im .* 
                           net["genAC"][eg_idx_net, QG])
                           for net in [jpc0, jpc1, jpc2]]
    S012[:, eg_idx_net] = transpose(hcat(S012_temp...))  # Convert to array shape equivalent to Python
    Sabc, Vabc = SVabc_from_SV012(S012, V012 ./ sqrt(3), n_res=n_res_eg, idx=eg_idx_net)

    pA = vec(real.(Sabc[1, :]))
    pB = vec(real.(Sabc[2, :]))
    pC = vec(real.(Sabc[3, :]))
    qA = vec(imag.(Sabc[1, :]))
    qB = vec(imag.(Sabc[2, :]))
    qC = vec(imag.(Sabc[3, :]))

    # Copy result indices
    jpc_3ph["res_ext_grid_3ph"].bus =  [el_data.bus for el_data in case.ext_grids]
    jpc_3ph["res_ext_grid_3ph"].p_a_mw = pA
    jpc_3ph["res_ext_grid_3ph"].p_b_mw = pB
    jpc_3ph["res_ext_grid_3ph"].p_c_mw = pC
    jpc_3ph["res_ext_grid_3ph"].q_a_mvar = qA
    jpc_3ph["res_ext_grid_3ph"].q_b_mvar = qB
    jpc_3ph["res_ext_grid_3ph"].q_c_mvar = qC

    # Get bus values for pq_bus
    b = [el_data.bus for el_data in case.ext_grids]
    
    return b, pA, qA, pB, qB, pC, qC
end

"""
    get_pp_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2, b, pA, qA, pB, qB, pC, qC)

Calculate and store three-phase generator results.
"""
function get_pp_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2, b, pA, qA, pB, qB, pC, qC)
    pA_gen, qA_gen, pB_gen, qB_gen, pC_gen, qC_gen = _get_p_q_gen_results_3ph(case, jpc_3ph, ppc0, ppc1, ppc2)
    _get_v_gen_results_3ph(net, ppc0, ppc1, ppc2)

    ac = net["_options"]["ac"]

    net["res_gen_3ph"].index = net["gen"].index
    b = vcat(b, net["gen"].bus)  # In Julia use vcat instead of np.hstack

    pA = vcat(pA, pA_gen)
    pB = vcat(pB, pB_gen)
    pC = vcat(pC, pC_gen)
    
    if ac
        qA = vcat(qA, qA_gen)
        qB = vcat(qB, qB_gen)
        qC = vcat(qC, qC_gen)
    end

    return b, pA, qA, pB, qB, pC, qC
end

"""
    _get_p_q_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2)

Calculate active and reactive power for generators in three-phase system.
"""
function _get_p_q_gen_results_3ph(case, jpc_3ph, jpc0, jpc1, jpc2)
    
    # indices of in service gens in the ppc
    # if any(gen_is_mask)
    #     gen_idx_ppc = gen_lookup[gen_is_idx]
    # else
    #     gen_idx_ppc = Int64[]
    # end

    # read results from ppc for these buses
    n_res_gen = length(case.gensAC)
    gen_idx_ppc = findall(jpc1["genAC"][:, GEN_STATUS] .== 1)  # Get indices of all in-use generators
    # 2 ext_grids Fix: Instead of the generator index, bus indices of the generators are used
    gen_bus_idx_ppc = Int64.(real(jpc1["genAC"][gen_idx_ppc, GEN_BUS]))

    V012 = zeros(ComplexF64, 3, n_res_gen)
    V012[:, gen_is_idx] = [ppc["busAC"][gen_bus_idx_ppc, VM] .* 
                          exp.(1im .* deg2rad.(ppc["busAC"][gen_bus_idx_ppc, VA])) 
                          for ppc in [jpc0, jpc1, jpc2]]

    S012 = zeros(ComplexF64, 3, n_res_gen)
    S012[:, gen_is_idx] = [-(ppc["gen"][gen_idx_ppc, PG] .+ 1im .* ppc["gen"][gen_idx_ppc, QG]) 
                          for ppc in [ppc0, ppc1, ppc2]]
                          
    I012 = zeros(ComplexF64, 3, n_res_gen)
    I012[:, gen_is_idx] = I_from_SV_elementwise(S012[:, gen_is_idx], V012[:, gen_is_idx])

    Vabc = sequence_to_phase(V012)
    Iabc = sequence_to_phase(I012)
    Sabc = S_from_VI_elementwise(Vabc, Iabc) .* 1e3
    
    pA, pB, pC = [vec(real(Sabc[i, :])) for i in 1:3]
    qA, qB, qC = [vec(imag(Sabc[i, :])) for i in 1:3]

    net["res_gen_3ph"]["p_a_mw"] = pA
    net["res_gen_3ph"]["p_b_mw"] = pB
    net["res_gen_3ph"]["p_c_mw"] = pC
    net["res_gen_3ph"]["q_a_mvar"] = qA
    net["res_gen_3ph"]["q_b_mvar"] = qB
    net["res_gen_3ph"]["q_c_mvar"] = qC

    return pA, qA, pB, qB, pC, qC
end

"""
    SVabc_from_SV012(S012, V012; n_res=nothing, idx=nothing)

Convert sequence components of power and voltage to phase components.
"""
function SVabc_from_SV012(S012, V012; n_res=nothing, idx=nothing)
    if isnothing(n_res)
        n_res = size(S012, 2)
    end
    
    if isnothing(idx)
        idx = trues(n_res)
    end
    
    I012 = zeros(ComplexF64, (3, n_res))
    I012[:, idx] = I_from_SV_elementwise(S012[:, idx], V012[:, idx])
    
    Vabc = sequence_to_phase(V012)
    Iabc = sequence_to_phase(I012)
    Sabc = S_from_VI_elementwise(Vabc, Iabc)
    
    return Sabc, Vabc
end

"""
    get_bus_results_3ph(case, jpc_3ph, bus_pq)

Store three-phase bus results in the appropriate fields.
"""
function get_bus_results_3ph(case, jpc_3ph, bus_pq)

    # update index in res bus bus
    # jpc["res_bus"].index = jpc["bus"].index
    jpc_3ph["res_bus_3ph"][:,RES_3PH_BUS] = jpc_3ph["busAC_1"][:, BUS_I]
    # write sum of p and q values to bus
    jpc_3ph["res_bus_3ph"][:,RES_3PH_PA_MW] = bus_pq[:, 1]  # Julia indices start at 1
    jpc_3ph["res_bus_3ph"][:,RES_3PH_PB_MW]= bus_pq[:, 3]  # Python's 2 corresponds to Julia's 3
    jpc_3ph["res_bus_3ph"][:,RES_3PH_PC_MW] = bus_pq[:, 5]  # Python's 4 corresponds to Julia's 5
    

    jpc_3ph["res_bus_3ph"][:,RES_3PH_QA_MVAR] = bus_pq[:, 2]  # Python's 1 corresponds to Julia's 2
    jpc_3ph["res_bus_3ph"][:,RES_3PH_QB_MVAR] = bus_pq[:, 4]  # Python's 3 corresponds to Julia's 4
    jpc_3ph["res_bus_3ph"][:,RES_3PH_QC_MVAR] = bus_pq[:, 6]  # Python's 5 corresponds to Julia's 6


    # Todo: OPF

    
    
    return nothing
end

"""
    _clean_up(net, res=true)

Clean up temporary data structures after power flow calculation.
"""
function _clean_up(net, res=true)
    
    # Most code is commented out, I kept these comments
    # mode = net.__internal_options["mode"]

    # set internal selected _is_elements to nothing. This way it is not stored (saves disk space)
    # net._is_elements = nothing

    #    mode = net["_options"]["mode"]
    #    if res
    #        res_bus = mode == "sc" ? net["res_bus_sc"] : 
    #                  mode == "pf_3ph" ? net["res_bus_3ph"] : 
    #                  net["res_bus"]
    #    end
    #    if length(net["trafo3w"]) > 0
    #        buses_3w = net["trafo3w"]["ad_bus"]
    #        deleteat!(net["bus"], buses_3w)
    #        select!(net["trafo3w"], Not(:ad_bus))
    #        if res
    #            deleteat!(res_bus, buses_3w)
    #        end
    #    end
    #
    #    if length(net["xward"]) > 0
    #        xward_buses = net["xward"]["ad_bus"]
    #        deleteat!(net["bus"], xward_buses)
    #        select!(net["xward"], Not(:ad_bus))
    #        if res
    #            deleteat!(res_bus, xward_buses)
    #        end
    #    end

    if length(net["dcline"]) > 0
        # Get generator indices corresponding to DC lines
        dc_gens = net.gen.index[(length(net.gen) - length(net.dcline) * 2):end]
        # Remove these generators
        net.gen = net.gen[setdiff(1:nrow(net.gen), dc_gens), :]
        if res
            net.res_gen = net.res_gen[setdiff(1:nrow(net.res_gen), dc_gens), :]
        end
    end
    
    return nothing
end

"""
    robust_process(net, jpc)

Process network and JPC data to ensure robust operation.
"""
function robust_process(net, jpc)

    bus = jpc["bus"]
    branch = jpc["branch"]
    gen = jpc["gen"]

    nb = size(bus, 1)
    
    internal_index = collect(1:nb)
    # Here internal_index is an array of integers from 1 to nb, representing indices of all buses
    # You can modify this index array as needed to select specific buses for processing
    bus_i = bus[:, BUS_I]

    # Create mapping from internal_index to bus_i
    idx2bus = Dict(i => bus_i[i] for i in internal_index)

    # Create mapping from bus_i to internal_index
    bus2idx = Dict(bus_i[i] => i for i in internal_index)

    net["bus2idx"] = bus2idx
    net["idx2bus"] = idx2bus
    # Update bus
    bus[:, BUS_I] = map(k -> bus2idx[k], bus[:, BUS_I])
    # Update generator and branch bus indices
    gen[:, GEN_BUS] = map(k -> bus2idx[k], gen[:, GEN_BUS])
    branch[:, F_BUS] = map(k -> bus2idx[k], branch[:, F_BUS])
    branch[:, T_BUS] = map(k -> bus2idx[k], branch[:, T_BUS])

    # Update external grid bus indices
    if haskey(net, "ext_grid")
        net["ext_grid"].bus = map(k -> bus2idx[k], net["ext_grid"].bus)
    end
    if haskey(net, "gen")
        net["gen"].bus = map(k -> bus2idx[k], net["gen"].bus)
    end
    if haskey(net, "load")
        net["load"].bus = map(k -> bus2idx[k], net["load"].bus)
    end
    if haskey(net, "sgen")
        net["sgen"].bus = map(k -> bus2idx[k], net["sgen"].bus)
    end
    if haskey(net, "storage")
        net["storage"].bus = map(k -> bus2idx[k], net["storage"].bus)
    end
    if haskey(net, "shunt")
        net["shunt"].bus = map(k -> bus2idx[k], net["shunt"].bus)
    end
    if haskey(net, "trafo")
        net["trafo"].hv_bus = map(k -> bus2idx[k], net["trafo"].hv_bus)
        net["trafo"].lv_bus = map(k -> bus2idx[k], net["trafo"].lv_bus)
    end
    if haskey(net, "trafo3w")
        net["trafo3w"].hv_bus = map(k -> bus2idx[k], net["trafo3w"].hv_bus)
        net["trafo3w"].mv_bus = map(k -> bus2idx[k], net["trafo3w"].mv_bus)
        net["trafo3w"].lv_bus = map(k -> bus2idx[k], net["trafo3w"].lv_bus)
    end
    if haskey(net, "line")
        net["line"].hv_bus = map(k -> bus2idx[k], net["line"].from_bus)
        net["line"].lv_bus = map(k -> bus2idx[k], net["line"].from_bus)
    end

    mask = branch[:, BR_STATUS] .== 1
    branch = branch[mask, :]
    jpc["bus"] = bus
    jpc["branch"] = branch
    jpc["gen"] = gen
end

"""
    get_full_branch_zero(jpc_3ph)

Add non-service branches to zero sequence branch data.
"""
function get_full_branch_zero(jpc_3ph)

    branch_not_in_service = jpc_3ph.branchAC_0[findall(jpc_3ph.branchAC_0[:, BR_STATUS] .== 0),:]
    branch0 = jpc_3ph.branchAC_0

    lb = size(branch_not_in_service, 1)
    connected_branches = zeros(Int64, lb, 4)
    branch_not_in_service = hcat(branch_not_in_service, connected_branches)
    if size(branch_not_in_service,1)>0
        branch0 = vcat(branch0, branch_not_in_service)
    end

    jpc_3ph["branchAC_0"] = branch0
end
