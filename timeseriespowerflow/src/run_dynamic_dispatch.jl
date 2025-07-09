using  JuMP
using  Ipopt
using  Gurobi

function run_dynamic_dispatch(new_jpc,
        Cld_ac, Cld_dc, 
        loadAC_PD, loadAC_QD,
        loadDC_PD,
        genAC_PG,
        Cgen_ac, 
        Cconv_ac, Cconv_dc, 
        η_rec, η_inv,
        Cpv_ac,
        Cpv_dc,
        pv_ac_p_mw_ratio,
        pv_ac_p_mw,
        pv_max_p_mw,
        pv_max_p_mw_ratio,
        Cstorage_ac,
        ess_initial_soc,
        ess_max_soc,
        ess_min_soc,
        ess_power_capacity_mw,
        ess_energy_capacity_mwh,
        ess_efficiency,
        day_price_line,
        num_hours=24)

    # Extract connection information
    branchAC = new_jpc.branchAC
    branchDC = new_jpc.branchDC
    converter = new_jpc.converter

    rac = branchAC[:,BR_R]
    rdc = branchDC[:,BR_R]
    rconv = zeros(size(converter, 1))
    # Extract bus count
    n_nodes = size(new_jpc.busAC, 1) + size(new_jpc.busDC, 1)
    A , branch_data = TimeSeriesPowerFlow.build_incidence_matrix_td(n_nodes, branchAC, branchDC, converter)

     # Initialize
    nbr = size(A, 1)  # Number of branches
    nc = size(converter, 1)  # Number of converters
    ng = size(genAC_PG, 1)  # Number of AC generators
    npv = size(new_jpc.pv, 1)  # Number of PV devices
    ns = size(new_jpc.storage, 1)  # Number of storage devices

    # Search for branch indices
    # dcbranch = filter(x -> x[3] == 2, branch_data)
    dcbranch_indices = findall(x -> x[3] == 2, branch_data)
    converter_indices = findall(x -> x[3] == 3, branch_data)
    acbranch_indices = findall(x -> x[3] == 1, branch_data)

    r = zeros(nbr)  # Initialize branch resistances
    r[acbranch_indices] = rac  # AC branch resistances
    r[dcbranch_indices] = rdc  # DC branch resistances
    r[converter_indices] = rconv  # Converter resistances

    # Optimization model configuration
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)  # Output verbosity
    set_optimizer_attribute(model, "max_iter", 3000)  # Maximum iterations
    set_optimizer_attribute(model, "tol", 1e-6)       # Convergence tolerance

    # model = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(model, "OutputFlag", 0)  # Show solving process
    # set_optimizer_attribute(model, "TimeLimit", 3600)  # Maximum solving time (seconds)
    # set_optimizer_attribute(model, "MIPGap", 1e-4)    # MIP relative gap


    # Define variables
    @variable(model, Pij[1:nbr,1:num_hours])  # Branch power flow
    @variable(model, Pgen[1:ng,1:num_hours])  # Inverter power flow
    @variable(model, Pij_inv[1:nc,1:num_hours] >= 0)  # Converter power flow
    @variable(model, Pij_rec[1:nc,1:num_hours] >= 0)  # Inverter power flow
    @variable(model, P_pv_mw[1:npv,1:num_hours] >= 0)  # PV power flow
    @variable(model, soc[1:ns,1:num_hours])  # Storage power flow
    @variable(model, ess_charge[1:ns,1:num_hours])  # Storage power flow
    @variable(model, ess_discharge[1:ns,1:num_hours])  # Storage power flow
    @variable(model, 0 <= ess_mode[1:ns, 1:num_hours] <= 1)  # Relaxed binary
    # @variable(model, ess_mode[1:ns, 1:num_hours], Bin)  # Strict binary variable

    # PV power constraints
    for t in 1:num_hours
        @constraint(model, 0 .<= P_pv_mw[:,t] .<= pv_max_p_mw[:].*pv_max_p_mw_ratio[:,t])  # PV power injection lower limit
    end

    # AC generator power constraints
    # for i in 1:ng
    #     @constraint(model, Pgen[i, :] >= 0)  # Generator power injection lower limit
    # end

    # Converter power constraints
    # for t in 1:num_hours
    #     @constraint(model, Pij_inv[:,t].*Pij_rec[:,t].==0)  # Inverter power injection lower limit
    # end

    # Converter active power constraints
    for t in 1:num_hours
         @constraint(model, Pij[converter_indices,t] .== 0)
    end
   

    # Storage power constraints
    for e in 1:ns
        for t in 1:num_hours
            # Charge/discharge power limits
            @constraint(model, ess_charge[e, t] >= 0)
            @constraint(model, ess_discharge[e, t] >= 0)
            @constraint(model, ess_charge[e, t] <= ess_power_capacity_mw[e])
            @constraint(model, ess_discharge[e, t] <= ess_power_capacity_mw[e])
            
            # Charge/discharge mutual exclusivity constraints
            @constraint(model, ess_discharge[e, t] <= ess_power_capacity_mw[e] * (1 - ess_mode[e, t]))
            @constraint(model, ess_charge[e, t] <= ess_power_capacity_mw[e] * ess_mode[e, t])
        end
        
        # Initial SOC constraint
        @constraint(model, soc[e, 1] == ess_initial_soc[e] - 
                        ess_discharge[e, 1] / ess_efficiency[e] / ess_energy_capacity_mwh[e] + 
                        ess_charge[e, 1] * ess_efficiency[e] / ess_energy_capacity_mwh[e])
        
        # SOC evolution over time
        for t in 2:num_hours
            @constraint(model, soc[e, t] == soc[e, t-1] - 
                            ess_discharge[e, t] / ess_efficiency[e] / ess_energy_capacity_mwh[e] + 
                            ess_charge[e, t] * ess_efficiency[e] / ess_energy_capacity_mwh[e])
        end
        
        # SOC upper and lower limits
        for t in 1:num_hours
            @constraint(model, soc[e, t] >= ess_min_soc[e])
            @constraint(model, soc[e, t] <= ess_max_soc[e])
        end
        
        # Cyclic constraint: Final SOC equals initial SOC
        @constraint(model, soc[e, num_hours] == ess_initial_soc[e])
    end

   # Power balance constraints
    for t in 1:num_hours
        @constraint(model,A' * Pij[:,t] + Cld_ac*loadAC_PD[:,t] + Cld_dc*loadDC_PD[:,t] -Cgen_ac*Pgen[:,t] -Cpv_ac*(pv_ac_p_mw.*pv_ac_p_mw_ratio[:,t]) + Cconv_ac*(Pij_inv[:,t]-η_rec.*Pij_rec[:,t]) + Cconv_dc*(Pij_rec[:,t]-η_inv.*Pij_inv[:,t]) - Cpv_dc * P_pv_mw[:,t] + Cstorage_ac*(ess_charge[:,t]-ess_discharge[:,t]) .== 0)
    end

    original_objective = sum(day_price_line[t,2] * Pgen[g, t] for g in 1:ng for t in 1:num_hours)
    penalty_weight = 1000  # Penalty weight, needs adjustment
    # Objective function: Minimize total power flow
    # @objective(model, Min, sum(day_price_line[t,2] * Pgen[g, t] for g in 1:ng for t in 1:num_hours))
    if nc > 0 && num_hours > 0
        complementarity_penalty = sum(Pij_inv[i,t] * Pij_rec[i,t] for i in 1:nc for t in 1:num_hours)
        @objective(model, Min, original_objective + penalty_weight * complementarity_penalty)
    else
        @objective(model, Min, original_objective)
        println("Warning: nc or num_hours is zero, skipping complementarity penalty term")
    end
    # Objective function: Minimize the sum of squares of AC branch power flows
    # @objective(model, Min, sum(Pij[i,t]^2 *r[i]  for i in acbranch_indices for t in 1:num_hours)  )
    

    # Solve the model
    optimize!(model)

    # Check solution status
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED
        if status == MOI.OPTIMAL
            println("Optimization successful, optimal solution found!")
        elseif status == MOI.LOCALLY_SOLVED
            println("Optimization successful, locally optimal solution found!")
        elseif status == MOI.ALMOST_LOCALLY_SOLVED
            println("Optimization successful, approximately optimal solution found!")
        end
        
        # # Output objective function value
        obj_value = objective_value(model)
        # println("Objective function value: ", obj_value)
        
        # # Output key variable values
        # println("\nGenerator output (Pgen):")
        Pgen_values = value.(Pgen)
        # display(Pgen_values)

        # println("\nConverter power (Pij_inv):")
        Pij_inv_values = value.(Pij_inv)
        # display(Pij_inv_values)

        # println("\nConverter power (Pij_rec):")
        Pij_rec_values = value.(Pij_rec)
        # display(Pij_rec_values)
        
        # println("\nPV output (P_pv_mw):")
        P_pv_values = value.(P_pv_mw)
        # display(P_pv_values)
        
        # println("\nStorage charging power (ess_charge):")
        ess_charge_values = value.(ess_charge)
        # display(ess_charge_values)
        
        # println("\nStorage discharging power (ess_discharge):")
        ess_discharge_values = value.(ess_discharge)
        # display(ess_discharge_values)
        
        # println("\nStorage SOC (soc):")
        soc_values = value.(soc)
        # display(soc_values)
        
        # Return results
        return Dict(
            "status" => status,
            "objective" => obj_value,
            "Pgen" => Pgen_values,
            "Pij_inv" => Pij_inv_values,
            "Pij_rec" => Pij_rec_values,
            "P_pv_mw" => P_pv_values,
            "ess_charge" => ess_charge_values,
            "ess_discharge" => ess_discharge_values,
            "soc" => soc_values
        )
    else
        println("Optimization not successfully completed. Status: ", status)
        println("Solver message: ", raw_status(model))
        
        # Return error status
        return Dict("status" => status, "error" => raw_status(model))
    end
end

function build_incidence_matrix_td(n_nodes, branchAC, branchDC, converter)
    # Calculate total number of branches
    n_branches = size(branchAC, 1) + size(branchDC, 1) + size(converter, 1)
    
    # Extract start and end nodes for all branches, and mark branch type
    # Type markers: 1=AC, 2=DC, 3=Converter
    branch_data = Vector{Tuple{Int, Int, Int, Int}}(undef, n_branches)
    
    # Add AC branches
    for i in 1:size(branchAC, 1)
        # Ensure smaller node number as starting node
        node1 = Int(branchAC[i, 1])
        node2 = Int(branchAC[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[i] = (from, to, 1, i)  # (start node, end node, type, original index)
    end
    
    # Add DC branches
    offset = size(branchAC, 1)
    for i in 1:size(branchDC, 1)
        # Ensure smaller node number as starting node
        node1 = Int(branchDC[i, 1])
        node2 = Int(branchDC[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[offset + i] = (from, to, 2, i)
    end
    
    # Add converters
    offset = size(branchAC, 1) + size(branchDC, 1)
    for i in 1:size(converter, 1)
        # Ensure smaller node number as starting node
        node1 = Int(converter[i, 1])
        node2 = Int(converter[i, 2])
        from = min(node1, node2)
        to = max(node1, node2)
        branch_data[offset + i] = (from, to, 3, i)
    end
    
    # Sort by starting node
    sort!(branch_data, by = x -> (x[1], x[2]))
    
    # Create incidence matrix
    A = zeros(n_branches, n_nodes)
    
    # Build incidence matrix based on sorted branches
    for (idx, (from, to, type, _)) in enumerate(branch_data)
        A[idx, from] = 1   # Outflow from node is positive
        A[idx, to] = -1    # Inflow to node is negative
    end
    
    return A, branch_data
end
