"""
    dcpfsoln(baseMVA::Float64, bus0::Matrix{Float64}, gen0::Matrix{Float64}, 
             branch0::Matrix{Float64}, load0::Matrix{Float64}, Ybus, Yf, Yt, V, ref, p)

Update power system state variables after a DC power flow solution.

# Arguments
- `baseMVA::Float64`: Base MVA value for the system
- `bus0::Matrix{Float64}`: Initial bus data matrix
- `gen0::Matrix{Float64}`: Initial generator data matrix
- `branch0::Matrix{Float64}`: Initial branch data matrix
- `load0::Matrix{Float64}`: Load data matrix
- `Ybus`: Bus admittance matrix
- `Yf`: From-bus branch admittance matrix
- `Yt`: To-bus branch admittance matrix
- `V`: Complex bus voltage vector solution
- `ref`: Reference (slack) bus indices
- `p`: P bus indices

# Returns
- `bus`: Updated bus data matrix
- `gen`: Updated generator data matrix
- `branch`: Updated branch data matrix with power flows

# Description
This function updates the power system state variables after a DC power flow solution
has been obtained. It performs the following operations:

1. Updates bus voltage magnitudes
2. Updates generator reactive power (Qg) for generators at PV and slack buses
   - Distributes reactive power proportionally among multiple generators at the same bus
   - Respects generator reactive power limits
3. Updates active power (Pg) for generators at slack buses
4. Calculates branch power flows
5. Expands the branch matrix if needed to store power flow results

The function handles special cases like multiple generators at the same bus and
generators with identical reactive power limits.
"""
function dcpfsoln(baseMVA::Float64, bus0::Matrix{Float64}, gen0::Matrix{Float64}, branch0::Matrix{Float64},load0::Matrix{Float64}, Ybus, Yf, Yt, V, ref, p)
    
    # Initialize return values
    bus = bus0
    gen = gen0
    branch = branch0
    load = load0
    ##----- update bus voltages -----
    bus[:, VM] = V

    ##----- update Qg for gens at PV/slack buses and Pg for slack bus(es) -----
    ## generator info
    on=findall((gen[:, GEN_STATUS].>0).*(bus[Int.(gen[:, GEN_BUS]), BUS_TYPE] .!= P))   # Which generators are on and not at PQ buses?
    off = findall(x -> x <= 0, gen[:, GEN_STATUS])  # Which generators are off?
    gbus = Int.(gen[on, GEN_BUS])  # What buses are they at?
    # Compute total injected bus powers   
    Sbus = V[gbus] .* conj.(Ybus[gbus, :] * V)
    # Update Qg for generators at PV/slack buses
    gen[off, QG] .= 0  # Zero out off-line Qg

    Ld=zeros(size(bus,1),size(load,2))
    Ld[:,1]=collect(1:size(bus,1))
    Ld[:,2]=collect(1:size(bus,1))
    Ld[:,3].=1
    Ld[Int.(load[:,LOAD_CND]),2:end] = load[:,2:end]

    Pd_gbus, Qd_gbus = total_load(bus[gbus, :],Ld[gbus,:]);
    gen[on, QG] = imag(Sbus) * baseMVA + Qd_gbus;   ## inj Q + local Qd
    if length(on) > 1
        # build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
        nb = size(bus, 1)
        ngon = size(on, 1)
        Cg = sparse((1:ngon), gbus, ones(ngon), ngon, nb)
        # divide Qg by number of generators at the bus to distribute equally
        ngg = Cg * sum(Cg, dims=1)'    # ngon x 1, number of gens at this gen's bus
        gen[on, QG] = gen[on, QG] ./ ngg
        #precision control
        #gen = round.(gen, digits=6)
        ## set finite proxy M for infinite limits (for ~ proportional splitting)
        ## equal to sum over all gens at bus of abs(Qg) plus any finite Q limits
        Qmin = gen[on, QMIN]
        Qmax = gen[on, QMAX]
        M = abs.(gen[on, QG])
        M[.!isinf.(Qmax)] = M[.!isinf.(Qmax)] .+ abs.(Qmax[.!isinf.(Qmax)])
        M[.!isinf.(Qmin)] = M[.!isinf.(Qmin)] .+ abs.(Qmin[.!isinf.(Qmin)])
        M = Cg * Cg' * M   # each gen gets sum over all gens at same bus
        #precision control
        #M = round.(M, digits=6)
        # replace +/- Inf limits with proxy +/- M
        Qmin[Qmin .==  Inf] =  M[Qmin .==  Inf]
        Qmin[Qmin .== -Inf] = -M[Qmin .== -Inf]
        Qmax[Qmax .==  Inf] =  M[Qmax .==  Inf]
        Qmax[Qmax .== -Inf] = -M[Qmax .== -Inf]
        # divide proportionally
        Cmin = sparse(vec(1:ngon), gbus, Qmin, ngon, nb)
        Cmax = sparse(vec(1:ngon), gbus, Qmax, ngon, nb);
        Qg_tot = Cg' * gen[on, QG]     # nb x 1 vector of total Qg at each bus
        Qg_min = sparse(sum(Cmin,dims=1)')            # nb x 1 vector of min total Qg at each bus
        Qg_max = sparse(sum(Cmax,dims=1)')            # nb x 1 vector of max total Qg at each bus
        eps = 2.2204e-16
        gen[on, QG] = Qmin .+ (Cg * ((Qg_tot .- Qg_min) ./ (Qg_max .- Qg_min .+ eps))) .* (Qmax .- Qmin)          # avoid div by 0
        # fix gens at buses with Qg range = 0 (use equal violation for all)
        # To do for the correction
        try
            ig = findall(abs.(Cg * (Qg_min .- Qg_max)) .< 10*eps)  # gens at buses with Qg range = 0
            if !isempty(ig)
                ig_int = [Tuple(i)[1] for i in ig]
                ib = findall(sum(Cg[ig_int, :], dims=1)'[:].>0)   # buses with Qg range = 0
                # total mismatch @ bus div by number of gens
                mis = sparse(ib, ones(Int64, length(ib)), vec(Array((Qg_tot[ib] .- Qg_min[ib]) ./ sum(Cg[:, ib]', dims=2)[:])), nb, 1)
                gen[on[ig], QG] = Qmin[ig] .+ Cg[getindex.(ig, 1), :] * mis
            end
        catch
            println("Error in the correction of the Qg range = 0")
        end
    end
    
     for k in eachindex(ref)
        refgen = findall(gbus .== ref[k])              # which is(are) the reference gen(s)?
        busm = bus[ref[k], :]
        busm = reshape(busm, :,1)
        busm = busm'
        loadm=Ld[ref[k],:]
        loadm = reshape(loadm, :,1)
        loadm = loadm'
        Pd_refk, Qd_refk = total_load(busm,loadm)
         gen[on[refgen[1]], PG] = real(Sbus[refgen[1]]) * baseMVA + Pd_refk[1, 1]  # inj P + local Pd
         if length(refgen) > 1       # more than one generator at this ref bus
             # subtract off what is generated by other gens at this bus
             gen[on[refgen[1]], PG] = gen[on[refgen[1]], PG] - sum(gen[on[refgen[2:end]], PG])
         end
         #precision control
         #gen = round.(gen, digits=6)
     end
    # update/compute branch power flows
    out = findall(branch[:, BR_STATUS] .== 0)      # out-of-service branches
    br = findall(branch[:, BR_STATUS].==1)            # in-service branches
    Sf = V[Int.(branch[br, F_BUS])] .* conj.(Yf[br, :] * V) * baseMVA  # complex power at "from" bus
    St = V[Int.(branch[br, T_BUS])] .* conj.(Yt[br, :] * V) * baseMVA  # complex power injected at "to" bus
    #branch = [branch zeros(size(branch, 1), 17-size(branch, 2))]
    # Determine the current size of the branch matrix
    (rows, cols) = size(branch)

    # Determine the number of columns to add
    cols_to_add = 18 - cols
    # If cols_to_add is greater than 0, add more columns
    if cols_to_add > 0
        # Add columns filled with zeros
        branch = [branch zeros(rows, cols_to_add)]
    end
    branch[br, [PF, QF, PT, QT]] = [real.(Sf) imag.(Sf) real.(St) imag.(St)]
    branch[out, [PF, QF, PT, QT]] = zeros(length(out), 4)
    #precision control
    #branch = round.(branch, digits=6)
    return bus, gen, branch
end
