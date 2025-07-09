function makeBdc(baseMVA, bus, branch)
    # constants
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
     ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = idx_brch()
    
    nb = size(bus, 1)          # number of buses
    nl = size(branch, 1)       # number of lines

    # check that bus numbers are equal to indices to bus (one set of bus numbers)
    if any(bus[:, BUS_I] .!= collect(1:nb))
        error("makeBdc: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering")
    end

    # for each branch, compute the elements of the branch B matrix and the phase
    # shift "quiescent" injections
    stat = branch[:, BR_STATUS]                    # ones at in-service branches
    b = stat ./ branch[:, BR_X]                    # series susceptance
    tap = ones(nl)                                 # default tap ratio = 1
    i = findall(branch[:, TAP] .!= 0)              # indices of non-zero tap ratios
    tap[i] = branch[i, TAP]                        # assign non-zero tap ratios
    b = b ./ tap

    # build connection matrix Cft = Cf - Ct for line and from - to buses
    f = branch[:, F_BUS]                           # list of "from" buses
    t = branch[:, T_BUS]                           # list of "to" buses
    i = vcat(1:nl, 1:nl)                           # double set of row indices
    Cft = sparse(i, vcat(f, t), vcat(ones(nl), -ones(nl)), nl, nb)    # connection matrix

    # build Bf such that Bf * Va is the vector of real branch powers injected
    # at each branch's "from" bus
    Bf = sparse(i, vcat(f, t), vcat(b, -b), nl, nb)    # = spdiagm(0 => b) * Cft

    # build Bbus
    Bbus = Cft' * Bf

    # build phase shift injection vectors
    Pfinj = b .* (-branch[:, SHIFT] * pi/180)      # injected at the from bus ...
    # Ptinj = -Pfinj                               # ... and extracted at the to bus
    Pbusinj = Cft' * Pfinj                         # Pbusinj = Cf * Pfinj + Ct * Ptinj

    return Bbus, Bf, Pbusinj, Pfinj
end