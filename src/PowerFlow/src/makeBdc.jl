"""
    makeBdc(baseMVA, bus, branch)

Build the DC power flow matrices Bbus and Bf from bus and branch data.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `branch`: Branch data matrix with columns representing branch parameters

# Returns
- `Bbus`: Nodal susceptance matrix (n×n sparse matrix)
- `Bf`: Branch susceptance matrix (nl×n sparse matrix) such that Bf*Va gives branch power injections at "from" buses
- `Pbusinj`: Vector of bus real power injections due to phase shifters
- `Pfinj`: Vector of branch real power injections at "from" buses due to phase shifters

# Description
This function builds the matrices necessary for DC power flow analysis:
- Bbus is the nodal susceptance matrix, where Bbus*Va = P gives the bus power injections
- Bf is the branch susceptance matrix, where Bf*Va gives the branch power injections at "from" buses
- Pbusinj is the vector of bus power injections due to phase shifters
- Pfinj is the vector of branch power injections at "from" buses due to phase shifters

The DC power flow model assumes:
1. Flat voltage profile (all voltage magnitudes = 1.0 p.u.)
2. Small angle differences between connected buses
3. Negligible branch resistance compared to reactance

# Notes
- Buses must be numbered consecutively starting from 1
- Branch status, tap ratios, and phase shifters are taken into account
- Use ext2int() to convert to internal bus numbering if needed

# Constants Used (assumed to be defined elsewhere)
- BUS_I: Column index for bus number in bus matrix
- BR_STATUS: Column index for branch status in branch matrix
- BR_X: Column index for branch reactance in branch matrix
- TAP: Column index for tap ratio in branch matrix
- F_BUS: Column index for "from" bus number in branch matrix
- T_BUS: Column index for "to" bus number in branch matrix
- SHIFT: Column index for phase shift in branch matrix
"""
function makeBdc(baseMVA, bus, branch)
    
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
    