"""
    makeYbus(baseMVA, bus, branch)

Build the bus admittance matrix and branch admittance matrices.

# Arguments
- `baseMVA`: Base MVA for the system
- `bus`: Bus data matrix with columns representing bus parameters
- `branch`: Branch data matrix with columns representing branch parameters

# Returns
- `Ybus`: Bus admittance matrix
- `Yf`: Branch admittance matrix for "from" end of branches
- `Yt`: Branch admittance matrix for "to" end of branches

# Description
This function builds the bus admittance matrix (Ybus) and branch admittance matrices (Yf and Yt)
for a given power system network. The admittance matrices are essential components for power flow
and other power system analyses.

The function:
1. Computes branch series admittances and line charging susceptances
2. Handles tap ratios and phase shifters
3. Incorporates bus shunt admittances
4. Builds connection matrices between branches and buses
5. Constructs the complete bus admittance matrix

# Notes
- All admittance values are in per-unit on system MVA base
- Requires buses to be numbered consecutively (internal ordering)
- Branch status values determine which branches are in service
"""
function makeYbus(baseMVA, bus, branch)
    ## define named indices into bus, branch matrices
    # constants
    nb = size(bus, 1)          # number of buses
    nl = size(branch, 1)       # number of lines
    # check that bus numbers are equal to indices to bus (one set of bus numbers)

    if any(bus[:, BUS_I] .!= (1:nb))
        error("makeYbus: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering")
    end
    # for each branch, compute the elements of the branch admittance matrix where
    Stat = branch[:, BR_STATUS]                    # ones at in-service branches
    Ys = Stat ./ (branch[:, BR_R] .+ 1im * branch[:, BR_X])  # series admittance
    Bc = Stat .* branch[:, BR_B]                           # line charging susceptance
    tap = ones(nl,1)                              # default tap ratio = 1
    i = findall(branch[:, TAP] .!= 0)                       # indices of non-zero tap ratios
    tap[i] = branch[i, TAP]                        # assign non-zero tap ratios
    tap = tap .* exp.(1im*pi/180 * branch[:, SHIFT]) # add phase shifters
    Ytt = Ys + 1im*Bc/2
    Yff = Ytt ./ (tap .* conj.(tap))
    Yft = - Ys ./ conj.(tap)
    Ytf = - Ys ./ tap
     # compute shunt admittance
     Ysh = (bus[:, GS] .+ 1im * bus[:, BS]) / baseMVA  # vector of shunt admittances
     # bus indices
    f = branch[:, F_BUS]                           # list of "from" buses
    t = branch[:, T_BUS]                           # list of "to" buses
    ## build connection matrices
    Cf = sparse(1:nl, f, vec(ones(nl, 1)), nl, nb);      ## connection matrix for line & from buses
    Ct = sparse(1:nl, t, vec(ones(nl, 1)), nl, nb);      ## connection matrix for line & to buses
    ## build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    ## at each branch's "from" bus, and Yt is the same for the "to" bus end
    Yf = sparse(1:nl, 1:nl, vec(Yff), nl, nl) * Cf + sparse(1:nl, 1:nl, vec(Yft), nl, nl) * Ct
    Yt = sparse(1:nl, 1:nl, vec(Ytf), nl, nl) * Cf + sparse(1:nl, 1:nl, vec(Ytt), nl, nl) * Ct
    ## build Ybus
    Ybus = Cf' * Yf + Ct' * Yt +             ## branch admittances
        sparse(1:nb, 1:nb, Ysh, nb, nb);    ## shunt admittance
    # end

    return Ybus, Yf, Yt
end

"""
    branch_vectors(branch, nl)

Compute branch admittance vectors for power flow calculations.

# Arguments
- `branch`: Branch data matrix with columns representing branch parameters
- `nl`: Number of branches (lines)

# Returns
- `Ytt`: Self-admittance at "to" bus end
- `Yff`: Self-admittance at "from" bus end
- `Yft`: Mutual admittance from "from" bus to "to" bus
- `Ytf`: Mutual admittance from "to" bus to "from" bus

# Description
This function computes the four components of branch admittance matrices needed for
power flow calculations:
1. Self-admittances at each end of the branch (Ytt, Yff)
2. Mutual admittances between the ends of the branch (Yft, Ytf)

The function handles:
- Branch status (in-service or out-of-service)
- Series impedance conversion to admittance
- Line charging susceptance
- Transformer tap ratios and phase shifters

# Notes
- Currently commented out in the code
- Used as a helper function for building branch admittance matrices
- All admittance values are in per-unit
"""
# function  branch_vectors(branch, nl)
#     (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
#     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
#     (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
#     ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
    
#     Stat = branch[:, BR_STATUS]                    # ones at in-service branches
#     Ysf = Stat ./ (branch[:, BR_R] .+ 1im * branch[:, BR_X])  # series admittance
#     Yst = Ysf
#     Bcf = Stat.* branch[:, BR_B] 
#     Bct = Bcf
#     tap = ones(nl,1)                              # default tap ratio = 1
#     i = findall(branch[:, TAP] .!= 0)                       # indices of non-zero tap ratios
#     tap[i] = real.(branch[i, TAP])                        # assign non-zero tap ratios
#     tap = tap .* exp.((1im * pi / 180) .* branch[:, SHIFT]) # add phase shifters
#     Ytt = Yst + Bct / 2
#     Yft = (Ysf+ Bcf/2) ./ (tap .* conj.(tap))
#     Yft = - Ysf ./ conj(tap)
#     Ytf = - Yst ./ tap

#     return Ytt, Yff, Yft, Ytf
# end
