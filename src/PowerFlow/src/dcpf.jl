"""
    dcpf(B, Pbus, Va0, ref, pv, pq)

Solve a DC power flow problem.

# Arguments
- `B`: Full system B matrix (susceptance matrix)
- `Pbus`: Vector of bus real power injections
- `Va0`: Initial vector of bus voltage angles in radians
- `ref`: Index of reference (slack) bus
- `pv`: Vector of PV bus indices
- `pq`: Vector of PQ bus indices

# Returns
- `Va`: Vector of bus voltage angles in radians
- `success`: Boolean indicating whether the solution was successful

# Description
This function solves the DC power flow equations to determine bus voltage angles.
The DC power flow is a linearized approximation of the AC power flow where:
1. Voltage magnitudes are assumed to be 1.0 per unit
2. Angle differences are assumed to be small
3. Line resistances are neglected

The function solves the linear system of equations:
    B[pv,pq],[pv,pq] * Va[pv,pq] = Pbus[pv,pq] - B[pv,pq],ref * Va0[ref]

A solution is considered unsuccessful if any voltage angle exceeds an arbitrary
threshold value (currently set to 1e5), which typically indicates numerical issues
or an unsolvable system.
"""
function dcpf(B, Pbus, Va0, ref, pv, pq)
    # DCPF  Solves a DC power flow.
    # Solves for the bus voltage angles at all but the reference bus, given the full system
    # B matrix and the vector of bus real power injections, the initial
    # vector of bus voltage angles (in radians), and column vectors with
    # the lists of bus indices for the swing bus, PV buses, and PQ buses,
    # respectively. Returns a vector of bus voltage angles in radians.

    # constant
    Va_threshold = 1e5     # arbitrary threshold on |Va| for declaring failure

    # initialize result vector
    Va = Va0
    success = true    # successful by default

    # update angles for non-reference buses
    Va[[pv; pq]] = B[[pv; pq], [pv; pq]] \ (Pbus[[pv; pq]] - B[[pv; pq], ref] * Va0[ref])

    # check for presence of *any* warning
    if maximum(abs.(Va)) > Va_threshold
        success = false
    end

    return Va, success
end
