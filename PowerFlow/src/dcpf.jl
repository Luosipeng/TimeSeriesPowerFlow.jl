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