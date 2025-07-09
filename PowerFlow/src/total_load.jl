"""
    total_load(bus, load)

Calculate the total fixed and dispatchable load at each bus in the system, considering ZIP load model.

# Arguments
- `bus`: Bus data matrix containing bus information
- `load`: Load data matrix containing load model percentages

# Returns
- `Pd`: Vector of real power demand at each bus
- `Qd`: Vector of reactive power demand at each bus (if want_Q=1)

This function computes the total load at each bus using the ZIP load model, which represents
loads as a combination of constant impedance (Z), constant current (I), and constant power (P)
components. The percentages for each component are specified in the load matrix.
"""
function total_load(bus, load)
    # define constants

    nb = size(bus, 1) # number of buses

    # default options
    want_Q = 1

    # fixed load at each bus, & initialize dispatchable
        pw_1=load[:,LOADP_PERCENT]
        pw_2=load[:,LOADI_PERCENT]
        pw_3=load[:,LOADZ_PERCENT]
        Sd = makeSdzip(1, bus, pw_1, pw_2, pw_3)
        Vm = bus[:, VM]
        Sbusd = Sd.p + Sd.i .* Vm + Sd.z .* Vm.^2
        Pdf = real(Sbusd) # real power
        if want_Q==1
            Qdf = imag(Sbusd) # reactive power
        end

    # dispatchable load at each bus
    
        Pdd = zeros(nb, 1)
        if want_Q==1
            Qdd = zeros(nb, 1)
        end

    # compute load sums
    
        Pd = (Pdf + Pdd) .* (bus[:, BUS_TYPE] .!= NONE)
        if want_Q==1
            Qd = (Qdf + Qdd) .* (bus[:, BUS_TYPE] .!= NONE)
        end
    return Pd, Qd
    end

"""
    total_load(bus)

Calculate the total fixed and dispatchable load at each bus in the system using default ZIP parameters.

# Arguments
- `bus`: Bus data matrix containing bus information

# Returns
- `Pd`: Vector of real power demand at each bus
- `Qd`: Vector of reactive power demand at each bus (if want_Q=1)

This is an overloaded version of `total_load` that uses default ZIP model parameters.
It computes the total load at each bus using the standard ZIP load model with default
percentages for constant impedance (Z), constant current (I), and constant power (P) components.
"""
function total_load(bus)
    nb = size(bus, 1) # number of buses

    # default options
    want_Q = 1

    # fixed load at each bus, & initialize dispatchable
        Sd = makeSdzip(1, bus)
        Vm = bus[:, VM]
        Sbusd = Sd.p + Sd.i .* Vm + Sd.z .* Vm.^2
        Pdf = real(Sbusd) # real power
        if want_Q==1
            Qdf = imag(Sbusd) # reactive power
        end

    # dispatchable load at each bus
    
        Pdd = zeros(nb, 1)
        if want_Q==1
            Qdd = zeros(nb, 1)
        end

    # compute load sums
    
        Pd = (Pdf + Pdd) .* (bus[:, BUS_TYPE] .!= NONE)
        if want_Q==1
            Qd = (Qdf + Qdd) .* (bus[:, BUS_TYPE] .!= NONE)
        end
    return Pd, Qd
    end
