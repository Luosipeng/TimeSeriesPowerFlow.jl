"""
Process PV AC systems and integrate them into the power system model.
This function handles AC-side photovoltaic systems by adding them to the generation data.
"""
function process_pv_acsystem(pv_acsystem, jpc)
    # Get AC-side PV systems that are in service
    mask = pv_acsystem[:, PV_AC_IN_SERVICE] .== 1
    pv_acsystem = pv_acsystem[mask, :]

    if isempty(pv_acsystem)
        # If there are no AC-side PV systems, return the original jpc
        return jpc
    end

    # Create a new genAC array to store AC-side PV system data
    genAC = zeros(size(pv_acsystem, 1), 26)
    genAC[:, 1] = pv_acsystem[:, PV_AC_BUS]  # bus
    genAC[:, 2] = pv_acsystem[:, PV_AC_INVERTER_PAC]   # Pg
    genAC[:, 3] = pv_acsystem[:, PV_AC_INVERTER_QAC]   # Qg
    genAC[:, 4] .= pv_acsystem[:, PV_AC_INVERTER_QAC_MAX] # Qg_max
    genAC[:, 5] .= pv_acsystem[:, PV_AC_INVERTER_QAC_MIN] # Qg_min
    genAC[:, 6] .= 1.0  # Voltage setpoint
    genAC[:, 7] .= 100.0  # Base MVA
    genAC[:, 8] .= 1.0  # Status (on)
    genAC[:, 9] .= 9999.0  # Maximum active power

    # Set bus type to PV for voltage control mode inverters
    for i in eachindex(pv_acsystem[:,PV_AC_INVERTER_MODE])
       if pv_acsystem[i,PV_AC_INVERTER_MODE] == 1
            bus_index = Int(findfirst(jpc["busAC"][:, BUS_I] .== pv_acsystem[i, PV_AC_BUS]))
            jpc["busAC"][bus_index, BUS_TYPE] = PV
       end
    end

    # Append the new PV generators to the existing genAC data
    jpc["genAC"] = vcat(jpc["genAC"], genAC)
    return jpc
end
