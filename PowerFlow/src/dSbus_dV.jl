function dSbus_dV(Ybus, V, vcart=0)
    n = length(V)
    Ibus = Ybus * V

    diagV = Diagonal(V)
    diagIbus = Diagonal(Ibus)

    if vcart == 0
        diagVnorm = Diagonal(V ./ abs.(V))
    end

    if vcart == 1
        dSbus_dV1 = conj(diagIbus) + diagV * conj(Ybus)  # dSbus/dVr
        dSbus_dV2 = 1im * (conj(diagIbus) - diagV * conj(Ybus))  # dSbus/dVi
    else
        dSbus_dV1 = 1im * diagV * conj(diagIbus - Ybus * diagV)  # dSbus/dVa
        dSbus_dV2 = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm  # dSbus/dVm
    end

    return dSbus_dV1, dSbus_dV2
end