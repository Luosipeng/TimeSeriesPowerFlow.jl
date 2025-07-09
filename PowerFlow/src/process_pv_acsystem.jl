function process_pv_acsystem(pv_acsystem, jpc)
    # 获取服务中的交流侧光伏系统
    mask = pv_acsystem[:, PV_AC_IN_SERVICE] .== 1
    pv_acsystem = pv_acsystem[mask, :]

    if isempty(pv_acsystem)
        # 如果没有交流侧光伏系统，返回原始的jpc
        return jpc
    end

    # 创建一个新的genAC数组来存储交流侧光伏系统数据
    genAC = zeros(size(pv_acsystem, 1), 26)
    genAC[:, 1] = pv_acsystem[:, PV_AC_BUS]  # bus
    genAC[:, 2] = pv_acsystem[:, PV_AC_INVERTER_PAC]   # Pg
    genAC[:, 3] = pv_acsystem[:, PV_AC_INVERTER_QAC]   #
    genAC[:, 4] .= pv_acsystem[:, PV_AC_INVERTER_QAC_MAX] # Qg_max
    genAC[:, 5] .= pv_acsystem[:, PV_AC_INVERTER_QAC_MIN] # Qg_min
    genAC[:, 6] .= 1.0
    genAC[:, 7] .= 100.0
    genAC[:, 8] .= 1.0
    genAC[:, 9] .= 9999.0

    for i in eachindex(pv_acsystem[:,PV_AC_INVERTER_MODE])
       if pv_acsystem[i,PV_AC_INVERTER_MODE] == 1
            bus_index = Int(findfirst(jpc["busAC"][:, BUS_I] .== pv_acsystem[i, PV_AC_BUS]))
            jpc["busAC"][bus_index, BUS_TYPE] = PV
       end
    end

    jpc["genAC"] = vcat(jpc["genAC"], genAC)
    return jpc
end