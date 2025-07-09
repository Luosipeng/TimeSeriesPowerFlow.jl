function acdccompare(busAC,dict_bus,dict_new)
    #find the original buses ID
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,
        VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =  PowerFlow.idx_bus();
    bus=busAC
    reversed_dict_new = Dict(value => key for (key, value) in dict_new)
    reversed_dict = Dict(value => key for (key, value) in dict_bus)

    bus[:,BUS_I]=map(k -> reversed_dict_new[k],bus[:,BUS_I])
    bus_ID=map(k -> reversed_dict[k],bus[:,BUS_I])

    #Output an Excel file
    data = [bus_ID bus[:,8]]       # 将两个向量组合成矩阵

    # 写入 Excel 文件
    # 合并数据并写入到 Excel
    XLSX.openxlsx("output.xlsx", mode="w") do workbook
        # 创建一个新的工作表
        sheet = workbook["Sheet1"]  # 如果工作表不存在，会自动创建
        # 写入表头
        sheet[1, 1] = "ID"
        sheet[1, 2] = "Voltage"
        sheet[1, 3] = "Angle"
    for i in eachindex(bus_ID)
        sheet[i + 1, 1] = bus_ID[i]  # 第一列写入名字
        sheet[i + 1, 2] = bus[i, VM]  # 第二列写入电压
        sheet[i + 1, 3] = bus[i, VA]  # 第三列写入角度
    end
    end

    println("Excel file 'output.xlsx' created successfully!")
    
end