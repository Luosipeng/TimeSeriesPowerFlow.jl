
function loadflow_result_ETAP(file_path,IEEE_sheet_path)
    # file_path = "C:/ETAP 2200/NEWTEST/Untitled.LF1S - Load Flow Report.xlsx"
    sheets_data = Dict{String, DataFrame}()

    # 打开 Excel 文件并读取工作表
    XLSX.openxlsx(file_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            data = sheet[:]
            sheets_data[sheet_name] = DataFrame(data, :auto)
        end
    end



    # 提取数据
    result_reort = sheets_data["Sheet1"]
    bus_ID=result_reort[:,1]
    bus_voltage=result_reort[:,2]
    bus_angle=result_reort[:,3]


    # Create a compared dictory
    Compared_dict_voltage=Dict(bus_ID.=>bus_voltage)
    Compared_dict_angle=Dict(bus_ID.=>bus_angle)

    #Compare the IEEE forrmat result and ETAP format result
    # IEEE_sheet_path = "C:/Users/13733/Desktop/DistributionPowerFlow-luosipeng/output.xlsx"
    sheets_data = Dict{String, DataFrame}()

    # 打开 Excel 文件并读取工作表
    XLSX.openxlsx(IEEE_sheet_path) do wb
        for sheet_name in XLSX.sheetnames(wb)
            sheet = XLSX.getsheet(wb, sheet_name)
            data = sheet[:]
            sheets_data[sheet_name] = DataFrame(data, :auto)
        end
    end

    # 提取数据
    IEEE_result_reort = sheets_data["Sheet1"]
    IEEE_busID=IEEE_result_reort[2:end,1]
    IEEE_bus_voltage=IEEE_result_reort[2:end,2]
    IEEE_bus_angle=IEEE_result_reort[2:end,3]

    #Calculation
    misvoltage=IEEE_bus_voltage.-map(k -> Compared_dict_voltage[k],IEEE_busID)
    misangle=IEEE_bus_angle.-map(k -> Compared_dict_angle[k],IEEE_busID)
    misv=misvoltage./(map(k -> Compared_dict_voltage[k],IEEE_busID))
    misa=misangle./(map(k -> Compared_dict_angle[k],IEEE_busID))
    # mis=abs.(misvoltage)
    # mis=abs.(misvoltage)
    
    # average_mis=sum(mis)/length(mis)
    return misv,misa
end