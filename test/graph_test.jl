using Graphs  # 或 LightGraphs
using XLSX
using DataFrames
#This function is used to test the graphs package
# file_path = "C:/Users/13733/Desktop/NANRUI.xlsx"
    ##读取整个表格内容
    # 创建一个字典，用于存储每个工作表的 DataFrame
    file_path = "C:/Users/13733/Desktop/ONELINE_DATA.xlsx"
    sheets_data = Dict{String, DataFrame}()

    # 打开 Excel 文件并处理每个工作表
    XLSX.openxlsx(file_path) do wb
        for sheet_name in XLSX.sheetnames(wb)  # 遍历所有工作表名称
            sheet = XLSX.getsheet(wb, sheet_name)  # 获取当前工作表
            data = sheet[:]  # 读取整个表格内容
            sheets_data[sheet_name] = DataFrame(data, :auto)  # 存储为 DataFrame
        end
    end
    #用户定义
    #TODO:
    baseMVA=100.0
    baseKV=10
    ##分离文件
    #节点母线索引
    ##===================================idx=========================##
    EquipmentID=1;#序号
    Voltage=2;#电压
    Initial_Voltage=3;#潮流初始值（标幺值）
    In_Service=4;#是否投入工作
    State=5;
    Phase=7;#相数
    Wire_connection=8;#连接线
    Min_DF=9;#最小需求指数
    Max_DF=10;#最大需求指数
    Bus_Type=11;#节点类型，7：slack节点；0：PQ节点；2：PV节点
    #TODO:
    Cont_Rating=13;#连续额定容量
    Bracing_Symm=14;#母线对称短路电流的有效值
    Bracing_Asymm=15;# 非对称短路条件下的有效电流

    #节点母线数据
    countrow=0;
    for i in 11:size(sheets_data["Bus"],1)
        if(ismissing(sheets_data["Bus"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Bus"],2)
        if(ismissing(sheets_data["Bus"][10,i]))
            countcol=i-1
            break
        end
    end
    bus_data=sheets_data["Bus"][11:countrow,1:countcol]

    #发电机索引
    ##===================================idx=========================##
    GenID=1;#序号
    Gen_connected_element=2;#电压
    Gen_inservice=3;#潮流初始值（标幺值）
    Gen_state=4;#是否投入工作
    Gen_controlmode=5;
    Gen_power_rating=7;#相数
    Gen_power_rating_unit=8;#连接线
    Gen_apparent_power_rating=9;#最小需求指数
    Gen_apparent_power_rating_unit=10;#最大需求指数
    Gen_voltage=11;#连续额定容量
    Gen_PF=12;#母线对称短路电流的有效值
    Gen_EFF=13;# 非对称短路条件下的有效电流
    Gen_pole=14;
    Gen_Xd11=15;
    Gen_Xd11_Ra=16;
    Gen_X0=17;
    Gen_X0_R0=18;
    Gen_X2=19;
    Gen_X2_R2=20;
    Gen_Rdc=21;
    Gen_modeltype=22;
    Gen_Xd=23;
    Gen_Xq=24;
    Gen_Tdo1=25;
    Gen_Xdu=26;
    Gen_Xqu=27;

    #发电机数据
    countrow=0;
    for i in 11:size(sheets_data["Synchronous Generator"],1)
        if(ismissing(sheets_data["Synchronous Generator"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Synchronous Generator"],2)
        if(ismissing(sheets_data["Synchronous Generator"][10,i]))
            countcol=i-1
            break
        end
    end
    gen_data=sheets_data["Synchronous Generator"][11:countrow,1:countcol]

    #线缆索引
    ##===================================idx=========================##
    Cable_equipmentID=1;#序号
    Cable_nophase=2;#电压
    Cable_inservice=3;#潮流初始值（标幺值）
    Cable_state=4;#是否投入工作
    Cable_phase=6;
    Cable_Felement=8;#相数
    Cable_Telement=9;#连接线
    Cable_Freq=11;#最小需求指数
    Cable_length=24;#最大需求指数
    Cable_length_unit=25;#节点类型，7：slack节点；0：PQ节点；2：PV节点
    Cable_Tol=26;#连续额定容量
    Cable_Tmin=27;
    Cable_Tmax=28
    Cable_r1=29;# 
    Cable_x1=30;
    Cable_y1=31;
    Cable_r0=32;
    Cable_x0=33;
    Cable_y0=34;
    Cable_Type=35;
    Cable_perunit_length_value=36;
    Cable_perunit_length_unit=37;

    #线缆数据
    countrow=0;
    for i in 11:size(sheets_data["Cable"],1)
        if(ismissing(sheets_data["Cable"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Cable"],2)
        if(ismissing(sheets_data["Cable"][10,i]))
            countcol=i-1
            break
        end
    end
    cable_data=sheets_data["Cable"][11:countrow,1:countcol]

    #传输线索引
    ##===================================idx=========================##
    Line_equipmentID=1;#序号
    Line_length=2;#电压
    Line_lengthunit=3;
    Line_Tol=4;
    Line_inservice=5;
    Line_state=6;#是否投入工作
    Line_phase=7;
    Line_Felement=9;#相数
    Line_Telement=10;#连接线
    Line_baseT1=36;
    Line_baseT2=37;
    Line_Tmin=38;#母线对称短路电流的有效值
    Line_Tmax=39
    Line_r11=40;# 非对称短路条件下的有效电流
    Line_r12=41;
    Line_x1=42;
    Line_y1=43;
    Line_r21=44;# 非对称短路条件下的有效电流
    Line_r22=45;
    Line_x2=46;
    Line_y2=47;
    Line_r01=48;
    Line_r02=49;
    Line_x0=50;
    Line_y0=51;
    Line_Type=52;
    Line_perunit_length_value=53;
    Line_perunit_length_unit=54;

    #传输线数据
    countrow=0;
    for i in 11:size(sheets_data["Transmission Line"],1)
        if(ismissing(sheets_data["Transmission Line"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Transmission Line"],2)
        if(ismissing(sheets_data["Transmission Line"][10,i]))
            countcol=i-1
            break
        end
    end
    transline_data=sheets_data["Transmission Line"][11:countrow,1:countcol]

    #变压器索引
    ##===================================idx=========================##
    Trans_equipmentID=1;#序号
    Trans_inservice=2;#电压
    Trans_state=3;
    Trans_Pelement=5;
    Trans_Selement=6;
    Trans_phase=7;#是否投入工作
    Trans_centertap=8;
    Trans_coretype=9;#相数
    Trans_temprise=10;#连接线
    prisecrated_Power=11;
    prisecrated_Power_unit=12;
    Trans_rating=13;
    KVA_T1=14;
    KVA_T2=15;
    KVA_T3=16;
    Trans_Pvoltage=17;
    Trans_Pvoltage_unit=18;
    Trans_Pfla=19;
    Trans_Svoltage=20;
    Trans_Svoltage_unit=21;
    Trans_Sfla=22;
    Trans_Pos_value=23;
    Trans_Pos_unit=24;
    Trans_Pos_XRrating=25;
    Trans_zero_value=26;
    Trans_zero_unit=27;
    Trans_zero_XRrating=28;
    Trans_ZT1=29;
    Trans_XRT1=30;
    Trans_ZT2=31;
    Trans_XRT2=32;
    Trans_ZT3=33;
    Trans_XRT3=34;
    Trans_Ptap=35;
    Trans_Stap=36;

    #变压器
    countrow=0;
    for i in 11:size(sheets_data["Two-Winding Transformer"],1)
        if(ismissing(sheets_data["Two-Winding Transformer"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Two-Winding Transformer"],2)
        if(ismissing(sheets_data["Two-Winding Transformer"][10,i]))
            countcol=i-1
            break
        end
    end
    transformer_data=sheets_data["Two-Winding Transformer"][11:countrow,1:countcol]

    #负荷索引
    ##===================================idx=========================##
    MotorID=1;#序号
    ConectedID=2;#所在节点序号
    datatype=3;#数据类型
    priority=4;#优先级
    configstatus=5;
    load_inservice=6;#相数
    load_state=7;#连接线
    load_voltage=9;#最小需求指数
    load_kva=10;#最大需求指数
    load_unit=11;#节点类型，7：slack节点；0：PQ节点；2：PV节点
    load_pf=12;#连续额定容量
    load_phase=13;#母线对称短路电流的有效值
    load_contdf=14;# 非对称短路条件下的有效电流
    load_intermdf=15;
    load_grounding=17;
    motor_load=18;
    LRC=19;
    X_R=20;

    #集成负荷
    countrow=0;
    for i in 11:size(sheets_data["Lumped Load"],1)
        if(ismissing(sheets_data["Lumped Load"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Lumped Load"],2)
        if(ismissing(sheets_data["Lumped Load"][10,i]))
            countcol=i-1
            break
        end
    end
    Load_data=sheets_data["Lumped Load"][11:countrow,1:countcol]

     #电网索引
    ##===================================idx=========================##
    Utility_ID=1;
    Utility_connected_ID=2;
    Utility_Inservice=3;
    Utility_state=4;
    Utility_rated_voltage=6;
    Utility_control_mode=7;

    countrow=0;
    for i in 11:size(sheets_data["Utility"],1)
        if(ismissing(sheets_data["Utility"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Utility"],2)
        if(ismissing(sheets_data["Utility"][10,i]))
            countcol=i-1
            break
        end
    end
    Utility_data=sheets_data["Utility"][11:countrow,1:countcol]

     #逆变器索引
    ##===================================idx=========================##
    Inverter_ID=1;
    Inverter_DCConElement=2;
    Inverter_ACConElement=3;
    Inverter_DCbus=4;
    Inverter_ACbus=5;
    Inverter_inservice=6;
    Inverter_state=7;
    Inverter_operation_mode=11;
    Inverter_Output_connection=12;
    Inverter_DC_powerrating=13;
    Inverter_DC_Power_Unit=14;
    Inverter_DC_voltagev_rate=15;
    Inverter_DC_Vmax=16;
    Inverter_DC_Vmin=17;
    Inverter_AC_powerrating=18;
    Inverter_AC_Power_Unit=19;
    Inverter_AC_voltagekv_rate=20;
    Inverter_AC_Power_factor=21;
    Inverter_AC_grounding=22;

    countrow=0;
    for i in 11:size(sheets_data["Inverter"],1)
        if(ismissing(sheets_data["Inverter"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["Inverter"],2)
        if(ismissing(sheets_data["Inverter"][10,i]))
            countcol=i-1
            break
        end
    end
    Inverter_data=sheets_data["Inverter"][11:countrow,1:countcol]

     #HVCB索引
    ##===================================idx=========================##
    HVCB_ID=1;
    HVCB_FROM_ELEMENT=2;
    HVCB_TO_ELEMENT=3;
    HVCB_INSERVICE=4;
    HVCB_STATE=5;
    HVCB_STATUS=6;
    HVCB_DESCRIPTION=7;

    countrow=0;
    for i in 11:size(sheets_data["HVCB"],1)
        if(ismissing(sheets_data["HVCB"][i,1]))
            countrow=i-1
            break
        end
        countrow=i
    end
    countcol=0;
    for i in 1:size(sheets_data["HVCB"],2)
        if(ismissing(sheets_data["HVCB"][10,i]))
            countcol=i-1
            break
        end
    end
    HVCB_data=sheets_data["HVCB"][11:countrow,1:countcol]
##===========================IEEE Format Bus Data====================
    #Define indexing for IEEE bus data
    BUS_I=1;#BUS ID
    TYPE=2;
    PD=3;
    QD=4;
    GS=5;
    BS=6;
    AREA=7;
    VM=8;
    VA=9;
    BASEKV=10;
    ZONE=11;
    VMAX=12;
    VMIN=13;

    #Transfer data
    inservice=findall(bus_data[:,In_Service].=="Yes")#Detect the in_Serviced bus
    bus=zeros(length(inservice),13)#Initialize the Bus data
    bus[:,BUS_I]=collect(1:length(inservice))#Assign a new ID for all buses

    PQindex=findall(bus_data[:,Bus_Type].==0)#Searching all PQ buses
    bus[PQindex,TYPE].=1#define PQ buses in IEEE format

    Load_inservice=findall(Load_data[:,load_inservice].=="Yes")#Find the in-serviced load
    load_data=Load_data[Load_inservice,:]

    #Define a ID index 
    bus_index=collect(1:size(bus_data[inservice,:],1))
    dict_bus=Dict(zip(bus_data[inservice,EquipmentID],bus_index))

    #Find slack buses
    geninservice=findall(gen_data[:,Gen_inservice].=="Yes")
    gen_data=gen_data[geninservice,:]

    pv_index=findall(gen_data[:,Gen_controlmode].=="Voltage Control")
    pv_bus=map(k -> dict_bus[k],gen_data[pv_index,Gen_connected_element])
    bus[pv_bus,TYPE].=2
    pv_index=findall(Inverter_data[:,Inverter_operation_mode].=="Voltage Control")
    pv_bus=map(k -> dict_bus[k],Inverter_data[pv_index,Inverter_ACbus])
    bus[pv_bus,TYPE].=2



    # slack_index=findall(gen_data[:,Gen_controlmode].=="Swing")
    # slack_bus=map(k -> dict_bus[k],gen_data[slack_index,Gen_connected_element])
    slack_index=findall(Utility_data[:,Utility_control_mode].=="Swing")
    slack_bus=map(k -> dict_bus[k],Utility_data[slack_index,Utility_connected_ID])
    bus[slack_bus,TYPE].=3

    #Assign load data to bus 
    load_connected_bus=map(k -> dict_bus[k],load_data[:,ConectedID])#Searching the connected bus ID
    bus[load_connected_bus,PD]=0.001.*load_data[:,load_kva].*parse.(Float64,load_data[:,load_pf])./100#Assign the active power load for buses
    bus[load_connected_bus,QD]=0.001.*Load_data[:,load_kva].*sin.(acos.(parse.(Float64,load_data[:,load_pf])./100))#Assign the reactive power load for buses

    #TODO:Define the area property

    #Assign voltage to bus 
    # bus[:,VM]=bus_data[inservice,Voltage].*parse.(Float64,bus_data[inservice,Initial_Voltage])./100
    bus[:,VM].=1.0
    bus[:,VA].=0

    #BASEKV
    bus[:,BASEKV].=bus_data[inservice,Voltage].*parse.(Float64,bus_data[inservice,Initial_Voltage])./100
    #Assign zone information
    bus[:,ZONE].=1

    #Assign the Vmax and Vmin 
    bus[:,VMAX].=1.05
    bus[:,VMIN].=0.8

    ##===========================IEEE Format Branch Data====================
    #Define indexing for IEEE branch data
    FBUS=1;
    TBUS=2;
    R=3;
    X=4;
    B=5;
    RATEA=6;
    RATEB=7;
    RATEC=8;
    RATIO=9;
    ANGLE=10;
    STATUS=11;
    ANGMIN=12;
    ANGMAX=13;

    #Transfer data
    cableinservice=findall(cable_data[:,Cable_inservice].=="Yes")#Find all the in-serviced cable
    translineinservice=findall(transline_data[:,Line_inservice].=="Yes")#Find all the in-serviced transline
    transforminservice=findall(transformer_data[:,Trans_inservice].=="Yes")#Find all the in-serviced transform

    #Select the in-serviced branch
    cable=cable_data[cableinservice,:]
    transline=transline_data[translineinservice,:]
    transformer=transformer_data[transforminservice,:]

    #Delete the extra cable branch
    cable = filter(row -> row[Cable_Felement] !== missing, cable)
    cable = filter(row -> row[Cable_Telement] !== missing, cable)

    #Delete the extra transline branch
    transline = filter(row -> row[Line_Felement] !== missing, transline)
    transline = filter(row -> row[Line_Telement] !== missing, transline)

    #Delete the extra transform branch
    transformer = filter(row -> row[Trans_Pelement] !== missing, transformer)
    transformer = filter(row -> row[Trans_Selement] !== missing, transformer)

    #Initialize the matrix
    branch_cable=zeros(size(cable,1),13)
    branch_transline=zeros(size(transline,1),13)
    branch_transformer=zeros(size(transformer,1),13)

    #Assign the ID
    branch_cable[:,FBUS]=map(k -> dict_bus[k],cable[:,Cable_Felement])
    branch_cable[:,TBUS]=map(k -> dict_bus[k],cable[:,Cable_Telement])

    branch_transline[:,FBUS]=map(k -> dict_bus[k],transline[:,Line_Felement])
    branch_transline[:,TBUS]=map(k -> dict_bus[k],transline[:,Line_Telement])

    branch_transformer[:,FBUS]=map(k -> dict_bus[k],transformer[:,Trans_Pelement])
    branch_transformer[:,TBUS]=map(k -> dict_bus[k],transformer[:,Trans_Selement])

    #Assign the R\X\B to the branch_cable 
    #TODO:
    branch_cable[:,R]=cable[:,Cable_length].*cable[:,Cable_r1]./cable[:,Cable_perunit_length_value]
    branch_cable[:,X]=cable[:,Cable_length].*cable[:,Cable_x1]./cable[:,Cable_perunit_length_value]
    branch_cable[:,B]=cable[:,Cable_length].*cable[:,Cable_y1]./cable[:,Cable_perunit_length_value]

    #define the rateA rateB and rateC
    branch_cable[:,RATEA].=100#MVA
    branch_cable[:,RATEB].=100#MVA
    branch_cable[:,RATEC].=100#MVA

    branch_cable[:,RATIO].=0
    branch_cable[:,ANGLE].=0
    branch_cable[:,STATUS].=1
    branch_cable[:,ANGMIN].=-180
    branch_cable[:,ANGMAX].=180


    #Assign the R\X\B to the branch_transline
    #TODO:
    branch_transline[:,R]=transline[:,Line_length].*parse.(Float64,transline[:,Line_r11]).*0.000621371
    branch_transline[:,X]=transline[:,Line_length].*parse.(Float64,transline[:,Line_x1]).*0.000621371
    branch_transline[:,B]=transline[:,Line_length].*parse.(Float64,transline[:,Line_y1]).*0.000621371

    #define the rateA rateB and rateC
    branch_transline[:,RATEA].=100#MVA
    branch_transline[:,RATEB].=100#MVA
    branch_transline[:,RATEC].=100#MVA

    branch_transline[:,RATIO].=0
    branch_transline[:,ANGLE].=0
    branch_transline[:,STATUS].=1
    #TODO:
    branch_transline[:,ANGMIN].=-180
    branch_transline[:,ANGMAX].=180

    #Assign the R\X\B to the branch_transformer 
    branch_transformer[:,R]=0.01.*transformer[:,Trans_Pos_value].*transformer[:,Trans_Pvoltage].^2 .*inv.(transformer[:,prisecrated_Power]).*inv.(sqrt.(1 .+transformer[:,Trans_Pos_XRrating].^2))#based on the primary side voltage
    branch_transformer[:,X]=0.01.*transformer[:,Trans_Pos_value].*transformer[:,Trans_Pvoltage].^2 .*inv.(transformer[:,prisecrated_Power]).*inv.(sqrt.(1 .+transformer[:,Trans_Pos_XRrating].^2)).*transformer[:,Trans_Pos_XRrating]#based on the primary side voltage
    branch_transformer[:,B].=0

    #define the rateA rateB and rateC
    branch_transformer[:,RATEA].=100#MVA
    branch_transformer[:,RATEB].=100#MVA
    branch_transformer[:,RATEC].=100#MVA

    branch_transformer[:,RATIO].=0
    branch_transformer[:,ANGLE].=0
    branch_transformer[:,STATUS].=1
    branch_transformer[:,ANGMIN].=-180
    branch_transformer[:,ANGMAX].=180

    #The branch impedance is normalized to per unit.
    #For cables
    cable_basekv=bus[Int.(branch_cable[:,FBUS]),BASEKV]
    n=findall(cable_basekv[:].==10.0)
    cable_basekv[n]=cable_basekv[n].*1.05
    branch_cable[:,R]=branch_cable[:,R].*baseMVA.*inv.(cable_basekv.^2)
    branch_cable[:,X]=branch_cable[:,X].*baseMVA.*inv.(cable_basekv.^2)
    branch_cable[:,B]=branch_cable[:,B].*baseMVA.*inv.(cable_basekv.^2)

    #For transmission lines
    trans_basekv=bus[Int.(branch_transline[:,FBUS]),BASEKV]
    n=findall(trans_basekv[:].==10.0)
    trans_basekv[n]=trans_basekv[n].*1.05
    branch_transline[:,R]=branch_transline[:,R].*baseMVA.*inv.(trans_basekv.^2)
    branch_transline[:,X]=branch_transline[:,X].*baseMVA.*inv.(trans_basekv.^2)
    branch_transline[:,B]=branch_transline[:,B].*baseMVA.*inv.(trans_basekv.^2)

    #For transformer
    transformer_basekv=transformer[:,Trans_Pvoltage]
    n=findall(transformer_basekv[:].==10.0)
    transformer_basekv[n]=transformer_basekv[n].*1.05
    branch_transformer[:,R]=branch_transformer[:,R].*baseMVA.*inv.(transformer_basekv.^2)
    branch_transformer[:,X]=branch_transformer[:,X].*baseMVA.*inv.(transformer_basekv.^2)
    branch_transformer[:,B]=branch_transformer[:,B].*baseMVA.*inv.(transformer_basekv.^2)

    branch=[branch_cable;branch_transline;branch_transformer]

    HVCB_data=HVCB_data[findall(HVCB_data[:,HVCB_INSERVICE].=="Yes"),:]
    HVCB_data = HVCB_data[.!ismissing.(HVCB_data[:, HVCB_FROM_ELEMENT]), :]
    HVCB_data = HVCB_data[.!ismissing.(HVCB_data[:, HVCB_TO_ELEMENT]), :]
    HVCB_data[:,HVCB_FROM_ELEMENT]=map(k -> dict_bus[k],HVCB_data[:,HVCB_FROM_ELEMENT])
    HVCB_data[:,HVCB_TO_ELEMENT]=map(k -> dict_bus[k],HVCB_data[:,HVCB_TO_ELEMENT])

    #构造图中的点
    g = SimpleGraph(size(bus,1))  # 创建一个包含10个节点的无向图
    # 添加边（表示输电线路）
    nb=size(branch,1)
    for i in 1:nb
        add_edge!(g, Int(branch[i,FBUS]), Int(branch[i,TBUS]))
    end

    nc=size(HVCB_data,1)
    for i in 1: nc
        if(HVCB_data[i,HVCB_STATUS]=="Closed")
            add_edge!(g, Int(HVCB_data[i,HVCB_FROM_ELEMENT]), Int(HVCB_data[i,HVCB_TO_ELEMENT]))
        end
    end

    # 使用 connected_components 函数查找连通分量
components = connected_components(g)
println("连通分量: ", components)

# 识别孤岛（假设孤岛网络的节点数小于一定阈值）
threshold = 10  # 设置孤岛网络的节点数上限
isolated_subgraphs = [c for c in components if length(c) < threshold]
println("孤岛网络: ", isolated_subgraphs)

#删除孤岛节点
ng=length(isolated_subgraphs)
delete_bus_isolated=Int[]
for i in 1:ng
    delete_bus_isolated=vcat(delete_bus_isolated,isolated_subgraphs[i])
end
if !isempty(delete_bus_isolated)
    delete_bus_isolated_set = unique(delete_bus_isolated)  # 去重，避免重复删除
    bus = bus[findall(x -> !(x in delete_bus_isolated_set), bus[:, BUS_I]), :]
end
#删除孤岛支路
delete_branch_isolated=Int[]
nb=size(branch,1)
for i in 1:ng
    if(length(isolated_subgraphs[i])>=2)
        for j in 1:length(isolated_subgraphs[i])
            for k in 1:length(isolated_subgraphs[i])
                for m in 1:nb
                    if((branch[m,FBUS]==isolated_subgraphs[i][j]&&branch[m,TBUS]==isolated_subgraphs[i][k])||(branch[m,TBUS]==isolated_subgraphs[i][j]&&branch[m,FBUS]==isolated_subgraphs[i][k]))
                        push!(delete_branch_isolated,m)
                    end
                end
            end
        end
    end
end
if !isempty(delete_branch_isolated)
    delete_branch_isolated_set=unique(delete_branch_isolated)
    remaining_rows = setdiff(1:size(branch, 1), delete_branch_isolated_set)
    branch = branch[remaining_rows, :]
end

#删除HVCB
# delete_HVCB_isolated=Int[]
# nc=size(HVCB_data,1)
# for i in 1:ng
#     if(length(isolated_subgraphs[i])>=2)
#         for j in 1:length(isolated_subgraphs[i])
#             for k in 1:length(isolated_subgraphs[i])
#                 for m in 1:nc
#                     if((HVCB_data[m,HVCB_FROM_ELEMENT]==isolated_subgraphs[i][j]&&HVCB_data[m,HVCB_TO_ELEMENT]==isolated_subgraphs[i][k])||(HVCB_data[m,HVCB_TO_ELEMENT]==isolated_subgraphs[i][j]&&HVCB_data[m,HVCB_FROM_ELEMENT]==isolated_subgraphs[i][k]))
#                         push!(delete_HVCB_isolated,m)
#                     end
#                 end
#             end
#         end
#     end
# end
# if !isempty(delete_HVCB_isolated)
#     delete_HVCB_isolated_set=unique(delete_HVCB_isolated)
#     remaining_rows = setdiff(1:size(HVCB_data, 1), delete_HVCB_isolated_set)
#     HVCB_data = HVCB_data[remaining_rows, :]
# end