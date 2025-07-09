"""
    load_julia_power_data(file_path::String)

从Excel文件加载电力系统数据，并将其转换为JuliaPowerCase结构体。

参数:
- `file_path::String`: Excel文件的路径

返回:
- `JuliaPowerCase`: 包含所有电力系统组件的结构体
"""
function load_julia_power_data(file_path::String)
    # 验证文件是否存在
    if !isfile(file_path)
        error("文件不存在: $file_path")
    end
    
    # 创建空的JuliaPowerCase结构体
    case = JuliaPowerCase()
    
    try
        # 读取Excel文件
        @info "正在读取Excel文件: $file_path"
        xf = XLSX.readxlsx(file_path)
        
        # 获取所有工作表名称
        sheet_names = XLSX.sheetnames(xf)
        @info "发现工作表: $(join(sheet_names, ", "))"
        
        # 定义组件加载函数和对应的工作表名
        component_loaders = [
            ("bus", load_buses!),
            ("dcbus", load_buses!),
            ("cable", load_lines!),
            ("xline", load_lines!),
            ("dcimpedance", load_dclines!),
            ("sgen", load_static_generators!),
            ("lumpedload", load_loads!),
            ("dclumpload", load_dcloads!),
            ("xform2w", load_trafo!),
            ("xform3w", load_trafo3ws!),
            ("gen", load_generators!),
            ("battery", load_storages!),
            ("inverter", load_converters!),
            ("vpp", load_virtual_power_plants!),
            ("util", load_ext_grids!),
            ("hvcb", load_switches! ),
            ("pvarray", load_pv_arrays!),
            ("pvarray", load_ac_pv_system!),
            # ("switch", load_switches!), # To do
            # ("equipment_carbon", load_carbons!), # To do 
            # ("carbon_time_series", load_carbon_time_series!), # To do
            # ("carbon_scenario", load_carbon_scenarios!), # To do 
            # ("charging_station", load_charging_stations!), # To do
            # ("charger", load_chargers!), # To do
            # ("v2g_service", load_v2g_services!), # To do
            # ("mobile_storage", load_mess!), # To do 
            # ("load_time_series", load_time_series!), # To do
            # ("res_time_series", res_time_series!), # To do
        ]
        
        # 加载各组件数据
        for (sheet_name, loader_func) in component_loaders
            if sheet_name in sheet_names
                try
                    @info "正在加载 $sheet_name 数据..."
                    loader_func(case, file_path, sheet_name)
                catch e
                    @error "加载 $sheet_name 数据时出错" exception=(e, catch_backtrace())
                end
            else
                @debug "$sheet_name 工作表不存在，跳过"
            end
        end
        
        # 验证加载的数据
        validate_case(case)
        
    catch e
        @error "加载电力系统数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
    
    # 返回填充好的JuliaPowerCase结构体
    return case
end

"""
    validate_case(case::JuliaPowerCase)

验证加载的电力系统案例数据的完整性和一致性。
"""
function validate_case(case::JuliaPowerCase)
    # 检查是否有母线数据
    if isempty(case.busesAC)
        @warn "警告: 未加载任何母线数据"
    end
    
    # 检查线路和母线的关联
    if !isempty(case.branchesAC)
        bus_indices = Set(bus.index for bus in case.busesAC)
        for line in case.branchesAC
            if !(line.from_bus in bus_indices) || !(line.to_bus in bus_indices)
                @warn "警告: 线路 $(line.name) (ID: $(line.index)) 连接到不存在的母线"
            end
        end
    end
    
    # 可以添加更多验证逻辑...
    
    @info "电力系统案例数据验证完成"
end

"""
    safe_get_value(cell, default_value, type_converter=identity)

安全地从单元格获取值，提供类型转换和默认值。

参数:
- `cell`: Excel单元格值
- `default_value`: 如果单元格为空或转换失败时的默认值
- `type_converter`: 类型转换函数或目标类型

返回:
- 转换后的值或默认值
"""
function safe_get_value(cell, default_value, type_converter=identity)
    if ismissing(cell) || cell === nothing || (typeof(cell) <: AbstractString && isempty(strip(string(cell))))
        return default_value
    else
        try
            # 检查type_converter是类型还是函数
            if type_converter isa DataType
                # 如果是类型，创建一个适当的转换函数
                if type_converter <: Number && typeof(cell) <: AbstractString
                    return parse(type_converter, cell)
                else
                    return convert(type_converter, cell)
                end
            else
                # 如果是函数，直接使用
                return type_converter(cell)
            end
        catch e
            @debug "值转换失败: $cell 转换为 $(typeof(default_value)) 类型" exception=e
            return default_value
        end
    end
end



"""
    parse_bool(value)

从各种类型解析布尔值。

参数:
- `value`: 要解析的值

返回:
- 解析后的布尔值
"""
function parse_bool(value)
    if typeof(value) <: Bool
        return value
    elseif typeof(value) <: AbstractString
        lowercase_value = lowercase(strip(value))
        if lowercase_value in ["true", "yes", "1", "t", "y"]
            return true
        elseif lowercase_value in ["false", "no", "0", "f", "n"]
            return false
        else
            @debug "无法解析布尔值: $value 默认为false"
            return false
        end
    elseif typeof(value) <: Number
        return value != 0
    else
        @debug "无法解析布尔值类型: $(typeof(value)) 默认为false"
        return false
    end
end

"""
    load_buses!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载母线数据并添加到电力系统案例中。
同时创建母线名称到整数ID的映射以及区域名称到整数ID的映射。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含母线数据的工作表名称
"""
function load_buses!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        # 确保数据不为空
        if isempty(df)
            @info "母线表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        if sheet_name == "bus"
            required_columns = [:index, :id, :nominalkv]
        elseif sheet_name == "dcbus"
            required_columns = [:index, :id, :nominalv]
        end
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "母线表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # (1) 将df[:name]映射为从1开始的连续整数
        name_to_id = Dict{String, Int}()
        for (i, name) in enumerate(df[:, :id])
            if !ismissing(name) && !isempty(strip(string(name)))
                name_str = string(name)
                if !haskey(name_to_id, name_str)
                    name_to_id[name_str] = i
                end
            end
        end
        
        
        # (2) 将df[:zone]映射为从1开始的连续整数向量
        zone_to_id = Dict{String, Int}()
        
        # 检查是否存在zone列
        has_zone_column = any(col -> lowercase(string(col)) == "zone", names(df))
        
        if has_zone_column
            # 首先为每个唯一的区域分配连续ID
            unique_zones = []
            for zone in df[:, :zone]
                zone_str = safe_get_value(zone, "", String)
                if !isempty(zone_str) && !(zone_str in unique_zones)
                    push!(unique_zones, zone_str)
                end
            end
            
            # 为唯一的区域分配从1开始的连续ID
            for (i, zone) in enumerate(unique_zones)
                zone_to_id[zone] = i
            end
        else
            @info "母线表格中没有zone列，跳过区域映射创建"
        end
        
        #(3) 将df[:area]映射为整数ID
        area_to_id = Dict{String, Int}()
        has_area_column = any(col -> lowercase(string(col)) == "area", names(df))
        if has_area_column
            # 首先为每个唯一的区域分配连续ID
            unique_areas = []
            for area in df[:, :area]
                area_str = safe_get_value(area, "", String)
                if !isempty(area_str) && !(area_str in unique_areas)
                    push!(unique_areas, area_str)
                end
            end
            
            # 为唯一的区域分配从1开始的连续ID
            for (i, area) in enumerate(unique_areas)
                area_to_id[area] = i
            end
        else
            @info "母线表格中没有area列，跳过区域映射创建"
        end

        # 将映射保存到case中
        if sheet_name == "bus"
            case.bus_name_to_id = name_to_id
        elseif sheet_name == "dcbus"
            # 对于直流母线，使用不同的映射
            case.busdc_name_to_id = name_to_id
        end
        case.zone_to_id = zone_to_id
        case.area_to_id = area_to_id
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # Based on the index information to check whether it exists or not
                if index <= 0
                    @warn "行 $i: 无效的母线索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                bus_id = name_to_id[name]
                zone = haskey(row, :zone) ? safe_get_value(row[:zone], "", String) : ""
                area = haskey(row, :area) ? safe_get_value(row[:area], "", String) : ""
                if sheet_name == "bus"
                    vn_kv = safe_get_value(row[:nominalkv], 0.0, Float64)
                elseif sheet_name == "dcbus"
                    vn_kv = safe_get_value(row[:nominalv], 0.0, Float64)/1000.0
                end

                if has_zone_column
                    zone_id = zone_to_id[zone]
                else
                    zone_id = 1
                end
                if has_area_column
                    area_id = area_to_id[area]
                else
                    area_id = 1
                end
                # 验证电压等级是否合理
                if vn_kv <= 0.0
                    @warn "行 $i: 母线 $name (ID: $index) 的电压等级无效 ($vn_kv kV)，使用默认值 0.4 kV"
                    vn_kv = 0.4
                end
                
                max_vm_pu = haskey(row, :max_vm_pu) ? safe_get_value(row[:max_vm_pu], 1.05, Float64) : 1.05
                min_vm_pu = haskey(row, :min_vm_pu) ? safe_get_value(row[:min_vm_pu], 0.95, Float64) : 0.95
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 验证电压限制是否合理
                if min_vm_pu >= max_vm_pu
                    @warn "行 $i: 母线 $name (ID: $index) 的电压限制无效 (min: $min_vm_pu, max: $max_vm_pu)，使用默认值"
                    min_vm_pu = 0.95
                    max_vm_pu = 1.05
                end
                
                # 创建Bus对象并添加到case中
                if sheet_name == "bus"
                    push!(case.busesAC, Bus(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id))
                elseif sheet_name == "dcbus"
                    push!(case.busesDC, BusDC(index, name, vn_kv, max_vm_pu, min_vm_pu, in_service, bus_id, area_id, zone_id))
                    
                end
                processed_rows += 1

            catch e
                @error "处理母线数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "母线数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        if sheet_name == "bus"
            @info "创建了 $(length(case.bus_name_to_id)) 个母线名称到ID的映射"
        else
            @info "创建了 $(length(case.busdc_name_to_id)) 个母线名称到ID的映射"
        end
        if has_zone_column
            @info "创建了 $(length(case.zone_to_id)) 个区域名称到ID的映射"
        end
        
    catch e
        @error "加载母线数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end



"""
    load_lines!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载线路数据并添加到电力系统案例中。
使用case.bus_name_to_id将母线名称映射为整数ID。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含线路数据的工作表名称
"""
function load_lines!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取线路数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "线路表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        if sheet_name == "cable"
            required_columns = [:index, :frombus, :tobus, :lengthvalue, :cablelengthunit, :ohmsperlengthunit, :ohmsperlengthvalue,:rposvalue, :xposvalue,:rzerovalue, :xzerovalue]
        elseif sheet_name == "xline"
            required_columns = [:index, :frombus, :tobus, :length, :lengthunit, :perlength, :perlengthunit, :rpos, :xpos,:rzero, :xzero]
        end
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "线路表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 检查bus_name_to_id映射是否存在
        if !isdefined(case, :bus_name_to_id) || isempty(case.bus_name_to_id)
            @warn "case.bus_name_to_id映射不存在或为空，无法将母线名称映射为ID"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的线路索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从数据中获取母线名称
                from_bus_name = safe_get_value(row[:frombus], "", String)
                to_bus_name = safe_get_value(row[:tobus], "", String)
                
                # 使用bus_name_to_id映射将母线名称转换为整数ID
                from_bus = 0
                to_bus = 0
                
                if haskey(case.bus_name_to_id, from_bus_name)
                    from_bus = case.bus_name_to_id[from_bus_name]
                else
                    @warn "行 $i: 线路 $name (ID: $index) 的起始母线名称 '$from_bus_name' 在映射中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.bus_name_to_id, to_bus_name)
                    to_bus = case.bus_name_to_id[to_bus_name]
                else
                    @warn "行 $i: 线路 $name (ID: $index) 的终止母线名称 '$to_bus_name' 在映射中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if from_bus <= 0 || to_bus <= 0
                    @warn "行 $i: 线路 $name (ID: $index) 连接到无效的母线ID (from: $from_bus, to: $to_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(bus.bus_id for bus in case.busesAC)
                if !(from_bus in bus_indices) || !(to_bus in bus_indices)
                    @warn "行 $i: 线路 $name (ID: $index) 连接到不存在的母线ID (from: $from_bus, to: $to_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证线路不是自环
                if from_bus == to_bus
                    @warn "行 $i: 线路 $name (ID: $index) 连接到相同的母线 ($from_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                
                if sheet_name == "cable"
                    length_value = haskey(row, :lengthvalue) ? safe_get_value(row[:lengthvalue], 0.0, Float64) : 0.0
                    length_unit = safe_get_value(row[:cablelengthunit], "", String)
                else sheet_name == "xline"
                    length_value = haskey(row, :length) ? safe_get_value(row[:length], 0.0, Float64) : 0.0
                    length_unit = safe_get_value(row[:lengthunit], "", String)
                end

                if length_unit =="0"
                    length_km = length_value * 0.0003048 # 英尺转换为公里
                elseif length_unit == "1"
                    length_km = length_value * 1.60934 # 英里转换为公里
                elseif length_unit == "2"
                    length_km = length_value * 0.001 # 米转换为公里
                elseif length_unit == "3"
                    length_km = length_value * 1.0 # 公里转换为公里
                end
                
                # 验证长度是否合理
                if length_km < 0.0
                    @warn "行 $i: 线路 $name (ID: $index) 的长度无效 ($length_km km)，设置为0"
                    length_km = 0.0
                end
                if sheet_name == "cable"
                    OhmsPerLengthUnit = safe_get_value(row[:ohmsperlengthunit], "", String)
                    OhmsPerLengthValue = safe_get_value(row[:ohmsperlengthvalue], 0.0, Float64)
                    RPosValue = safe_get_value(row[:rposvalue], 0.0, Float64)
                    XPosValue = safe_get_value(row[:xposvalue], 0.0, Float64)
                    CPosValue = safe_get_value(row[:yposvalue], 0.0, Float64)
                    RZeroValue = safe_get_value(row[:rzerovalue], 0.0, Float64)
                    XZeroValue = safe_get_value(row[:xzerovalue], 0.0, Float64)
                    CZeroValue = safe_get_value(row[:yzerovalue], 0.0, Float64)
                elseif sheet_name == "xline"
                    OhmsPerLengthUnit = safe_get_value(row[:perlengthunit], "", String)
                    OhmsPerLengthValue = safe_get_value(row[:perlength], 0.0, Float64)
                    RPosValue = safe_get_value(row[:rpos], 0.0, Float64)
                    XPosValue = safe_get_value(row[:xpos], 0.0, Float64)
                    CPosValue = safe_get_value(row[:ypos], 0.0, Float64)
                    RZeroValue = safe_get_value(row[:rzero], 0.0, Float64)
                    XZeroValue = safe_get_value(row[:xzero], 0.0, Float64)
                    CZeroValue = safe_get_value(row[:yzero], 0.0, Float64)
                end

                if OhmsPerLengthUnit == "0"
                    ohmsperkm = OhmsPerLengthValue * 0.0003048 # 英尺转换为公里
                elseif OhmsPerLengthUnit == "1"
                    ohmsperkm = OhmsPerLengthValue * 1.60934 # 英里转换为公里
                elseif OhmsPerLengthUnit == "2"
                    ohmsperkm = OhmsPerLengthValue * 0.001 # 米转换为公里
                elseif OhmsPerLengthUnit == "3"
                    ohmsperkm = OhmsPerLengthValue * 1.0 # 公里转换为公里
                end

                r_ohm_per_km = RPosValue/ohmsperkm
                x_ohm_per_km = XPosValue/ohmsperkm
                c_nf_per_km = CPosValue/ohmsperkm
                r0_ohm_per_km = RZeroValue/ohmsperkm
                x0_ohm_per_km = XZeroValue/ohmsperkm
                c0_nf_per_km = CZeroValue/ohmsperkm
                g_us_per_km = CZeroValue/ohmsperkm
                
                
                # 验证阻抗是否合理
                if r_ohm_per_km < 0.0
                    @warn "行 $i: 线路 $name (ID: $index) 的电阻值无效 ($r_ohm_per_km Ω/km)，使用默认值 0.1"
                    r_ohm_per_km = 0.1
                end
                
                if x_ohm_per_km < 0.0
                    @warn "行 $i: 线路 $name (ID: $index) 的电抗值无效 ($x_ohm_per_km Ω/km)，使用默认值 0.1"
                    x_ohm_per_km = 0.1
                end
                
                max_i_ka = haskey(row, :max_i_ka) ? safe_get_value(row[:max_i_ka], 0.0, Float64) : 0.0
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                max_loading_percent = haskey(row, :max_loading_percent) ? safe_get_value(row[:max_loading_percent], 100.0, Float64) : 100.0
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # 验证并行数量是否合理
                if parallel <= 0
                    @warn "行 $i: 线路 $name (ID: $index) 的并行数量无效 ($parallel)，设置为1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 收集可靠性参数作为命名参数
                reliability_params = Dict{Symbol, Any}()
                
                # 定义可靠性参数列表及其默认值
                reliability_fields = [
                    (:mtbf_hours, 0.0),
                    (:mttr_hours, 0.0),
                    (:failure_rate_per_year, 0.0),
                    (:planned_outage_hours_per_year, 0.0),
                    (:forced_outage_rate, 0.0),
                    (:permanent_fault_rate_per_km_year, 0.0),
                    (:temporary_fault_rate_per_km_year, 0.0),
                    (:repair_time_permanent_hours, 0.0),
                    (:auto_reclosing_success_rate, 0.0)
                ]
                
                # 安全地提取可靠性参数
                for (field, default_value) in reliability_fields
                    field_str = String(field)
                    if haskey(row, Symbol(field_str))
                        reliability_params[field] = safe_get_value(row[Symbol(field_str)], default_value, Float64)
                    else
                        reliability_params[field] = default_value
                    end
                end
                
                # 创建Line对象并添加到case中，使用命名参数传递可靠性参数
                push!(case.branchesAC, Line(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                                      c_nf_per_km, r0_ohm_per_km, x0_ohm_per_km,
                                      c0_nf_per_km, g_us_per_km, max_i_ka, type, max_loading_percent,
                                      parallel, df, in_service; reliability_params...))
                
                processed_rows += 1
                
            catch e
                @error "处理线路数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "线路数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载线路数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_dclines!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载直流线路数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含直流线路数据的工作表名称
"""
function load_dclines!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取直流线路数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "直流线路表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :id, :frombus, :tobus, :rvalue, :lvalue]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "直流线路表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的直流线路索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                from_bus_name = safe_get_value(row[:frombus], "", String)
                to_bus_name = safe_get_value(row[:tobus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                from_bus = 0
                to_bus = 0
                
                if haskey(case.busdc_name_to_id, from_bus_name)
                    from_bus = case.busdc_name_to_id[from_bus_name]
                else
                    @warn "行 $i: 直流线路 $name (ID: $index) 的起始母线名称 '$from_bus_name' 在busdc_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.busdc_name_to_id, to_bus_name)
                    to_bus = case.busdc_name_to_id[to_bus_name]
                else
                    @warn "行 $i: 直流线路 $name (ID: $index) 的终止母线名称 '$to_bus_name' 在busdc_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if from_bus <= 0 || to_bus <= 0
                    @warn "行 $i: 直流线路 $name (ID: $index) 连接到无效的母线 (from: $from_bus, to: $to_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(bus.index for bus in case.busesDC)
                if !(from_bus in bus_indices) || !(to_bus in bus_indices)
                    @warn "行 $i: 直流线路 $name (ID: $index) 连接到不存在的母线 (from: $from_bus, to: $to_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证线路不是自环
                if from_bus == to_bus
                    @warn "行 $i: 直流线路 $name (ID: $index) 连接到相同的母线 ($from_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                if sheet_name == "dcimpedance"
                    r = safe_get_value(row[:rvalue], 0.0, Float64)
                    x = safe_get_value(row[:lvalue], 0.0, Float64)
                end
                
                # if x == 0.0
                #     @warn "行 $i: 直流线路 $name (ID: $index) 的电抗值无效 ($x)，使用默认值 1e-6"
                #     x = 1e-6
                # end
                
                # 从行数据中提取其他字段值
                if sheet_name == "dcimpedance"
                    length_km = 1.0
                else
                    length_km = haskey(row, :length_km) ? safe_get_value(row[:length_km], 0.0, Float64) : 0.0
                end
                
                # 验证长度是否合理
                if length_km < 0.0
                    @warn "行 $i: 直流线路 $name (ID: $index) 的长度无效 ($length_km km)，设置为0"
                    length_km = 0.0
                end
                if sheet_name == "dcimpedance"
                    r_ohm_per_km = r
                    x_ohm_per_km = x
                else
                    r_ohm_per_km = safe_get_value(row[:rvalue], 0.0, Float64)
                end
                
                # 验证阻抗是否合理
                if r_ohm_per_km < 0.0
                    @warn "行 $i: 直流线路 $name (ID: $index) 的电阻值无效 ($r_ohm_per_km Ω/km)，使用默认值 0.1"
                    r_ohm_per_km = 0.1
                end
                
                g_us_per_km = haskey(row, :g_us_per_km) ? safe_get_value(row[:g_us_per_km], 0.0, Float64) : 0.0
                max_i_ka = haskey(row, :max_i_ka) ? safe_get_value(row[:max_i_ka], 0.0, Float64) : 0.0
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                max_loading_percent = haskey(row, :max_loading_percent) ? safe_get_value(row[:max_loading_percent], 100.0, Float64) : 100.0
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # 验证并行数量是否合理
                if parallel <= 0
                    @warn "行 $i: 直流线路 $name (ID: $index) 的并行数量无效 ($parallel)，设置为1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : false
                
                # 创建LineDC对象并添加到case中
                push!(case.branchesDC, LineDC(index, name, from_bus, to_bus, length_km, r_ohm_per_km, x_ohm_per_km,
                                          g_us_per_km, max_i_ka, type, max_loading_percent,
                                          parallel, df, in_service))
                
                processed_rows += 1
                
            catch e
                @error "处理直流线路数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "直流线路数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载直流线路数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_static_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载静态发电机数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含静态发电机数据的工作表名称
"""
function load_static_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取静态发电机数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "静态发电机表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :bus, :p_mw, :q_mvar]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "静态发电机表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的静态发电机索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "行 $i: 静态发电机 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 静态发电机 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(b.index for b in case.buses)
                if !(bus in bus_indices)
                    @warn "行 $i: 静态发电机 $name (ID: $index) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                p_mw = safe_get_value(row[:p_mw], 0.0, Float64)
                q_mvar = safe_get_value(row[:q_mvar], 0.0, Float64)
                scaling = haskey(row, :scaling) ? safe_get_value(row[:scaling], 1.0, Float64) : 1.0
                
                # 验证缩放因子是否合理
                if scaling < 0.0
                    @warn "行 $i: 静态发电机 $name (ID: $index) 的缩放因子无效 ($scaling)，设置为1.0"
                    scaling = 1.0
                end
                
                max_p_mw = haskey(row, :max_p_mw) ? safe_get_value(row[:max_p_mw], 0.0, Float64) : 0.0
                min_p_mw = haskey(row, :min_p_mw) ? safe_get_value(row[:min_p_mw], 0.0, Float64) : 0.0
                max_q_mvar = haskey(row, :max_q_mvar) ? safe_get_value(row[:max_q_mvar], 0.0, Float64) : 0.0
                min_q_mvar = haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0
                
                # 验证功率限制是否合理
                if max_p_mw < min_p_mw
                    @warn "行 $i: 静态发电机 $name (ID: $index) 的有功功率限制无效 (min: $min_p_mw, max: $max_p_mw)，交换值"
                    max_p_mw, min_p_mw = min_p_mw, max_p_mw
                end
                
                if max_q_mvar < min_q_mvar
                    @warn "行 $i: 静态发电机 $name (ID: $index) 的无功功率限制无效 (min: $min_q_mvar, max: $max_q_mvar)，交换值"
                    max_q_mvar, min_q_mvar = min_q_mvar, max_q_mvar
                end
                
                k = haskey(row, :k) ? safe_get_value(row[:k], 0.0, Float64) : 0.0
                rx = haskey(row, :rx) ? safe_get_value(row[:rx], 0.0, Float64) : 0.0
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], false)) : false
                
                # 创建StaticGenerator对象并添加到case中
                push!(case.static_generators, StaticGenerator(index, name, bus, p_mw, q_mvar, scaling, max_p_mw, min_p_mw,
                                                           max_q_mvar, min_q_mvar, k, rx, in_service, type, controllable))
                
                processed_rows += 1
                
            catch e
                @error "处理静态发电机数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "静态发电机数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载静态发电机数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_loads!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载负荷数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含负荷数据的工作表名称
"""
function load_loads!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取负荷数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "负荷表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :bus, :mva, :pf]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "负荷表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的负荷索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "行 $i: 负荷 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 负荷 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(b.index for b in case.busesAC)
                if !(bus in bus_indices)
                    @warn "行 $i: 负荷 $name (ID: $index) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                s_mva = safe_get_value(row[:mva], 0.0, Float64)/1000
                p_mw = s_mva * safe_get_value(row[:pf], 0.0, Float64)/100
                q_mvar = s_mva * sqrt(1 - (safe_get_value(row[:pf], 0.0, Float64)/100)^2)
                
                # 验证负荷模型参数
                const_p_percent = haskey(row, :mtloadpercent) ? safe_get_value(row[:mtloadpercent], 100.0, Float64) : 100.0
                const_z_percent = 100.0 - const_p_percent
                const_i_percent = haskey(row, :const_i_percent) ? safe_get_value(row[:const_i_percent], 0.0, Float64) : 0.0
                
                
                # 验证百分比总和是否为100%
                total_percent = const_z_percent + const_i_percent + const_p_percent
                if abs(total_percent - 100.0) > 1e-6
                    @warn "行 $i: 负荷 $name (ID: $index) 的负荷模型百分比总和不为100% ($total_percent%)，进行归一化"
                    if total_percent > 0
                        const_z_percent = const_z_percent * 100.0 / total_percent
                        const_i_percent = const_i_percent * 100.0 / total_percent
                        const_p_percent = const_p_percent * 100.0 / total_percent
                    else
                        const_z_percent = 0.0
                        const_i_percent = 0.0
                        const_p_percent = 100.0
                    end
                end
                
                scaling = haskey(row, :uniformscalepq) ? safe_get_value(row[:uniformscalepq], 1.0, Float64) : 1.0
                
                # 验证缩放因子是否合理
                if scaling < 0.0
                    @warn "行 $i: 负荷 $name (ID: $index) 的缩放因子无效 ($scaling)，设置为1.0"
                    scaling = 1.0
                end
                
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                type = haskey(row, :LoadType) ? safe_get_value(row[:type], "", String) : "wye"
                
                if type == "0"
                    type = "wye"
                elseif type == "1"
                    type = "delta"
                end
                
                # 创建Load对象并添加到case中
                push!(case.loadsAC, Load(index, name, bus, p_mw, q_mvar, const_z_percent, const_i_percent,
                                     const_p_percent, scaling, in_service, type))
                
                processed_rows += 1
                
            catch e
                @error "处理负荷数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "负荷数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载负荷数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_dcloads!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载负荷数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含负荷数据的工作表名称
"""
function load_dcloads!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取负荷数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "负荷表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :bus, :kw]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "负荷表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的负荷索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus = case.busdc_name_to_id[bus_name]
                else
                    @warn "行 $i: 负荷 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 负荷 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(b.index for b in case.busesDC)
                if !(bus in bus_indices)
                    @warn "行 $i: 负荷 $name (ID: $index) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                p_mw = safe_get_value(row[:kw], 0.0, Float64)/1000
                
                # 验证负荷模型参数
                const_p_percent = haskey(row, :mtloadpercent) ? safe_get_value(row[:mtloadpercent], 100.0, Float64) : 100.0
                const_z_percent = haskey(row, :staticloadpercent) ? safe_get_value(row[:staticloadpercent], 100.0, Float64) : 100.0
                const_i_percent = 100.0 - const_z_percent - const_p_percent
                
                
                # 验证百分比总和是否为100%
                total_percent = const_z_percent + const_i_percent + const_p_percent
                if abs(total_percent - 100.0) > 1e-6
                    @warn "行 $i: 负荷 $name (ID: $index) 的负荷模型百分比总和不为100% ($total_percent%)，进行归一化"
                    if total_percent > 0
                        const_z_percent = const_z_percent * 100.0 / total_percent
                        const_i_percent = const_i_percent * 100.0 / total_percent
                        const_p_percent = const_p_percent * 100.0 / total_percent
                    else
                        const_z_percent = 0.0
                        const_i_percent = 0.0
                        const_p_percent = 100.0
                    end
                end
                
                # scaling = haskey(row, :uniformscalepq) ? safe_get_value(row[:uniformscalepq], 1.0, Float64) : 1.0
                scaling = 1.0
                
                # 验证缩放因子是否合理
                if scaling < 0.0
                    @warn "行 $i: 负荷 $name (ID: $index) 的缩放因子无效 ($scaling)，设置为1.0"
                    scaling = 1.0
                end
                
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 创建Load对象并添加到case中
                push!(case.loadsDC, LoadDC(index, name, bus, p_mw, const_z_percent, const_i_percent,
                                     const_p_percent, scaling, in_service))
                
                processed_rows += 1
                
            catch e
                @error "处理负荷数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "负荷数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载负荷数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_trafo!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载两卷变压器数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含两卷变压器数据的工作表名称
"""
function load_trafo!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取两卷变压器数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "两卷变压器表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :frombus, :tobus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "两卷变压器表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的变压器索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                hv_bus_name = safe_get_value(row[:frombus], "", String)
                lv_bus_name = safe_get_value(row[:tobus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                hv_bus = 0
                lv_bus = 0
                
                if haskey(case.bus_name_to_id, hv_bus_name)
                    hv_bus = case.bus_name_to_id[hv_bus_name]
                else
                    @warn "行 $i: 变压器 $name (ID: $index) 的高压侧母线名称 '$hv_bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.bus_name_to_id, lv_bus_name)
                    lv_bus = case.bus_name_to_id[lv_bus_name]
                else
                    @warn "行 $i: 变压器 $name (ID: $index) 的低压侧母线名称 '$lv_bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if hv_bus <= 0 || lv_bus <= 0
                    @warn "行 $i: 变压器 $name (ID: $index) 连接到无效的母线 (HV: $hv_bus, LV: $lv_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(bus.index for bus in case.busesAC)
                if !(hv_bus in bus_indices) || !(lv_bus in bus_indices)
                    @warn "行 $i: 变压器 $name (ID: $index) 连接到不存在的母线 (HV: $hv_bus, LV: $lv_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证高低压母线不同
                if hv_bus == lv_bus
                    @warn "行 $i: 变压器 $name (ID: $index) 的高低压母线相同 ($hv_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                sn_mva = haskey(row, :ansimva) ? safe_get_value(row[:ansimva], 0.0, Float64)/1000 : 0.0
                
                # 验证额定容量是否合理
                if sn_mva <= 0.0
                    @warn "行 $i: 变压器 $name (ID: $index) 的额定容量无效 ($sn_mva MVA)，使用默认值 10.0 MVA"
                    sn_mva = 10.0
                end
                
                vn_hv_kv = haskey(row, :primkv) ? safe_get_value(row[:primkv], 0.0, Float64) : 0.0
                vn_lv_kv = haskey(row, :seckv) ? safe_get_value(row[:seckv], 0.0, Float64) : 0.0
                
                # 验证额定电压是否合理
                if vn_hv_kv <= 0.0
                    # 尝试从母线数据获取额定电压
                    hv_bus_obj = findfirst(b -> b.index == hv_bus, case.buses)
                    if hv_bus_obj !== nothing
                        vn_hv_kv = case.buses[hv_bus_obj].vn_kv
                        @warn "行 $i: 变压器 $name (ID: $index) 的高压侧额定电压无效，使用母线电压 $vn_hv_kv kV"
                    else
                        vn_hv_kv = 110.0
                        @warn "行 $i: 变压器 $name (ID: $index) 的高压侧额定电压无效，使用默认值 $vn_hv_kv kV"
                    end
                end
                
                if vn_lv_kv <= 0.0
                    # 尝试从母线数据获取额定电压
                    lv_bus_obj = findfirst(b -> b.index == lv_bus, case.buses)
                    if lv_bus_obj !== nothing
                        vn_lv_kv = case.buses[lv_bus_obj].vn_kv
                        @warn "行 $i: 变压器 $name (ID: $index) 的低压侧额定电压无效，使用母线电压 $vn_lv_kv kV"
                    else
                        vn_lv_kv = 10.0
                        @warn "行 $i: 变压器 $name (ID: $index) 的低压侧额定电压无效，使用默认值 $vn_lv_kv kV"
                    end
                end
                
                # 验证高低压侧电压大小关系
                if vn_hv_kv <= vn_lv_kv
                    @warn "行 $i: 变压器 $name (ID: $index) 的高压侧电压 ($vn_hv_kv kV) 不大于低压侧电压 ($vn_lv_kv kV)，已标记但继续处理"
                end
                
                # 阻抗参数
                z_percent = haskey(row, :ansiposz) ? safe_get_value(row[:ansiposz], 0.0, Float64)/100 : 0.0
                x_r = haskey(row, :ansiposxr) ? safe_get_value(row[:ansiposxr], 0.0, Float64) : 0.0

                z0_percent = haskey(row, :ansizeroz) ? safe_get_value(row[:ansizeroz], 0.0, Float64)/100 : 0.0
                x0_r0 = haskey(row, :ansizeroxoverr) ? safe_get_value(row[:ansizeroxoverr], 0.0, Float64) : 0.0

                # 验证阻抗参数是否合理
                if z_percent <= 0.0
                    @warn "行 $i: 变压器 $name (ID: $index) 的阻抗占比无效 ($z_percent%)，使用默认值 10.0%"
                    z_percent = 10.0/100
                end
                
                if x_r < 0
                    @warn "行 $i: 变压器 $name (ID: $index) 的短路电阻无效 ($x_r%)，使用默认值 20.0"
                    x_r = 20.0
                end
                
                if z0_percent <= 0.0
                    @warn "行 $i: 变压器 $name (ID: $index) 的零序阻抗无效 ($z0_percent%)，使用默认值 10.0"
                    z0_percent = 10.0/100
                end

                if x0_r0 < 0.0
                    @warn "行 $i: 变压器 $name (ID: $index) 的零序电阻无效 ($x0_r0%)，使用默认值 20.0"
                    x0_r0 = 20.0
                end
                
                # 变压器接头
                tap_neutral = haskey(row, :centertap) ? safe_get_value(row[:centertap], 0.0, Float64) : 0.0
                prim_tap = haskey(row, :primpercenttap) ? safe_get_value(row[:primpercenttap], 0.0, Float64) : 0.0
                sec_tap = haskey(row, :secpercenttap) ? safe_get_value(row[:secpercenttap], 0.0, Float64) : 0.0

                prim_tap_min = haskey(row, :primminpercentfixedtap) ? safe_get_value(row[:primminpercentfixedtap], 0.0, Float64) : 0.0
                prim_tap_max = haskey(row, :primmaxpercentfixedtap) ? safe_get_value(row[:primmaxpercentfixedtap], 0.0, Float64) : 0.0

                sec_tap_min = haskey(row, :secminpercentfixedtap) ? safe_get_value(row[:secminpercentfixedtap], 0.0, Float64) : 0.0
                sec_tap_max = haskey(row, :secmaxpercentfixedtap) ? safe_get_value(row[:secmaxpercentfixedtap], 0.0, Float64) : 0.0
                
                
                # 验证变压器接头是否合理
                if prim_tap > prim_tap_max
                    @warn "行 $i: 变压器 $name (ID: $index) 的接头超过上限，无效 ($prim_tap)，设置为0"
                    prim_tap = 0.0
                end
                
                if prim_tap < prim_tap_min
                    @warn "行 $i: 变压器 $name (ID: $index) 的接头低于下限，无效 ($prim_tap)，设置为0"
                    prim_tap = 0.0
                end
                
                if sec_tap > sec_tap_max
                    @warn "行 $i: 变压器 $name (ID: $index) 的接头超过上限，无效 ($sec_tap)，设置为0"
                    sec_tap = 0.0
                end

                if sec_tap < sec_tap_min
                    @warn "行 $i: 变压器 $name (ID: $index) 的接头低于下限，无效 ($sec_tap)，设置为0"
                    sec_tap = 0.0
                end

                # 相移角度
                phaseshifthl = haskey(row, :phaseshifthl) ? safe_get_value(row[:phaseshifthl], 0.0, Float64) : 0.0
                phaseshiftps = haskey(row, :phaseshiftps) ? safe_get_value(row[:phaseshiftps], 0.0, Float64) : 0.0

                #接线方式
                vectororwinding = haskey(row, :vectororwinding) ? safe_get_value(row[:vectororwinding], "", String) : ""
                primconnectionbutton = haskey(row, :primconnectionbutton) ? safe_get_value(row[:primconnectionbutton], "", String) : ""
                secconnectionbutton = haskey(row, :secconnectionbutton) ? safe_get_value(row[:secconnectionbutton], "", String) : ""
                primneutralconn = haskey(row, :primneutralconn) ? safe_get_value(row[:primneutralconn], "", String) : ""
                secneutralconn = haskey(row, :secneutralconn) ? safe_get_value(row[:secneutralconn], "", String) : ""

                # 首先检查是否为矢量绕组模式
                if vectororwinding == "1"
                    # 初始化主边和副边的连接类型
                    prim_conn_type = ""
                    sec_conn_type = ""
                    
                    # 确定主边连接类型
                    if primconnectionbutton == "1"
                        prim_conn_type = "D"  # D型连接
                    elseif primconnectionbutton == "0"
                        prim_conn_type = "Y"  # Y型连接
                        # 检查主边中性点是否接地
                        if primneutralconn == "1"
                            prim_conn_type = "Yn"  # 中性点接地的Y型连接
                        end
                    end
                    
                    # 确定副边连接类型
                    if secconnectionbutton == "1"
                        sec_conn_type = "d"  # D型连接
                    elseif secconnectionbutton == "0"
                        sec_conn_type = "y"  # Y型连接
                        # 检查副边中性点是否接地
                        if secneutralconn == "1"
                            sec_conn_type = "yn"  # 中性点接地的Y型连接
                        end
                    end
                    
                    # 组合形成最终的vector_group
                    vector_group = prim_conn_type * sec_conn_type
                else
                    # 如果不是矢量绕组模式，可以设置为默认值或其他处理
                    vector_group = ""
                end

                
                # 其他参数
                parallel = haskey(row, :parallel) ? safe_get_value(row[:parallel], 1, Int) : 1
                
                # 验证并行数量是否合理
                if parallel <= 0
                    @warn "行 $i: 变压器 $name (ID: $index) 的并行数量无效 ($parallel)，设置为1"
                    parallel = 1
                end
                
                df = haskey(row, :df) ? safe_get_value(row[:df], 1.0, Float64) : 1.0
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                
                # 创建Transformer2Wetap对象并添加到case中
                # 在 load_trafo! 函数中修改创建 Transformer2Wetap 对象的代码
                push!(case.transformers_2w_etap, Transformer2Wetap(index, name, "", hv_bus, lv_bus, sn_mva, vn_hv_kv, vn_lv_kv,z_percent,
                        x_r, z0_percent, x0_r0, Int(round(tap_neutral)), prim_tap=prim_tap, sec_tap=sec_tap, prim_tap_min=Int(round(prim_tap_min)), prim_tap_max=Int(round(prim_tap_max)), # 转换为 Int
                        sec_tap_min=Int(round(sec_tap_min)),sec_tap_max=Int(round(sec_tap_max)), phaseshifthl=phaseshifthl, phaseshiftps=phaseshiftps, vector_group=vector_group, parallel=parallel, 
                        df=df, in_service=in_service))
                
                processed_rows += 1
                
            catch e
                @error "处理两卷变压器数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "两卷变压器数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载两卷变压器数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end



# 加载三卷变压器数据
"""
    load_trafo3ws!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载三卷变压器数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含三卷变压器数据的工作表名称
"""
function load_trafo3ws!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取三卷变压器数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "三卷变压器表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :hv_bus, :mv_bus, :lv_bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "三卷变压器表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取基本字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的三卷变压器索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                std_type = haskey(row, :std_type) ? safe_get_value(row[:std_type], "", String) : ""
                
                # 三个绕组的母线连接
                hv_bus = safe_get_value(row[:hv_bus], 0, Int)
                mv_bus = safe_get_value(row[:mv_bus], 0, Int)
                lv_bus = safe_get_value(row[:lv_bus], 0, Int)
                
                # 验证母线索引是否有效
                if hv_bus <= 0 || mv_bus <= 0 || lv_bus <= 0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 连接到无效的母线 (HV: $hv_bus, MV: $mv_bus, LV: $lv_bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线是否存在
                bus_indices = Set(bus.index for bus in case.buses)
                if !(hv_bus in bus_indices) || !(mv_bus in bus_indices) || !(lv_bus in bus_indices)
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 连接到不存在的母线，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证三个母线不同
                if hv_bus == mv_bus || hv_bus == lv_bus || mv_bus == lv_bus
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的三个绕组连接到相同的母线，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 额定容量
                sn_hv_mva = haskey(row, :sn_hv_mva) ? safe_get_value(row[:sn_hv_mva], 0.0, Float64) : 0.0
                sn_mv_mva = haskey(row, :sn_mv_mva) ? safe_get_value(row[:sn_mv_mva], 0.0, Float64) : 0.0
                sn_lv_mva = haskey(row, :sn_lv_mva) ? safe_get_value(row[:sn_lv_mva], 0.0, Float64) : 0.0
                
                # 验证额定容量是否合理
                if sn_hv_mva <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧额定容量无效 ($sn_hv_mva MVA)，使用默认值 10.0 MVA"
                    sn_hv_mva = 10.0
                end
                
                if sn_mv_mva <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧额定容量无效 ($sn_mv_mva MVA)，使用默认值 10.0 MVA"
                    sn_mv_mva = 10.0
                end
                
                if sn_lv_mva <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧额定容量无效 ($sn_lv_mva MVA)，使用默认值 10.0 MVA"
                    sn_lv_mva = 10.0
                end
                
                # 额定电压
                vn_hv_kv = haskey(row, :vn_hv_kv) ? safe_get_value(row[:vn_hv_kv], 0.0, Float64) : 0.0
                vn_mv_kv = haskey(row, :vn_mv_kv) ? safe_get_value(row[:vn_mv_kv], 0.0, Float64) : 0.0
                vn_lv_kv = haskey(row, :vn_lv_kv) ? safe_get_value(row[:vn_lv_kv], 0.0, Float64) : 0.0
                
                # 验证额定电压是否合理
                if vn_hv_kv <= 0.0
                    # 尝试从母线数据获取额定电压
                    hv_bus_obj = findfirst(b -> b.index == hv_bus, case.buses)
                    if hv_bus_obj !== nothing
                        vn_hv_kv = case.buses[hv_bus_obj].vn_kv
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧额定电压无效，使用母线电压 $vn_hv_kv kV"
                    else
                        vn_hv_kv = 110.0
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧额定电压无效，使用默认值 $vn_hv_kv kV"
                    end
                end
                
                if vn_mv_kv <= 0.0
                    # 尝试从母线数据获取额定电压
                    mv_bus_obj = findfirst(b -> b.index == mv_bus, case.buses)
                    if mv_bus_obj !== nothing
                        vn_mv_kv = case.buses[mv_bus_obj].vn_kv
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧额定电压无效，使用母线电压 $vn_mv_kv kV"
                    else
                        vn_mv_kv = 35.0
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧额定电压无效，使用默认值 $vn_mv_kv kV"
                    end
                end
                
                if vn_lv_kv <= 0.0
                    # 尝试从母线数据获取额定电压
                    lv_bus_obj = findfirst(b -> b.index == lv_bus, case.buses)
                    if lv_bus_obj !== nothing
                        vn_lv_kv = case.buses[lv_bus_obj].vn_kv
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧额定电压无效，使用母线电压 $vn_lv_kv kV"
                    else
                        vn_lv_kv = 10.0
                        @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧额定电压无效，使用默认值 $vn_lv_kv kV"
                    end
                end
                
                # 验证电压等级大小关系
                if !(vn_hv_kv > vn_mv_kv && vn_mv_kv > vn_lv_kv)
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的电压等级关系不正确 (HV: $vn_hv_kv kV, MV: $vn_mv_kv kV, LV: $vn_lv_kv kV)，已标记但继续处理"
                end
                
                # 短路电压百分比
                vk_hv_percent = haskey(row, :vk_hv_percent) ? safe_get_value(row[:vk_hv_percent], 0.0, Float64) : 0.0
                vk_mv_percent = haskey(row, :vk_mv_percent) ? safe_get_value(row[:vk_mv_percent], 0.0, Float64) : 0.0
                vk_lv_percent = haskey(row, :vk_lv_percent) ? safe_get_value(row[:vk_lv_percent], 0.0, Float64) : 0.0
                
                # 验证短路电压是否合理
                if vk_hv_percent <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧短路电压无效 ($vk_hv_percent%)，使用默认值 6.0%"
                    vk_hv_percent = 6.0
                end
                
                if vk_mv_percent <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧短路电压无效 ($vk_mv_percent%)，使用默认值 6.0%"
                    vk_mv_percent = 6.0
                end
                
                if vk_lv_percent <= 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧短路电压无效 ($vk_lv_percent%)，使用默认值 6.0%"
                    vk_lv_percent = 6.0
                end
                
                # 短路损耗百分比
                vkr_hv_percent = haskey(row, :vkr_hv_percent) ? safe_get_value(row[:vkr_hv_percent], 0.0, Float64) : 0.0
                vkr_mv_percent = haskey(row, :vkr_mv_percent) ? safe_get_value(row[:vkr_mv_percent], 0.0, Float64) : 0.0
                vkr_lv_percent = haskey(row, :vkr_lv_percent) ? safe_get_value(row[:vkr_lv_percent], 0.0, Float64) : 0.0
                
                # 验证短路损耗是否合理
                if vkr_hv_percent < 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧短路损耗无效 ($vkr_hv_percent%)，设置为0"
                    vkr_hv_percent = 0.0
                end
                
                if vkr_mv_percent < 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧短路损耗无效 ($vkr_mv_percent%)，设置为0"
                    vkr_mv_percent = 0.0
                end
                
                if vkr_lv_percent < 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧短路损耗无效 ($vkr_lv_percent%)，设置为0"
                    vkr_lv_percent = 0.0
                end
                
                # 验证短路阻抗关系
                if vkr_hv_percent > vk_hv_percent
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧短路损耗 ($vkr_hv_percent%) 大于短路电压 ($vk_hv_percent%)，调整为 $(vk_hv_percent * 0.9)%"
                    vkr_hv_percent = vk_hv_percent * 0.9
                end
                
                if vkr_mv_percent > vk_mv_percent
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧短路损耗 ($vkr_mv_percent%) 大于短路电压 ($vk_mv_percent%)，调整为 $(vk_mv_percent * 0.9)%"
                    vkr_mv_percent = vk_mv_percent * 0.9
                end
                
                if vkr_lv_percent > vk_lv_percent
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧短路损耗 ($vkr_lv_percent%) 大于短路电压 ($vk_lv_percent%)，调整为 $(vk_lv_percent * 0.9)%"
                    vkr_lv_percent = vk_lv_percent * 0.9
                end
                
                # 铁损和空载电流
                pfe_kw = haskey(row, :pfe_kw) ? safe_get_value(row[:pfe_kw], 0.0, Float64) : 0.0
                i0_percent = haskey(row, :i0_percent) ? safe_get_value(row[:i0_percent], 0.0, Float64) : 0.0
                
                # 验证铁损和空载电流是否合理
                if pfe_kw < 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的铁损无效 ($pfe_kw kW)，设置为0"
                    pfe_kw = 0.0
                end
                
                if i0_percent < 0.0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的空载电流无效 ($i0_percent%)，设置为0"
                    i0_percent = 0.0
                end
                
                # 相移角度
                shift_mv_degree = haskey(row, :shift_mv_degree) ? safe_get_value(row[:shift_mv_degree], 0.0, Float64) : 0.0
                shift_lv_degree = haskey(row, :shift_lv_degree) ? safe_get_value(row[:shift_lv_degree], 0.0, Float64) : 0.0
                
                # 分接头参数
                tap_side = haskey(row, :tap_side) ? safe_get_value(row[:tap_side], "", String) : ""
                
                # 验证分接头侧是否有效
                if !isempty(tap_side) && !(tap_side in ["hv", "mv", "lv"])
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的分接头侧无效 ($tap_side)，设置为hv"
                    tap_side = "hv"
                end
                
                tap_neutral = haskey(row, :tap_neutral) ? safe_get_value(row[:tap_neutral], 0, Int) : 0
                tap_min = haskey(row, :tap_min) ? safe_get_value(row[:tap_min], 0, Int) : 0
                tap_max = haskey(row, :tap_max) ? safe_get_value(row[:tap_max], 0, Int) : 0
                
                # 验证分接头位置范围
                if tap_min > tap_max && tap_max != 0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的分接头位置范围无效 (min: $tap_min, max: $tap_max)，交换值"
                    tap_min, tap_max = tap_max, tap_min
                end
                
                tap_step_percent = haskey(row, :tap_step_percent) ? safe_get_value(row[:tap_step_percent], 0.0, Float64) : 0.0
                tap_step_degree = haskey(row, :tap_step_degree) ? safe_get_value(row[:tap_step_degree], 0.0, Float64) : 0.0
                tap_at_star_point = haskey(row, :tap_at_star_point) ? parse_bool(safe_get_value(row[:tap_at_star_point], false)) : false
                tap_pos = haskey(row, :tap_pos) ? safe_get_value(row[:tap_pos], 0, Int) : 0
                
                # 验证分接头位置是否在范围内
                if tap_pos < tap_min && tap_min != 0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的当前分接头位置 ($tap_pos) 小于最小值 ($tap_min)，设置为最小值"
                    tap_pos = tap_min
                elseif tap_pos > tap_max && tap_max != 0
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的当前分接头位置 ($tap_pos) 大于最大值 ($tap_max)，设置为最大值"
                    tap_pos = tap_max
                end
                
                # 运行状态
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                
                # 技术参数
                vector_group_hv_mv = haskey(row, :vector_group_hv_mv) ? safe_get_value(row[:vector_group_hv_mv], "", String) : ""
                vector_group_hv_lv = haskey(row, :vector_group_hv_lv) ? safe_get_value(row[:vector_group_hv_lv], "", String) : ""
                vector_group_mv_lv = haskey(row, :vector_group_mv_lv) ? safe_get_value(row[:vector_group_mv_lv], "", String) : ""
                hv_connection = haskey(row, :hv_connection) ? safe_get_value(row[:hv_connection], "", String) : ""
                mv_connection = haskey(row, :mv_connection) ? safe_get_value(row[:mv_connection], "", String) : ""
                lv_connection = haskey(row, :lv_connection) ? safe_get_value(row[:lv_connection], "", String) : ""
                thermal_capacity_mw = haskey(row, :thermal_capacity_mw) ? safe_get_value(row[:thermal_capacity_mw], 0.0, Float64) : 0.0
                cooling_type = haskey(row, :cooling_type) ? safe_get_value(row[:cooling_type], "", String) : ""
                oil_volume_liters = haskey(row, :oil_volume_liters) ? safe_get_value(row[:oil_volume_liters], 0.0, Float64) : 0.0
                winding_material = haskey(row, :winding_material) ? safe_get_value(row[:winding_material], "", String) : ""
                
                # 验证接线方式是否有效
                valid_connections = ["Y", "D", "Z", "y", "d", "z", ""]
                if !(hv_connection in valid_connections)
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的高压侧接线方式无效 ($hv_connection)，设置为Y"
                    hv_connection = "Y"
                end
                
                if !(mv_connection in valid_connections)
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的中压侧接线方式无效 ($mv_connection)，设置为Y"
                    mv_connection = "Y"
                end
                
                if !(lv_connection in valid_connections)
                    @warn "行 $i: 三卷变压器 $name (ID: $index) 的低压侧接线方式无效 ($lv_connection)，设置为Y"
                    lv_connection = "Y"
                end
                
                # 创建额外参数的字典
                kwargs = Dict{Symbol, Any}(
                    :tap_side => tap_side,
                    :tap_neutral => tap_neutral,
                    :tap_min => tap_min,
                    :tap_max => tap_max,
                    :tap_step_percent => tap_step_percent,
                    :tap_step_degree => tap_step_degree,
                    :tap_at_star_point => tap_at_star_point,
                    :tap_pos => tap_pos,
                    :in_service => in_service,
                    :vector_group_hv_mv => vector_group_hv_mv,
                    :vector_group_hv_lv => vector_group_hv_lv,
                    :vector_group_mv_lv => vector_group_mv_lv,
                    :hv_connection => hv_connection,
                    :mv_connection => mv_connection,
                    :lv_connection => lv_connection,
                    :thermal_capacity_mw => thermal_capacity_mw,
                    :cooling_type => cooling_type,
                    :oil_volume_liters => oil_volume_liters,
                    :winding_material => winding_material
                )
                
                # 创建Transformer3W对象并添加到case中
                try
                    trafo3w = Transformer3W(
                        index, name, std_type, hv_bus, mv_bus, lv_bus, 
                        sn_hv_mva, sn_mv_mva, sn_lv_mva, 
                        vn_hv_kv, vn_mv_kv, vn_lv_kv,
                        vk_hv_percent, vk_mv_percent, vk_lv_percent,
                        vkr_hv_percent, vkr_mv_percent, vkr_lv_percent,
                        pfe_kw, i0_percent, shift_mv_degree, shift_lv_degree; 
                        kwargs...
                    )
                    
                    push!(case.trafo3ws, trafo3w)
                    processed_rows += 1
                catch e
                    @error "创建三卷变压器对象时出错" exception=(e, catch_backtrace()) transformer_data=(index=index, name=name)
                    error_rows += 1
                end
                
            catch e
                @error "处理三卷变压器数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "三卷变压器数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载三卷变压器数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载发电机数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含发电机数据的工作表名称
"""
function load_generators!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取发电机数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "发电机表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :name, :bus, :p_mw]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "发电机表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取基本字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的发电机索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:name], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "行 $i: 发电机 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 发电机 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.index == bus, case.buses)
                    @warn "行 $i: 发电机 $name (ID: $index) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                p_mw = safe_get_value(row[:p_mw], 0.0, Float64)
                vm_pu = haskey(row, :vm_pu) ? safe_get_value(row[:vm_pu], 1.0, Float64) : 1.0
                
                # 验证电压值是否合理
                if vm_pu <= 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的电压设定值无效 ($vm_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_pu = 1.0
                end
                
                sn_mva = haskey(row, :sn_mva) ? safe_get_value(row[:sn_mva], 0.0, Float64) : 0.0
                
                # 验证额定容量是否合理
                if sn_mva < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的额定容量无效 ($sn_mva MVA)，设置为 0 MVA"
                    sn_mva = 0.0
                end
                
                min_p_mw = haskey(row, :min_p_mw) ? safe_get_value(row[:min_p_mw], 0.0, Float64) : 0.0
                max_p_mw = haskey(row, :max_p_mw) ? safe_get_value(row[:max_p_mw], 0.0, Float64) : 0.0
                min_q_mvar = haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0
                max_q_mvar = haskey(row, :max_q_mvar) ? safe_get_value(row[:max_q_mvar], 0.0, Float64) : 0.0
                
                # 验证功率限制是否合理
                if max_p_mw < min_p_mw && max_p_mw != 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的有功功率限制无效 (min: $min_p_mw MW, max: $max_p_mw MW)，交换值"
                    max_p_mw, min_p_mw = min_p_mw, max_p_mw
                end
                
                if max_q_mvar < min_q_mvar && max_q_mvar != 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的无功功率限制无效 (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar)，交换值"
                    max_q_mvar, min_q_mvar = min_q_mvar, max_q_mvar
                end
                
                scaling = haskey(row, :scaling) ? safe_get_value(row[:scaling], 1.0, Float64) : 1.0
                
                # 验证缩放因子是否合理
                if scaling < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的缩放因子无效 ($scaling)，设置为默认值 1.0"
                    scaling = 1.0
                end
                
                slack = haskey(row, :slack) ? parse_bool(safe_get_value(row[:slack], false)) : false
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # 经济参数
                cost_a = haskey(row, :cost_a) ? safe_get_value(row[:cost_a], 0.0, Float64) : 0.0
                cost_b = haskey(row, :cost_b) ? safe_get_value(row[:cost_b], 0.0, Float64) : 0.0
                cost_c = haskey(row, :cost_c) ? safe_get_value(row[:cost_c], 0.0, Float64) : 0.0
                startup_cost = haskey(row, :startup_cost) ? safe_get_value(row[:startup_cost], 0.0, Float64) : 0.0
                shutdown_cost = haskey(row, :shutdown_cost) ? safe_get_value(row[:shutdown_cost], 0.0, Float64) : 0.0
                
                # 验证成本系数是否合理
                if cost_a < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的二次成本系数无效 ($cost_a)，设置为 0"
                    cost_a = 0.0
                end
                
                if cost_b < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的一次成本系数无效 ($cost_b)，设置为 0"
                    cost_b = 0.0
                end
                
                # 技术参数
                min_up_time = haskey(row, :min_up_time) ? safe_get_value(row[:min_up_time], 0, Int) : 0
                min_down_time = haskey(row, :min_down_time) ? safe_get_value(row[:min_down_time], 0, Int) : 0
                ramp_up_mw_per_min = haskey(row, :ramp_up_mw_per_min) ? safe_get_value(row[:ramp_up_mw_per_min], 0.0, Float64) : 0.0
                ramp_down_mw_per_min = haskey(row, :ramp_down_mw_per_min) ? safe_get_value(row[:ramp_down_mw_per_min], 0.0, Float64) : 0.0
                startup_time = haskey(row, :startup_time) ? safe_get_value(row[:startup_time], 0, Int) : 0
                
                # 验证技术参数是否合理
                if min_up_time < 0
                    @warn "行 $i: 发电机 $name (ID: $index) 的最小启动时间无效 ($min_up_time)，设置为 0"
                    min_up_time = 0
                end
                
                if min_down_time < 0
                    @warn "行 $i: 发电机 $name (ID: $index) 的最小停机时间无效 ($min_down_time)，设置为 0"
                    min_down_time = 0
                end
                
                if ramp_up_mw_per_min < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的爬坡上升速率无效 ($ramp_up_mw_per_min MW/min)，设置为 0"
                    ramp_up_mw_per_min = 0.0
                end
                
                if ramp_down_mw_per_min < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的爬坡下降速率无效 ($ramp_down_mw_per_min MW/min)，设置为 0"
                    ramp_down_mw_per_min = 0.0
                end
                
                if startup_time < 0
                    @warn "行 $i: 发电机 $name (ID: $index) 的启动时间无效 ($startup_time)，设置为 0"
                    startup_time = 0
                end
                
                # 可靠性参数
                mtbf_hours = haskey(row, :mtbf_hours) ? safe_get_value(row[:mtbf_hours], 0.0, Float64) : 0.0
                mttr_hours = haskey(row, :mttr_hours) ? safe_get_value(row[:mttr_hours], 0.0, Float64) : 0.0
                failure_rate_per_year = haskey(row, :failure_rate_per_year) ? safe_get_value(row[:failure_rate_per_year], 0.0, Float64) : 0.0
                planned_outage_hours_per_year = haskey(row, :planned_outage_hours_per_year) ? safe_get_value(row[:planned_outage_hours_per_year], 0.0, Float64) : 0.0
                forced_outage_rate = haskey(row, :forced_outage_rate) ? safe_get_value(row[:forced_outage_rate], 0.0, Float64) : 0.0
                
                # 验证可靠性参数是否合理
                if mtbf_hours < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的平均故障间隔时间无效 ($mtbf_hours 小时)，设置为 0"
                    mtbf_hours = 0.0
                end
                
                if mttr_hours < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的平均修复时间无效 ($mttr_hours 小时)，设置为 0"
                    mttr_hours = 0.0
                end
                
                if failure_rate_per_year < 0.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的年故障率无效 ($failure_rate_per_year)，设置为 0"
                    failure_rate_per_year = 0.0
                end
                
                if forced_outage_rate < 0.0 || forced_outage_rate > 1.0
                    @warn "行 $i: 发电机 $name (ID: $index) 的强制停运率无效 ($forced_outage_rate)，设置为 0"
                    forced_outage_rate = 0.0
                end
                
                # 创建Generator对象并添加到case中
                generator = Dict(
                    :index => index,
                    :name => name,
                    :bus => bus,
                    :p_mw => p_mw,
                    :vm_pu => vm_pu,
                    :sn_mva => sn_mva,
                    :min_p_mw => min_p_mw,
                    :max_p_mw => max_p_mw,
                    :min_q_mvar => min_q_mvar,
                    :max_q_mvar => max_q_mvar,
                    :scaling => scaling,
                    :slack => slack,
                    :in_service => in_service,
                    :type => type,
                    :controllable => controllable,
                    :cost_a => cost_a,
                    :cost_b => cost_b,
                    :cost_c => cost_c,
                    :startup_cost => startup_cost,
                    :shutdown_cost => shutdown_cost,
                    :min_up_time => min_up_time,
                    :min_down_time => min_down_time,
                    :ramp_up_mw_per_min => ramp_up_mw_per_min,
                    :ramp_down_mw_per_min => ramp_down_mw_per_min,
                    :startup_time => startup_time,
                    :mtbf_hours => mtbf_hours,
                    :mttr_hours => mttr_hours,
                    :failure_rate_per_year => failure_rate_per_year,
                    :planned_outage_hours_per_year => planned_outage_hours_per_year,
                    :forced_outage_rate => forced_outage_rate
                )
                
                push!(case.generators, generator)
                processed_rows += 1
                
            catch e
                @error "处理发电机数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "发电机数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
        # 检查是否有平衡节点
        if any(gen -> gen[:slack], case.generators)
            @info "系统中存在 $(count(gen -> gen[:slack], case.generators)) 个平衡节点"
        else
            @warn "系统中没有平衡节点，可能导致潮流计算无法收敛"
        end
        
    catch e
        @error "加载发电机数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
从Excel文件中加载储能设备数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含储能设备数据的工作表名称
"""
function load_storages!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取储能设备数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "储能表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "储能表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的储能设备索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus = case.busdc_name_to_id[bus_name]
                else
                    @warn "行 $i: 储能设备 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 储能设备 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.index == bus, case.busesDC)
                    @warn "行 $i: 储能设备 $name (ID: $index) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # ETAP储能设备特有参数
                ra = haskey(row, :ra) ? safe_get_value(row[:ra], 0.0, Float64) : 0.0025
                cell = haskey(row, :nrofcells) ? safe_get_value(row[:nrofcells], 0.0, Float64) : 0.0
                str = haskey(row, :nrofstrings) ? safe_get_value(row[:nrofstrings], 0.0, Float64) : 0.0
                package = haskey(row, :noofpacks) ? safe_get_value(row[:noofpacks], 0.0, Float64) : 0.0
                
                # 验证参数是否合理
                if ra < 0.0
                    @warn "行 $i: 储能设备 $name (ID: $index) 的ra参数无效 ($ra)，设置为默认值 0.0"
                    ra = 0.0
                end
                
                if cell < 0.0
                    @warn "行 $i: 储能设备 $name (ID: $index) 的cell参数无效 ($cell)，设置为默认值 0.0"
                    cell = 0.0
                end
                
                if str < 0.0
                    @warn "行 $i: 储能设备 $name (ID: $index) 的string参数无效 ($str)，设置为默认值 0.0"
                    str = 0.0
                end
                
                if package < 0.0
                    @warn "行 $i: 储能设备 $name (ID: $index) 的package参数无效 ($package)，设置为默认值 0.0"
                    package = 0.0
                end

                vpc = haskey(row, :vpc) ? safe_get_value(row[:vpc], 0.0, Float64) : 2.06
                voc = vpc * cell *package/1000
                if voc < 0.0
                    @warn "行 $i: 储能设备 $name (ID: $index) 的voc参数无效 ($voc)，设置为默认值 0.0"
                    voc = 0.0
                end
                
                # 其他参数
                in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
                type = haskey(row, :type) ? safe_get_value(row[:type], "", String) : ""
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # 创建Storageetap对象并添加到case中
                push!(case.storageetap, Storageetap(
                    index, name, bus, ra, cell, str, package, voc, in_service, type, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "处理储能设备数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "储能设备数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载储能设备数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_converters!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件加载换流器数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含换流器数据的工作表名称
"""
function load_converters!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取换流器数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "换流器表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :busid, :cznetwork]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "换流器表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的换流器索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                parent_pv = haskey(row, :parentname) ? safe_get_value(row[:parentname], "", String) : ""
                if parent_pv !== ""
                    @warn "行 $i: 换流器 $name (ID: $index) 接入了交流光伏系统，于光伏系统加载时处理"
                    continue
                end

                # 从母线名称映射到母线ID
                bus_ac_name = safe_get_value(row[:busid], "", String)
                bus_dc_name = safe_get_value(row[:cznetwork], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus_ac = 0
                bus_dc = 0
                
                if haskey(case.bus_name_to_id, bus_ac_name)
                    bus_ac = case.bus_name_to_id[bus_ac_name]
                else
                    @warn "行 $i: 换流器 $name (ID: $index) 的交流侧母线名称 '$bus_ac_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if haskey(case.busdc_name_to_id, bus_dc_name)
                    bus_dc = case.busdc_name_to_id[bus_dc_name]
                else
                    @warn "行 $i: 换流器 $name (ID: $index) 的直流侧母线名称 '$bus_dc_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus_ac <= 0 || bus_dc <= 0
                    @warn "行 $i: 换流器 $name (ID: $index) 连接到无效的母线 (AC: $bus_ac, DC: $bus_dc)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从行数据中提取其他字段值
                p_mw = safe_get_value(row[:gencat0ackw], 0.0, Float64)/1000
                q_mvar = (haskey(row, :gencat0kvar) ? safe_get_value(row[:gencat0kvar], 0.0, Float64) : 0.0)/1000
                vm_ac_pu = haskey(row, :vm_ac_pu) ? safe_get_value(row[:vm_ac_pu], 1.0, Float64) : 1.0
                vm_dc_pu = haskey(row, :vm_dc_pu) ? safe_get_value(row[:vm_dc_pu], 1.0, Float64) : 1.0
                
                # 验证电压值是否合理
                if vm_ac_pu <= 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的交流侧电压无效 ($vm_ac_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_ac_pu = 1.0
                end
                
                if vm_dc_pu <= 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的直流侧电压无效 ($vm_dc_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_dc_pu = 1.0
                end
                
                loss_percent = 1-(haskey(row, :dcpercenteff) ? safe_get_value(row[:dcpercenteff], 0.0, Float64) : 0.0)/100
                loss_mw = haskey(row, :loss_mw) ? safe_get_value(row[:loss_mw], 0.0, Float64) : 0.0
                
                # 验证损耗是否合理
                if loss_percent < 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的损耗百分比无效 ($loss_percent%)，设置为 0%"
                    loss_percent = 0.0
                end
                
                if loss_mw < 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的损耗功率无效 ($loss_mw MW)，设置为 0 MW"
                    loss_mw = 0.0
                end
                
                max_p_mw = (haskey(row, :kwmax) ? safe_get_value(row[:kwmax], 0.0, Float64) : 0.0)/1000
                min_p_mw = (haskey(row, :kwmin) ? safe_get_value(row[:kwmin], 0.0, Float64) : 0.0)/1000
                max_q_mvar = (haskey(row, :kvarmax) ? safe_get_value(row[:kvarmax], 0.0, Float64) : 0.0)/1000
                min_q_mvar = (haskey(row, :min_q_mvar) ? safe_get_value(row[:min_q_mvar], 0.0, Float64) : 0.0)/1000
                
                # 验证功率限制是否合理
                if min_p_mw > max_p_mw && max_p_mw != 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的有功功率限制无效 (min: $min_p_mw MW, max: $max_p_mw MW)，交换值"
                    min_p_mw, max_p_mw = max_p_mw, min_p_mw
                end
                
                if min_q_mvar > max_q_mvar && max_q_mvar != 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的无功功率限制无效 (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar)，交换值"
                    min_q_mvar, max_q_mvar = max_q_mvar, min_q_mvar
                end
                
                control_mode = haskey(row, :control_mode) ? safe_get_value(row[:control_mode], "", String) : ""
                droop_kv = haskey(row, :droop_kv) ? safe_get_value(row[:droop_kv], 0.0, Float64) : 0.05
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # 创建Converter对象并添加到case中
                push!(case.converters, Converter(
                    index, name, bus_ac, bus_dc, p_mw, q_mvar, 
                    vm_ac_pu, vm_dc_pu, loss_percent, loss_mw,
                    max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                    control_mode, droop_kv, in_service, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "处理换流器数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "换流器数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载换流器数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_virtual_power_plants!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件中加载虚拟电厂(VPP)数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含虚拟电厂数据的工作表名称
"""
function load_virtual_power_plants!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    # 使用DataFrame处理
    df = DataFrame(XLSX.readtable(file_path, sheet_name))
    
    # 确保数据不为空
    if isempty(df)
        @info "虚拟电厂表格为空"
        return
    end
    
    # 将列名转换为小写
    rename!(df, lowercase.(names(df)))
    
    # 遍历每一行数据
    for row in eachrow(df)
        try
            # 从行数据中提取基本字段值
            index = safe_get_value(row[:index], 0, Int)
            name = safe_get_value(row[:name], "", String)
            description = haskey(row, :description) ? safe_get_value(row[:description], "", String) : ""
            control_area = haskey(row, :control_area) ? safe_get_value(row[:control_area], "", String) : ""
            
            # 容量和能量参数
            capacity_mw = haskey(row, :capacity_mw) ? safe_get_value(row[:capacity_mw], 0.0, Float64) : 0.0
            energy_mwh = haskey(row, :energy_mwh) ? safe_get_value(row[:energy_mwh], 0.0, Float64) : 0.0
            
            # 响应和爬坡参数
            response_time_s = haskey(row, :response_time_s) ? safe_get_value(row[:response_time_s], 0.0, Float64) : 0.0
            ramp_rate_mw_per_min = haskey(row, :ramp_rate_mw_per_min) ? safe_get_value(row[:ramp_rate_mw_per_min], 0.0, Float64) : 0.0
            availability_percent = haskey(row, :availability_percent) ? safe_get_value(row[:availability_percent], 100.0, Float64) : 100.0
            
            # 运营信息
            operator = haskey(row, :operator) ? safe_get_value(row[:operator], "", String) : ""
            in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
            
            # 收集额外的kwargs参数
            kwargs = Dict{Symbol, Any}()
            
            # 资源信息
            if haskey(row, :resource_type)
                kwargs[:resource_type] = safe_get_value(row[:resource_type], "", String)
            end
            
            if haskey(row, :resource_id)
                kwargs[:resource_id] = safe_get_value(row[:resource_id], 0, Int)
            end
            
            if haskey(row, :capacity_share_percent)
                kwargs[:capacity_share_percent] = safe_get_value(row[:capacity_share_percent], 0.0, Float64)
            end
            
            if haskey(row, :control_priority)
                kwargs[:control_priority] = safe_get_value(row[:control_priority], 0, Int)
            end
            
            if haskey(row, :resource_response_time_s)
                kwargs[:resource_response_time_s] = safe_get_value(row[:resource_response_time_s], 0.0, Float64)
            end
            
            if haskey(row, :max_duration_h)
                kwargs[:max_duration_h] = safe_get_value(row[:max_duration_h], 0.0, Float64)
            end
            
            # 负荷信息
            if haskey(row, :timestamp) && !ismissing(row[:timestamp])
                # 处理时间戳，根据实际格式调整
                if isa(row[:timestamp], String)
                    kwargs[:timestamp] = DateTime(row[:timestamp], dateformat"yyyy-mm-dd HH:MM:SS")
                elseif isa(row[:timestamp], DateTime)
                    kwargs[:timestamp] = row[:timestamp]
                else
                    kwargs[:timestamp] = DateTime(now())
                end
            end
            
            if haskey(row, :p_mw)
                kwargs[:p_mw] = safe_get_value(row[:p_mw], 0.0, Float64)
            end
            
            if haskey(row, :q_mvar)
                kwargs[:q_mvar] = safe_get_value(row[:q_mvar], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_up_mw)
                kwargs[:flexibility_up_mw] = safe_get_value(row[:flexibility_up_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_down_mw)
                kwargs[:flexibility_down_mw] = safe_get_value(row[:flexibility_down_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_duration_h)
                kwargs[:flexibility_duration_h] = safe_get_value(row[:flexibility_duration_h], 0.0, Float64)
            end
            
            # 创建VirtualPowerPlant对象并添加到case中
            vpp = VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                                   response_time_s, ramp_rate_mw_per_min, availability_percent,
                                   operator, in_service; kwargs...)
            
            push!(case.virtual_power_plants, vpp)
            
        catch e
            @warn "处理虚拟电厂时出错: $e"
            @warn "问题行: $(row)"
        end
    end
    
    @info "已加载 $(length(case.virtual_power_plants)) 个虚拟电厂"
end


"""
    load_ext_grids!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件中加载外部电网数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含外部电网数据的工作表名称
"""
function load_ext_grids!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取外部电网数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "外部电网表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :bus]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "外部电网表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的外部电网索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus = case.bus_name_to_id[bus_name]
                else
                    @warn "行 $i: 外部电网 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus <= 0
                    @warn "行 $i: 外部电网 $name (ID: $index) 连接到无效的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.bus_id == bus, case.busesAC)
                    @warn "行 $i: 外部电网 $name (ID: $bus_id) 连接到不存在的母线 ($bus)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 电压参数
                vm_pu = haskey(row, :vm_pu) ? safe_get_value(row[:vm_pu], 1.0, Float64) : 1.0
                va_degree = haskey(row, :va_degree) ? safe_get_value(row[:va_degree], 0.0, Float64) : 0.0
                
                # 验证电压值是否合理
                if vm_pu <= 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的电压设定值无效 ($vm_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_pu = 1.0
                end
                
                # 运行状态
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 短路容量参数
                s_sc_max_mva = haskey(row, :s_sc_max_mva) ? safe_get_value(row[:s_sc_max_mva], 1000.0, Float64) : 1000.0
                s_sc_min_mva = haskey(row, :s_sc_min_mva) ? safe_get_value(row[:s_sc_min_mva], 1000.0, Float64) : 1000.0
                
                # 验证短路容量是否合理
                if s_sc_max_mva < s_sc_min_mva && s_sc_min_mva > 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的短路容量范围无效 (min: $s_sc_min_mva MVA, max: $s_sc_max_mva MVA)，交换值"
                    s_sc_max_mva, s_sc_min_mva = s_sc_min_mva, s_sc_max_mva
                end
                
                if s_sc_max_mva <= 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的最大短路容量无效 ($s_sc_max_mva MVA)，设置为默认值 1000.0 MVA"
                    s_sc_max_mva = 1000.0
                end
                
                if s_sc_min_mva <= 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的最小短路容量无效 ($s_sc_min_mva MVA)，设置为默认值 1000.0 MVA"
                    s_sc_min_mva = 1000.0
                end
                
                # 阻抗比参数
                posr = haskey(row, :posr) ? safe_get_value(row[:posr], 0.1, Float64) : 0.1
                posx = haskey(row, :posx) ? safe_get_value(row[:posx], 1.0, Float64) : 1.0
                zeror = haskey(row, :zeror) ? safe_get_value(row[:zeror], 0.1, Float64) : 0.1
                zerox = haskey(row, :zerox) ? safe_get_value(row[:zerox], 1.0, Float64) : 1.0
                rx_max = posr / posx
                rx_min = posr / posx
                r0x0_max = zeror / zerox
                x0x_max = zerox / posx
                
                # 验证阻抗比是否合理
                if rx_max < 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的最大R/X比无效 ($rx_max)，设置为默认值 0.1"
                    rx_max = 0.1
                end
                
                if rx_min < 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的最小R/X比无效 ($rx_min)，设置为默认值 0.1"
                    rx_min = 0.1
                end
                
                if rx_max < rx_min && rx_min > 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的R/X比范围无效 (min: $rx_min, max: $rx_max)，交换值"
                    rx_max, rx_min = rx_min, rx_max
                end
                
                if r0x0_max < 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的零序R0/X0比无效 ($r0x0_max)，设置为默认值 0.1"
                    r0x0_max = 0.1
                end
                
                if x0x_max < 0.0
                    @warn "行 $i: 外部电网 $name (ID: $index) 的零序X0/X比无效 ($x0x_max)，设置为默认值 1.0"
                    x0x_max = 1.0
                end
                
                # 可控性
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], false)) : false
                
                # 创建ExternalGrid对象并添加到case中
                push!(case.ext_grids, ExternalGrid(
                    index, name, bus, vm_pu, va_degree, in_service,
                    s_sc_max_mva, s_sc_min_mva, rx_max, rx_min, r0x0_max, x0x_max, controllable
                ))
                
                processed_rows += 1
                
            catch e
                @error "处理外部电网数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "外部电网数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载外部电网数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_switches!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件中加载开关数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含开关数据的工作表名称
"""
function load_switches!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取开关数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "开关表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :fromelement, :toelement]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "开关表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的开关索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                # 从母线名称映射到母线ID - 起始母线
                bus_from_name = safe_get_value(row[:fromelement], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus_from = 0
                
                if haskey(case.bus_name_to_id, bus_from_name)
                    bus_from = case.bus_name_to_id[bus_from_name]
                else
                    @warn "行 $i: 开关 $name (ID: $index) 的起始母线名称 '$bus_from_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 从母线名称映射到母线ID - 终止母线
                bus_to_name = safe_get_value(row[:toelement], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus_to = 0
                
                if haskey(case.bus_name_to_id, bus_to_name)
                    bus_to = case.bus_name_to_id[bus_to_name]
                else
                    @warn "行 $i: 开关 $name (ID: $index) 的终止母线名称 '$bus_to_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus_from <= 0
                    @warn "行 $i: 开关 $name (ID: $index) 连接到无效的起始母线 ($bus_from)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if bus_to <= 0
                    @warn "行 $i: 开关 $name (ID: $index) 连接到无效的终止母线 ($bus_to)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.index == bus_from, case.busesAC)
                    @warn "行 $i: 开关 $name (ID: $index) 连接到不存在的起始母线 ($bus_from)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if !any(b -> b.index == bus_to, case.busesAC)
                    @warn "行 $i: 开关 $name (ID: $index) 连接到不存在的终止母线 ($bus_to)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证起始和终止母线不同
                if bus_from == bus_to
                    @warn "行 $i: 开关 $name (ID: $index) 的起始母线和终止母线相同 ($bus_from)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 元件类型和ID
                element_type = haskey(row, :element_type) ? safe_get_value(row[:element_type], "l", String) : "l"
                element_id = haskey(row, :element_id) ? safe_get_value(row[:element_id], 0, Int) : 0
                
                # 验证元件类型是否有效
                valid_element_types = ["l", "t", "t3", "b"]
                if !(lowercase(element_type) in valid_element_types)
                    @warn "行 $i: 开关 $name (ID: $index) 的元件类型无效 ($element_type)，设置为默认值 'l'"
                    element_type = "l"
                end
                
                # 开关状态
                closed = haskey(row, :pdestatus) ? parse_bool(safe_get_value(row[:pdestatus], true)=="Closed") : true
                
                # 开关类型和参数
                type = haskey(row, :type) ? safe_get_value(row[:type], "CB", String) : "CB"
                z_ohm = haskey(row, :z_ohm) ? safe_get_value(row[:z_ohm], 0.0, Float64) : 0.0
                
                # 验证阻抗值是否合理
                if z_ohm < 0.0
                    @warn "行 $i: 开关 $name (ID: $index) 的阻抗值无效 ($z_ohm Ω)，设置为默认值 0.0 Ω"
                    z_ohm = 0.0
                end
                
                # 运行状态
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 创建Switch对象并添加到case中
                push!(case.hvcbs, HighVoltageCircuitBreaker(
                    index, name, bus_from, bus_to, element_type, element_id,
                    closed, type, z_ohm, in_service
                ))
                
                # 更新索引映射
                # case.switch_indices[index] = length(case.hvcbs)
                
                processed_rows += 1
                
            catch e
                @error "处理开关数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "开关数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载开关数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_pv_arrays!(case::JuliaPowerCase, file_path::String, sheet_name::String)
从Excel文件中加载光伏板(pv array)数据并添加到电力系统案例中。
参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含虚拟电厂数据的工作表名称
"""
function load_pv_arrays!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取光伏阵列数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        
        # 确保数据不为空
        if isempty(df)
            @info "光伏阵列表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :id, :bus, :numpanelseries, :numpanelparallel, :vmpp, :impp, :voc, :isc, :pvanumcells]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "光伏阵列表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的光伏阵列索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)

                inverterinclude = haskey(row, :inverterincluded) ? safe_get_value(row[:inverterincluded], "", String) : "0"
                if  inverterinclude == "1"
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的光伏阵列用于构建交流光伏系统，跳过此行"
                    continue
                end
                
                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus_id = 0
                
                if haskey(case.busdc_name_to_id, bus_name)
                    bus_id = case.busdc_name_to_id[bus_name]
                else
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus_id <= 0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 连接到无效的母线 ($bus_id)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.index == bus_id, case.busesDC)
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 连接到不存在的母线 ($bus_id)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 提取光伏阵列参数
                numpanelseries = safe_get_value(row[:numpanelseries], 0, Int)
                numpanelparallel = safe_get_value(row[:numpanelparallel], 0, Int)
                vmpp = safe_get_value(row[:vmpp], 0.0, Float64)
                impp = safe_get_value(row[:impp], 0.0, Float64)
                voc = safe_get_value(row[:voc], 0.0, Float64)
                isc = safe_get_value(row[:isc], 0.0, Float64)
                pvanumcells = safe_get_value(row[:pvanumcells], 0, Int)
                
                # 验证参数值是否有效
                if numpanelseries <= 0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的串联面板数量无效 ($numpanelseries)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if numpanelparallel <= 0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的并联面板数量无效 ($numpanelparallel)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if vmpp <= 0.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的最大功率点电压无效 ($vmpp V)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if impp <= 0.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的最大功率点电流无效 ($impp A)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if voc <= 0.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的开路电压无效 ($voc V)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if isc <= 0.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的短路电流无效 ($isc A)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if pvanumcells <= 0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的电池数量无效 ($pvanumcells)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                #温度参数
                temperature = haskey(row, :temperature) ? safe_get_value(row[:temperature], 25.0, Float64) : 25.0
                # 验证温度值是否合理
                if temperature < -40.0 || temperature > 85.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的温度值无效 ($temperature °C)，设置为默认值 25.0 °C"
                    temperature = 25.0
                end

                # 光照强度参数
                irradiance = haskey(row, :irradiance) ? safe_get_value(row[:irradiance], 1000.0, Float64) : 1000.0
                # 验证光照强度值是否合理
                if irradiance < 0.0 || irradiance > 2000.0
                    @warn "行 $i: 光伏阵列 $name (ID: $index) 的光照强度值无效 ($irradiance W/m²)，设置为默认值 1000.0 W/m²"
                    irradiance = 1000.0
                end

                # 额外参数
                α_isc = haskey(row, :aisctemp) ? safe_get_value(row[:aisctemp], 0.0, Float64) : 0.0

                β_voc = haskey(row, :bvoctemp) ? safe_get_value(row[:bvoctemp], 0.0, Float64) : 0.0
                # 运行状态
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 创建PVArray对象并添加到case中
                push!(case.pvarray, PVArray(
                    index, name, bus_id, numpanelseries, numpanelparallel,
                    vmpp, impp, voc, isc, pvanumcells,temperature,irradiance, α_isc, β_voc,
                    in_service
                ))
                
                # 更新索引映射（如果需要）
                # case.pvarray_indices[index] = length(case.pvarrays)
                
                processed_rows += 1
                
            catch e
                @error "处理光伏阵列数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "光伏阵列数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载光伏阵列数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_ac_pv_system!(case::JuliaPowerCase, file_path::String, sheet_name::String)
从Excel文件中加载交流光伏系统(AC PV System)数据并添加到电力系统案例中。
参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含交流光伏系统数据的工作表名称
"""
function load_ac_pv_system!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    try
        # 使用DataFrame处理
        @info "正在读取交流光伏系统数据..."
        df = DataFrame(XLSX.readtable(file_path, sheet_name))
        inverter = DataFrame(XLSX.readtable(file_path, "inverter"))
        inverter = filter(row -> !ismissing(row[:ParentName]), inverter)
        
        # 确保数据不为空
        if isempty(df)
            @info "交流光伏系统表格为空"
            return
        end
        
        # 将列名转换为小写
        rename!(df, lowercase.(names(df)))
        rename!(inverter, lowercase.(names(inverter)))
        
        # 验证必要的列是否存在
        required_columns = [:index, :id, :bus, :numpanelseries, :numpanelparallel, :vmpp, :impp, :voc, :isc, :pvanumcells]
        missing_columns = filter(col -> !(col in Symbol.(lowercase.(names(df)))), required_columns)
        
        if !isempty(missing_columns)
            @warn "交流光伏系统表格缺少必要列: $(join(missing_columns, ", "))"
            return
        end
        
        # 记录处理的行数和错误的行数
        processed_rows = 0
        error_rows = 0
        
        # 遍历每一行数据
        for (i, row) in enumerate(eachrow(df))
            try
                # 从行数据中提取字段值
                index = safe_get_value(row[:index], 0, Int)
                
                # 验证索引是否有效
                if index <= 0
                    @warn "行 $i: 无效的交流光伏系统索引 ($index)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                name = safe_get_value(row[:id], "", String)
                
                inverterinclude = haskey(row, :inverterincluded) ? safe_get_value(row[:inverterincluded], "", String) : "0"
                if  inverterinclude == "0"
                    # @warn "行 $i: 交流光伏系统 $name (ID: $index) 的光伏阵列不包含逆变器，跳过此行"
                    continue
                end

                # 从母线名称映射到母线ID
                bus_name = safe_get_value(row[:bus], "", String)
                
                # 使用case.bus_name_to_id字典将母线名称转换为整数ID
                bus_id = 0
                
                if haskey(case.bus_name_to_id, bus_name)
                    bus_id = case.bus_name_to_id[bus_name]
                else
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的母线名称 '$bus_name' 在bus_name_to_id字典中不存在，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 验证母线索引是否有效
                if bus_id <= 0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 连接到无效的母线 ($bus_id)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 检查母线是否存在
                if !any(b -> b.index == bus_id, case.busesAC)
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 连接到不存在的母线 ($bus_id)，跳过此行"
                    error_rows += 1
                    continue
                end

                # 提取对应的逆变器
                inverter_row = inverter[inverter.parentname .== name, :][1,:]
                
                # 提取功率参数
                p_mw = safe_get_value(inverter_row[:gencat0ackw], 0.0, Float64)/1000
                q_mvar = (haskey(inverter_row, :gencat0kvar) ? safe_get_value(inverter_row[:gencat0kvar], 0.0, Float64) : 0.0)/1000
                vm_ac_pu = haskey(inverter_row, :vm_ac_pu) ? safe_get_value(inverter_row[:vm_ac_pu], 1.0, Float64) : 1.0
                vm_dc_pu = haskey(inverter_row, :vm_dc_pu) ? safe_get_value(inverter_row[:vm_dc_pu], 1.0, Float64) : 1.0

                # 验证电压值是否合理
                if vm_ac_pu <= 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的交流侧电压无效 ($vm_ac_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_ac_pu = 1.0
                end
                
                if vm_dc_pu <= 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的直流侧电压无效 ($vm_dc_pu p.u.)，设置为默认值 1.0 p.u."
                    vm_dc_pu = 1.0
                end
                
                loss_percent = 1-(haskey(inverter_row, :dcpercenteff) ? safe_get_value(inverter_row[:dcpercenteff], 0.0, Float64) : 0.0)/100
                loss_mw = haskey(inverter_row, :loss_mw) ? safe_get_value(inverter_row[:loss_mw], 0.0, Float64) : 0.0
                
                # 验证损耗是否合理
                if loss_percent < 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的损耗百分比无效 ($loss_percent%)，设置为 0%"
                    loss_percent = 0.0
                end
                
                if loss_mw < 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的损耗功率无效 ($loss_mw MW)，设置为 0 MW"
                    loss_mw = 0.0
                end
                
                max_p_mw = (haskey(inverter_row, :kwmax) ? safe_get_value(inverter_row[:kwmax], 0.0, Float64) : 0.0)/1000
                min_p_mw = (haskey(inverter_row, :kwmin) ? safe_get_value(inverter_row[:kwmin], 0.0, Float64) : 0.0)/1000
                max_q_mvar = (haskey(inverter_row, :kvarmax) ? safe_get_value(inverter_row[:kvarmax], 0.0, Float64) : 0.0)/1000
                min_q_mvar = (haskey(inverter_row, :min_q_mvar) ? safe_get_value(inverter_row[:min_q_mvar], 0.0, Float64) : 0.0)/1000
                
                # 验证功率限制是否合理
                if min_p_mw > max_p_mw && max_p_mw != 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的有功功率限制无效 (min: $min_p_mw MW, max: $max_p_mw MW)，交换值"
                    min_p_mw, max_p_mw = max_p_mw, min_p_mw
                end
                
                if min_q_mvar > max_q_mvar && max_q_mvar != 0.0
                    @warn "行 $i: 换流器 $name (ID: $index) 的无功功率限制无效 (min: $min_q_mvar Mvar, max: $max_q_mvar Mvar)，交换值"
                    min_q_mvar, max_q_mvar = max_q_mvar, min_q_mvar
                end
                
                # 提取光伏阵列参数
                numpanelseries = safe_get_value(row[:numpanelseries], 0, Int)
                numpanelparallel = safe_get_value(row[:numpanelparallel], 0, Int)
                vmpp = safe_get_value(row[:vmpp], 0.0, Float64)
                impp = safe_get_value(row[:impp], 0.0, Float64)
                voc = safe_get_value(row[:voc], 0.0, Float64)
                isc = safe_get_value(row[:isc], 0.0, Float64)
                pvanumcells = safe_get_value(row[:pvanumcells], 0, Int)
                
                # 验证光伏参数值是否有效
                if numpanelseries <= 0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的串联面板数量无效 ($numpanelseries)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if numpanelparallel <= 0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的并联面板数量无效 ($numpanelparallel)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if vmpp <= 0.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的最大功率点电压无效 ($vmpp V)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if impp <= 0.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的最大功率点电流无效 ($impp A)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if voc <= 0.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的开路电压无效 ($voc V)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if isc <= 0.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的短路电流无效 ($isc A)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                if pvanumcells <= 0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的电池数量无效 ($pvanumcells)，跳过此行"
                    error_rows += 1
                    continue
                end
                
                # 环境参数
                irradiance = haskey(row, :irradiance) ? safe_get_value(row[:irradiance], 1000.0, Float64) : 1000.0
                # 验证光照强度值是否合理
                if irradiance < 0.0 || irradiance > 2000.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的光照强度值无效 ($irradiance W/m²)，设置为默认值 1000.0 W/m²"
                    irradiance = 1000.0
                end
                
                temperature = haskey(row, :temperature) ? safe_get_value(row[:temperature], 25.0, Float64) : 25.0
                # 验证温度值是否合理
                if temperature < -40.0 || temperature > 85.0
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的温度值无效 ($temperature °C)，设置为默认值 25.0 °C"
                    temperature = 25.0
                end
                
                # 温度系数参数
                α_isc = haskey(row, :aisctemp) ? safe_get_value(row[:aisctemp], 0.0, Float64) : 0.0
                β_voc = haskey(row, :bvoctemp) ? safe_get_value(row[:bvoctemp], 0.0, Float64) : 0.0
                
                # 控制参数
                control_mode = haskey(inverter_row, :acoperationmode) ? safe_get_value(inverter_row[:acoperationmode], "", String) : "Voltage Control"
                # 验证控制模式
                valid_control_modes = ["Voltage Control", "Swing", "Mvar Control", "PF Control"]
                if !(control_mode in valid_control_modes)
                    @warn "行 $i: 交流光伏系统 $name (ID: $index) 的控制模式无效 ($control_mode)，设置为默认值 'Voltage Control'"
                    control_mode = "Voltage Control"
                end
                
                controllable = haskey(row, :controllable) ? parse_bool(safe_get_value(row[:controllable], true)) : true
                
                # 运行状态
                in_service = haskey(row, :inservice) ? parse_bool(safe_get_value(row[:inservice], true)) : true
                
                # 创建ACPVSystem对象并添加到case中
                push!(case.ACPVSystems, ACPVSystem(
                    index, name, bus_id, p_mw, q_mvar, vm_ac_pu, vm_dc_pu,loss_percent, loss_mw,
                    max_p_mw, min_p_mw, max_q_mvar, min_q_mvar,
                    numpanelseries, numpanelparallel, vmpp, impp, voc, isc, pvanumcells,
                    irradiance, temperature, α_isc, β_voc,
                    control_mode, controllable, in_service
                ))
                
                processed_rows += 1
                
            catch e
                @error "处理交流光伏系统数据第 $i 行时出错" exception=(e, catch_backtrace()) row_data=row
                error_rows += 1
            end
        end
        
        @info "交流光伏系统数据加载完成: 成功处理 $processed_rows 行，错误 $error_rows 行"
        
    catch e
        @error "加载交流光伏系统数据时出错" exception=(e, catch_backtrace())
        rethrow(e)
    end
end


"""
    load_vpps!(case::JuliaPowerCase, file_path::String, sheet_name::String)

从Excel文件中加载虚拟电厂(VPP)数据并添加到电力系统案例中。

参数:
- `case::JuliaPowerCase`: 电力系统案例
- `file_path::String`: Excel文件路径
- `sheet_name::String`: 包含虚拟电厂数据的工作表名称
"""
function load_vpps!(case::JuliaPowerCase, file_path::String, sheet_name::String)
    # 使用DataFrame处理
    df = DataFrame(XLSX.readtable(file_path, sheet_name))
    
    # 确保数据不为空
    if isempty(df)
        @info "虚拟电厂表格为空"
        return
    end
    
    # 将列名转换为小写
    rename!(df, lowercase.(names(df)))
    
    # 遍历每一行数据
    for row in eachrow(df)
        try
            # 从行数据中提取基本字段值
            index = safe_get_value(row[:index], 0, Int)
            name = safe_get_value(row[:name], "", String)
            description = haskey(row, :description) ? safe_get_value(row[:description], "", String) : ""
            control_area = haskey(row, :control_area) ? safe_get_value(row[:control_area], "", String) : ""
            
            # 容量和能量参数
            capacity_mw = haskey(row, :capacity_mw) ? safe_get_value(row[:capacity_mw], 0.0, Float64) : 0.0
            energy_mwh = haskey(row, :energy_mwh) ? safe_get_value(row[:energy_mwh], 0.0, Float64) : 0.0
            
            # 响应和爬坡参数
            response_time_s = haskey(row, :response_time_s) ? safe_get_value(row[:response_time_s], 0.0, Float64) : 0.0
            ramp_rate_mw_per_min = haskey(row, :ramp_rate_mw_per_min) ? safe_get_value(row[:ramp_rate_mw_per_min], 0.0, Float64) : 0.0
            availability_percent = haskey(row, :availability_percent) ? safe_get_value(row[:availability_percent], 100.0, Float64) : 100.0
            
            # 运营信息
            operator = haskey(row, :operator) ? safe_get_value(row[:operator], "", String) : ""
            in_service = haskey(row, :in_service) ? parse_bool(safe_get_value(row[:in_service], true)) : true
            
            # 收集额外的kwargs参数
            kwargs = Dict{Symbol, Any}()
            
            # 资源信息
            if haskey(row, :resource_type)
                kwargs[:resource_type] = safe_get_value(row[:resource_type], "", String)
            end
            
            if haskey(row, :resource_id)
                kwargs[:resource_id] = safe_get_value(row[:resource_id], 0, Int)
            end
            
            if haskey(row, :capacity_share_percent)
                kwargs[:capacity_share_percent] = safe_get_value(row[:capacity_share_percent], 0.0, Float64)
            end
            
            if haskey(row, :control_priority)
                kwargs[:control_priority] = safe_get_value(row[:control_priority], 0, Int)
            end
            
            if haskey(row, :resource_response_time_s)
                kwargs[:resource_response_time_s] = safe_get_value(row[:resource_response_time_s], 0.0, Float64)
            end
            
            if haskey(row, :max_duration_h)
                kwargs[:max_duration_h] = safe_get_value(row[:max_duration_h], 0.0, Float64)
            end
            
            # 负荷信息
            if haskey(row, :timestamp) && !ismissing(row[:timestamp])
                # 处理时间戳，根据实际格式调整
                if isa(row[:timestamp], String)
                    kwargs[:timestamp] = DateTime(row[:timestamp], dateformat"yyyy-mm-dd HH:MM:SS")
                elseif isa(row[:timestamp], DateTime)
                    kwargs[:timestamp] = row[:timestamp]
                else
                    kwargs[:timestamp] = DateTime(now())
                end
            else
                kwargs[:timestamp] = DateTime(now())
            end
            
            if haskey(row, :p_mw)
                kwargs[:p_mw] = safe_get_value(row[:p_mw], 0.0, Float64)
            end
            
            if haskey(row, :q_mvar)
                kwargs[:q_mvar] = safe_get_value(row[:q_mvar], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_up_mw)
                kwargs[:flexibility_up_mw] = safe_get_value(row[:flexibility_up_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_down_mw)
                kwargs[:flexibility_down_mw] = safe_get_value(row[:flexibility_down_mw], 0.0, Float64)
            end
            
            if haskey(row, :flexibility_duration_h)
                kwargs[:flexibility_duration_h] = safe_get_value(row[:flexibility_duration_h], 0.0, Float64)
            end
            
            # 创建VirtualPowerPlant对象并添加到case中
            vpp = VirtualPowerPlant(index, name, description, control_area, capacity_mw, energy_mwh,
                                   response_time_s, ramp_rate_mw_per_min, availability_percent,
                                   operator, in_service; kwargs...)
            
            push!(case.vpps, vpp)
            
            # 更新索引映射
            case.vpp_indices[index] = length(case.vpps)
            
        catch e
            @warn "处理虚拟电厂时出错: $e"
            @warn "问题行: $(row)"
        end
    end
    
    @info "已加载 $(length(case.vpps)) 个虚拟电厂"
end
