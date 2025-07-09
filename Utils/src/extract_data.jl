"""
从Excel工作表中提取数据
    
参数:
- sheet_name::String: 工作表名称
- sheets_data::Dict{String, DataFrame}: 包含所有工作表数据的字典

返回:
- DataFrame: 提取的数据子集
"""
function extract_data(sheet_name::String, sheets_data::Dict{String, DataFrame})::DataFrame

    # 检查工作表是否存在
    if !haskey(sheets_data, sheet_name)
        # @warn "工作表 '$sheet_name' 不存在"
        return DataFrame()
    end
    
    # 获取当前工作表
    current_sheet::DataFrame = sheets_data[sheet_name]
    
    # 初始化计数器
    row_count::Int = 0
    col_count::Int = 0
    
    # 计算有效行数
    sheet_rows::Int = size(current_sheet, 1)
    for row_idx in 1:sheet_rows
        if ismissing(current_sheet[row_idx, 1])
            row_count = row_idx - 1
            break
        end
        row_count = row_idx
    end
    
    # 计算有效列数
    sheet_cols::Int = size(current_sheet, 2)
    for col_idx in 1:sheet_cols
        if ismissing(current_sheet[1, col_idx])
            col_count = col_idx - 1
            break
        end
        col_count = col_idx
    end
    
    # 检查是否存在有效数据
    if row_count == 0 || col_count == 0
        @warn "工作表 '$sheet_name' 中未找到有效数据"
        return DataFrame()
    end
    
    # 返回截取的数据（从第二行开始，跳过表头）
    return current_sheet[2:row_count, 1:col_count]
end