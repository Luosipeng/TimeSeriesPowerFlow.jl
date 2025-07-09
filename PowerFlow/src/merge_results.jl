# 简单合并孤岛的潮流计算结果
function merge_results(results)
    # 创建一个新的 JPC 对象
    merged_result = JPC()
    
    # 合并基本结果字段
    merged_result.success = all([r.success for r in results])
    merged_result.iterationsAC = maximum([r.iterationsAC for r in results])
    
    # 设置 baseMVA（假设所有结果使用相同的基准功率）
    if !isempty(results)
        merged_result.baseMVA = results[1].baseMVA
        merged_result.version = results[1].version
    end
    
    # 合并主要数据矩阵并按第一列排序
    for key in [:busAC, :branchAC, :genAC, :loadAC, :loadAC_flex, :loadAC_asymm, 
                :branch3ph, :busDC, :branchDC, :sgenAC, :storage, :sgenDC, 
                :converter, :ext_grid, :hvcb, :microgrid]
        # 收集所有结果中的相应矩阵
        matrices = []
        for r in results
            if hasproperty(r, key) && !isempty(getproperty(r, key))
                push!(matrices, getproperty(r, key))
            end
        end
        
        if !isempty(matrices)
            # 垂直连接所有矩阵
            combined = vcat(matrices...)
            
            # 按第一列排序（如果矩阵非空）
            if !isempty(combined)
                sorted_indices = sortperm(combined[:, 1])
                setproperty!(merged_result, key, combined[sorted_indices, :])
            else
                setproperty!(merged_result, key, combined)
            end
        end
    end
    
    area = length(results)
    return merged_result, area
end
