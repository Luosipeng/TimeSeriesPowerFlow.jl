function process_result(results, isolated, file_path)
    merged_result, area = PowerFlow.merge_results(results[1])
    execution_time = results[2]  # 获取执行时间（秒）
    PowerFlow.generate_matpower_report(merged_result,area,execution_time,isolated,file_path)
end