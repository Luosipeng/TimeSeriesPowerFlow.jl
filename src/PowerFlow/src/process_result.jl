"""
Process power flow calculation results and generate a report.
This function merges results, extracts execution time, and creates a MATPOWER format report.

Parameters:
- results: Array containing calculation results and timing information
- isolated: Information about isolated parts of the network
- file_path: Path where the report will be saved
"""
function process_result(results, isolated, file_path)
    merged_result, area = PowerFlow.merge_results(results[1])
    execution_time = results[2]  # Get execution time (seconds)
    PowerFlow.generate_matpower_report(merged_result, area, execution_time, isolated, file_path)
end
