using Documenter

# 清理 LOAD_PATH
filter!(x -> x ∉ ["../src/", ".."], LOAD_PATH)

# 重新添加源代码路径
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "..")

# 导入所有需要文档化的模块
using ComponentModel
using Utils
using PowerFlow
using TimeSeriesPowerFlow

# 生成文档
makedocs(
    clean = true,
    sitename = "电力系统建模与分析工具集",
    format = Documenter.HTML(
        prettyurls = false,
        # 移除重复的highlights配置，使用默认值
        assets = ["assets/custom.css"],  # 如果favicon.ico不存在，先移除它
        sidebar_sitename = false,
        footer = "电力系统建模与分析工具集"
    ),
    authors = "开发团队",
    modules = [ComponentModel, Utils, PowerFlow, TimeSeriesPowerFlow],
    pages = [
        "首页" => "index.md",
        "入门指南" => "getting_started.md",
        "API参考" => [
            "概述" => "api/overview.md",
            "ComponentModel模块" => "api/component_model.md",
            "PowerFlow模块" => "api/power_flow.md",
            "Utils模块" => "api/utils.md",
            "TimeSeriesPowerFlow模块" => "api/time_series_power_flow.md"
        ],
        "函数说明" => [
            "组件函数" => "functions/component_functions.md",
            "潮流计算函数" => "functions/power_flow_functions.md",
            "工具函数" => "functions/utils_functions.md",
            "时序潮流函数" => "functions/time_series_functions.md"
        ],
        "主程序说明" => "main.md",
        "示例" => [
            "基础示例" => "examples/basic.md",
            "潮流计算示例" => "examples/power_flow.md",
            "时序潮流示例" => "examples/time_series.md",
            "可视化示例" => "examples/visualization.md"
        ],
        "开发指南" => "development.md",
        "常见问题" => "faq.md"
    ],
    source = "src",
    build = "build",
    warnonly = true
)
