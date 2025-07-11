using Documenter
using TSPflow
# using TSPflow.PowerFlow
# using TSPflow.Utils
# using TSPflow.TimeSeriesPowerFlow
# 设置 DocMeta
DocMeta.setdocmeta!(TSPflow, :DocTestSetup, :(using TSPflow); recursive=true)

# 生成文档
makedocs(
    sitename = "TSPflow.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/Luosipeng/TSPflow.jl.git",
        canonical = "https://Luosipeng.github.io/TSPflow.jl",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,
    ),
    modules = [TSPflow],
    authors = "Luosipeng",
    repo = "github.com/Luosipeng/TSPflow.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Module" => [
            "Overview" => "modules/overview.md",
            "ComponentModel" => "modules/componentmodel.md",
            "Utils" => "modules/utils.md",
            "PowerFlow" => "modules/powerflow.md",
            "TimeSeriesPowerFlow" => "modules/timeseriespowerflow.md",
        ],
        "API" => [
            "ComponentModel API" => "api/componentmodel.md",
            "Utils API" => "api/utils.md",
            "PowerFlow API" => "api/powerflow.md",
            "TimeSeriesPowerFlow API" => "api/timeseriespowerflow.md",
            "TSPflow API" => "api/tspflow.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs, :cross_references], 
)

# 部署文档
deploydocs(
    repo = "github.com/Luosipeng/TSPflow.jl.git",
    devbranch = "master",
)
