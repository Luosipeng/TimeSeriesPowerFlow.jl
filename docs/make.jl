using Documenter
using TimeSeriesPowerFlow

# 设置 DocMeta
DocMeta.setdocmeta!(TimeSeriesPowerFlow, :DocTestSetup, :(using TimeSeriesPowerFlow); recursive=true)

# 生成文档
makedocs(
    sitename = "TimeSeriesPowerFlow.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/Luosipeng/TimeSeriesPowerFlow.jl.git",
        canonical = "https://Luosipeng.github.io/TimeSeriesPowerFlow.jl",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,
    ),
    modules = [TimeSeriesPowerFlow],
    authors = "Luosipeng",
    repo = "github.com/Luosipeng/TimeSeriesPowerFlow.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Module" => [
            "Overview" => "modules/overview.md",
            "ComponentModel" => "modules/componentmodel.md",
            "Utils" => "modules/utils.md",
            "PowerFlow" => "modules/powerflow.md",
            "TimeDomainPowerFlow" => "modules/timedomainpowerflow.md",
        ],
        "API" => [
            "ComponentModel API" => "api/componentmodel.md",
            "Utils API" => "api/utils.md",
            "PowerFlow API" => "api/powerflow.md",
            "TimeDomainPowerFlow API" => "api/timedomainpowerflow.md",
            "TimeSeriesPowerFlow API" => "api/timeseriespowerflow.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs, :cross_references], 
)

# 部署文档
deploydocs(
    repo = "github.com/Luosipeng/TimeSeriesPowerFlow.jl.git",
    devbranch = "master",
)
