using Documenter
using HyDistFlow

# 设置 DocMeta
DocMeta.setdocmeta!(HyDistFlow, :DocTestSetup, :(using HyDistFlow); recursive=true)

# 生成文档
makedocs(
    sitename = "HyDistFlow.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/Luosipeng/HyDistFlow.jl.git",
        canonical = "https://Luosipeng.github.io/HyDistFlow.jl",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,
    ),
    modules = [HyDistFlow],
    authors = "Luosipeng",
    repo = "github.com/Luosipeng/HyDistFlow.jl/blob/{commit}{path}#L{line}",
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
            "HyDistFlow API" => "api/hydistflow.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs, :cross_references], 
)

# 部署文档
deploydocs(
    repo = "github.com/Luosipeng/HyDistFlow.jl.git",
    devbranch = "master",
)
