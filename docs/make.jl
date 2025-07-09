using Documenter
using TSPF
# using TSPF.PowerFlow
# using TSPF.Utils
# using TSPF.TimeSeriesPowerFlow
# 设置 DocMeta
DocMeta.setdocmeta!(TSPF, :DocTestSetup, :(using TSPF); recursive=true)

# 生成文档
makedocs(
    sitename = "TSPF.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/Luosipeng/TSPF.jl.git",
        canonical = "https://Luosipeng.github.io/TSPF.jl",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,
    ),
    modules = [TSPF],
    authors = "Luosipeng",
    repo = "github.com/Luosipeng/TSPF.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Module" => [
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
            "TSPF API" => "api/tspf.md",
        ],
    ],
)

# 部署文档
deploydocs(
    repo = "github.com/Luosipeng/TSPF.jl.git",
    devbranch = "master",
)

# using .TSPF
# @doc ComponentModel
# @doc TSPF.ComponentModel.AbstractComponent
