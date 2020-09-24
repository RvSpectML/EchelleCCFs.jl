using EchelleCCFs
using Documenter

makedocs(;
    modules=[EchelleCCFs],
    authors="Eric Ford",
    repo="https://github.com/RvSpectML/EchelleCCFs.jl/blob/{commit}{path}#L{line}",
    sitename="EchelleCCFs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/EchelleCCFs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RvSpectML/EchelleCCFs.jl",
)
