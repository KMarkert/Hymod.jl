push!(LOAD_PATH,"../src/")

using Documenter, Hymod

pages = [
    "Home" => "index.md",
    "API" => "api.md",
]

makedocs(;
    modules = [Hymod],
    authors = "Kel Markert",
    repo = "https://github.com/KMarkert/Hymod.jl/blob/{commit}{path}#L{line}",
    sitename = "Hymod.jl",
    # format = Documenter.HTML(;
    #     prettyurls = get(ENV, "CI", "false") == "true",
    #     canonical = "https://deltares.github.io/Wflow.jl",
    #     assets = String[],
    # ),
    pages = pages,
)

deploydocs(;
    repo = "github.com/KMarkert/Hymod.jl.git",
    devbranch = "master"
)
