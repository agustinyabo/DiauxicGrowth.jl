# to run the documentation generation:
# julia --project=. docs/make.jl
pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, @__DIR__)

using Documenter

# For reproducibility
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force=true,
)
cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force=true,
)

repo_url = "github.com/agustinyabo/DiauxicGrowth.jl"

makedocs(;
    draft=false, # if draft is true, then the julia code from .md is not executed
    # to disable the draft mode in a specific markdown file, use the following:
    #=
    ```@meta
    Draft = false
    ```
    =#
    remotes=nothing,
    warnonly=:cross_references,
    sitename="Diauxic Growth",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["index.md", "implementation.md"],
        example_size_threshold=nothing,
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Introduction" => "index.md", "Implementation and Results" => "implementation.md"
    ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
