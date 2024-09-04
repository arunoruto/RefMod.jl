push!(LOAD_PATH, "../src/")

using TOML
using Documenter
using DocumenterCitations
using RefMod

PROJECT_TOML = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "github.com/arunoruto/$NAME.jl"

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib"),
    style=:authoryear,
)
makedocs(
    modules = [RefMod],
    authors=AUTHORS,
    sitename="$NAME.jl",
    format=Documenter.HTML(
        prettyurls=true,
        assets=String["assets/citations.css"],
        footer="[$NAME.jl](http://$GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
    ),
    plugins=[bib]
)
deploydocs(repo="$GITHUB.git")
