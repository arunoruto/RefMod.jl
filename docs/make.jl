push!(LOAD_PATH, "../src/")

using Documenter, DocumenterCitations, RefMod

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
makedocs(sitename = "RefMod", plugins = [bib])
deploydocs(repo = "github.com/arunoruto/refmod.jl.git")
