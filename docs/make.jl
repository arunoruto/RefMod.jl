push!(LOAD_PATH, "../src/")

using Documenter, RefMod

makedocs(sitename = "RefMod")
deploydocs(repo = "github.com/arunoruto/refmod.jl.git")