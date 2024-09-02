include("setup.jl")

@testset "RefMod Tests" begin
    # Building blocks
    include("gradient.jl")
    include("h-function.jl")
    include("phase-function.jl")
    include("coefficients.jl")
    include("roughness.jl")

    # Models
    include("amsa.jl")
end
