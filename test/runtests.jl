# include("setup.jl")

using TestItemRunner

@run_package_tests verbose=true
# @testset "RefMod Tests" begin
#     # Building blocks
#     # include("gradient.jl")
#     # include("h-function.jl")
#     # include("phase-function.jl")
#     # include("coefficients.jl")
#     # include("roughness.jl")

#     # Models
#     include("amsa.jl")
# end
