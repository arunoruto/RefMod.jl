using Test
# using FITSIO

# import RefMod

# f = FITS("data/hopper_imsa.fits")
# println(f)
# println(f["result"])

@testset "RefMod Tests" begin
    include("h-function.jl")
    include("phase-function.jl")
    include("coefficients.jl")
    include("roughness.jl")
end
