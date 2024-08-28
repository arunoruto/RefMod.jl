using FITSIO

include("../src/hapke/roughness.jl")

datadir = joinpath(@__DIR__, "data")
data = FITS(joinpath(datadir, "roughness.fits"))

@testset "Roughness" begin
    @testset "Roughness test" begin
        i_vec = read(data["incidence"])
        e_vec = read(data["emission"])
        N = size(i_vec)[1]
        normal = repeat(reshape([0.0, 0.0, 1.0], 1, 1, 3), N, N)
        tb = deg2rad(read_header(data["S"])["ROUGHNESS"])
        S, mu0, mu = microscopic_roughness(tb, i_vec, e_vec, normal)
        S = dropdims(S, dims=3)
        mu0 = dropdims(mu0, dims=3)
        mu = dropdims(mu, dims=3)

        @test isapprox(mu0, read(data["mu0"]))
        @test isapprox(mu, read(data["mu"]))
        @test isapprox(S, read(data["S"]))
    end
end
