using CSV, DataFrames

include("../src/hapke/functions.jl")

datadir = joinpath(@__DIR__, "data")

# Data obtained from https://doi.org/10.1007/s12648-023-02604-3
ground_truth = CSV.read(joinpath(datadir, "h-function.csv"), DataFrame)

@testset "chandrasekhar h function tests" begin
    @testset "level 1" begin
        h = h_function(ground_truth.mu, ground_truth.w, 1)
        @test isapprox(h, ground_truth.h; rtol=5e-2)
    end

    @testset "level 2" begin
        h = h_function(ground_truth.mu, ground_truth.w, 2)
        @test isapprox(h, ground_truth.h; rtol=8e-3)
    end
end
