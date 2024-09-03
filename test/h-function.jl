# Data obtained from https://doi.org/10.1007/s12648-023-02604-3


@testitem "level 1" begin
    using CSV, DataFrames
    include("../src/hapke/functions.jl")
    ground_truth = CSV.read(joinpath(@__DIR__, "data", "h-function.csv"), DataFrame)
    h = h_function(ground_truth.mu, ground_truth.w, 1)
    @test isapprox(h, ground_truth.h; rtol=5e-2)
end

@testitem "level 2" begin
    using CSV, DataFrames
    include("../src/hapke/functions.jl")
    ground_truth = CSV.read(joinpath(@__DIR__, "data", "h-function.csv"), DataFrame)
    h = h_function(ground_truth.mu, ground_truth.w, 2)
    @test isapprox(h, ground_truth.h; rtol=8e-3)
end

@testitem "Errors" begin
    using CSV, DataFrames
    include("../src/hapke/functions.jl")
    ground_truth = CSV.read(joinpath(@__DIR__, "data", "h-function.csv"), DataFrame)
    @test_throws ArgumentError h_function(ground_truth.mu, ground_truth.w, 3)
end

