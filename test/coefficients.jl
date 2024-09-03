@testitem "A coefficients" begin
    include("../src/hapke/legendre.jl")
    a_n = coef_a()
    gt = [0.0, -0.5, 0.0, 0.125, 0.0, -0.0625, 0.0, 0.0390625, 0.0, -0.02734375, 0.0, 0.0205078125, 0.0, -0.01611328125, 0.0, 0.013092041015625]
    @test isapprox(a_n, gt)
end

@testitem "b coefficients" begin
    include("../src/hapke/legendre.jl")
    b_n = coef_b(0.18, 1.1)
    gt = [1.1, 0.594, 0.1782, 0.0449064, 0.010392624, 0.0022863772799999998, 0.0004863748031999999, 0.00010101630527999999, 2.0607326277119994e-5, 4.1457091686911986e-6, 8.247779293501437e-7, 1.6259907750045692e-7, 3.181286298921983e-8, 6.184420565104335e-9, 1.195654642586838e-9, 2.300604450218812e-10]
    @test isapprox(b_n, gt)
end

@testitem "P function" begin
    include("../src/hapke/legendre.jl")
    p = function_p(LinRange(0, 1, 10), 0.18, 1.1)
    gt = [1.0, 0.9660547437620005, 0.9322352894770082, 0.8986648155318951, 0.8654613865764601, 0.8327357177287303, 0.8005892911272132, 0.7691128924542041, 0.73838560141173, 0.7084742372481021]
    @test isapprox(p, gt)
end

@testitem "P value" begin
    include("../src/hapke/legendre.jl")
    p = value_p(0.18, 1.1)
    gt = 1.1492107509695955
    @test isapprox(p, gt)
end
