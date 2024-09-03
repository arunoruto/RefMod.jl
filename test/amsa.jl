@testitem "AMSA Hopper" begin
    using FITSIO

    # include("../src/hapke/amsa.jl")
    # include("../src/gradient.jl")

    data = FITS(joinpath(@__DIR__, "data", "hopper_amsa.fits"))
    header = read_header(data["result"])

    albedo = read(data["albedo"])

    i = repeat(reshape([sind(header["I"]), 0.0, cosd(header["I"])], 1, 1, 3), size(albedo)...)
    e = repeat(reshape([sind(header["E"]), 0.0, cosd(header["E"])], 1, 1, 3), size(albedo)...)

    hs = header["HS"]
    bs0 = header["BS0"]
    tb = header["TB"]
    hc = header["HC"]
    bc0 = header["BC0"]

    n = normals(read(data["dtm"]), read_header(data["dtm"])["RES"], true)

    phase_function = Dict("b" => header["B"], "c" => header["C"])

    refl = amsa(albedo, i, e, n, phase_function, tb, hs, bs0, hc, bc0)
    result = read(data["result"])

    @test isapprox(refl, result)
    @test_throws ArgumentError amsa(albedo, i, e, n, Dict("a" => header["B"], "d" => header["C"]), tb, hs, bs0, hc, bc0)
end

