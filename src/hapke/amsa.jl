include("roughness.jl")
include("functions.jl")
include("legendre.jl")

"""
Calculates the reflectance using the AMSA (Advanced Modified Shadowing and Coherent Backscattering) model.

# Arguments
- `single_scattering_albedo`: The single scattering albedo.
- `incidence_direction`: The incidence direction.
- `emission_direction`: The emission direction.
- `surface_orientation`: The surface orientation.
- `phase_function`: The phase function.
- `roughness::Float64=0.0`: The roughness.
- `hs::Float64=0.0`: The shadowing parameter.
- `bs0::Float64=0.0`: The shadowing parameter.
- `hc::Float64=0.0`: The coherent backscattering parameter.
- `bc0::Float64=0.0`: The coherent backscattering parameter.

# Returns
- `Array{Float64}`: An array representing the reflectance values.

# Raises
- `ArgumentError`: If the phase function parameters are unsupported.
"""
function amsa(
    single_scattering_albedo,
    incidence_direction,
    emission_direction,
    surface_orientation,
    phase_function,
    roughness::Float64=0.0,
    hs::Float64=0.0,
    bs0::Float64=0.0,
    hc::Float64=0.0,
    bc0::Float64=0.0,
)

    along_dim = ndims(surface_orientation)
    length_along(x) = sqrt.(sum(abs2, x, dims=along_dim))

    # Normalize
    surface_orientation ./= length_along(surface_orientation)
    incidence_direction ./= length_along(incidence_direction)
    emission_direction ./= length_along(emission_direction)

    # Roguhness
    s, mu_0, mu = microscopic_roughness(
        roughness,
        incidence_direction,
        emission_direction,
        surface_orientation,
    )

    # Alpha angle
    cos_alpha, sin_alpha, _, _ = process_angle(incidence_direction, emission_direction)
    tan_alpha_2 = sin_alpha ./ (1 .+ cos_alpha)

    # Legendre
    if haskey(phase_function, "b") && haskey(phase_function, "c")
        # Double Hanyey-Greenstein
        b = phase_function["b"]
        c = phase_function["c"]
        p_g = double_henyey_greenstein(cos_alpha, phase_function["b"], phase_function["c"])

        p_mu_0 = function_p(mu_0, b, c)
        p_mu = function_p(mu, b, c)
        p = value_p(b, c)
    else
        throw(ArgumentError("Unsupported phase function parameters"))
    end

    # H function
    h_mu_0 = h_function_2(mu_0, single_scattering_albedo)
    h_mu = h_function_2(mu, single_scattering_albedo)

    # M term
    m = p_mu_0 .* (h_mu .- 1) + p_mu .* (h_mu_0 .- 1) + p .* (h_mu_0 .- 1) .* (h_mu .- 1)

    # Shadow-hiding effect
    b_sh = ones(size(tan_alpha_2))
    if (bs0 != 0) && (hs != 0)
        @. b_sh += bs0 / (1 + tan_alpha_2 / hs)
    end

    # Coherent backscattering effect
    b_cb = ones(size(tan_alpha_2))
    if (bc0 != 0) && (hc != 0)
        hc_2 = tan_alpha_2 ./ hc
        bc = @. 0.5 * (1 + (1 - exp(-hc_2)) / hc_2) / (1 + hc_2)^2
        b_cb += bc0 .* bc
    end

    albedo_independent = @. mu_0 / (mu_0 + mu) * s / (4 * pi) * b_cb
    refl = @. albedo_independent * single_scattering_albedo * (b_sh * p_g + m)

    return refl
end
