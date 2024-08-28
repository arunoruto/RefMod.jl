function process_angle(vec_a::AbstractArray{Float64}, vec_b::AbstractArray{Float64})
    last = ndims(vec_a)
    cos_t = sum(vec_a .* vec_b, dims=last)
    clamp!(cos_t, -1, 1)
    sin_t = sqrt.(1 .- cos_t .^ 2)
    cot_t = cos_t ./ sin_t
    theta = acos.(cos_t)

    # mask = any(isnan.(vec_a) .| isnan.(vec_b), dims=last)
    # cos_t = ifelse.(mask, 0.0, cos_t)
    # sin_t = ifelse.(mask, 0.0, sin_t)
    # cot_t = ifelse.(mask, 0.0, cot_t)
    # theta = ifelse.(mask, 0.0, theta)

    return cos_t, sin_t, cot_t, theta
end

function microscopic_roughness(
    roughness::Float64,
    incidence_direction::AbstractArray{Float64},
    emission_direction::AbstractArray{Float64},
    surface_orientation::AbstractArray{Float64},
)
    sum_abs(x) = sum(abs, x[.!isnan.(x)])
    along_dim = ndims(surface_orientation)

    # Normalize
    surface_orientation ./= sqrt.(sum(abs2, surface_orientation, dims=along_dim))
    incidence_direction ./= sqrt.(sum(abs2, incidence_direction, dims=along_dim))
    emission_direction ./= sqrt.(sum(abs2, emission_direction, dims=along_dim))

    # Angles
    cos_i, sin_i, cot_i, i = process_angle(incidence_direction, surface_orientation)
    cos_e, sin_e, cot_e, e = process_angle(emission_direction, surface_orientation)
    cos_g, _, _, _ = process_angle(incidence_direction, emission_direction)

    # Projections
    # projection_incidence = incidence_direction .- cos_i .* surface_orientation
    # projection_emission = emission_direction .- cos_e .* surface_orientation

    # projection_incidence_length = sqrt.(sum(abs2, projection_incidence, dims=along_dim))
    # projection_emission_length = sqrt.(sum(abs2, projection_emission, dims=along_dim))

    # projection_incidence_length[projection_incidence_length .== 0.0] .= 1.0
    # projection_emission_length[projection_emission_length .== 0.0] .= 1.0

    # projection_incidence ./= projection_incidence_length
    # projection_emission ./= projection_emission_length

    # Azimuth angle
    # cos_psi, sin_psi, _, psi = process_angle(projection_incidence, projection_emission)
    # sin_psi_div_2_sq = abs.(1 .- cos_psi) ./ 2
    cos_psi = (cos_g .- cos_i .* cos_e) ./ (sin_i .* sin_e)
    clamp!(cos_psi, -1, 1)
    sin_psi = sqrt.(1 .- cos_psi .^ 2)
    psi = acos.(cos_psi)
    sin_psi_div_2_sq = abs.(1 .- cos_psi) ./ 2

    # Roughness
    cot_rough = cot(roughness)

    # Masks
    ile = i .< e
    ige = i .>= e
    mask = (cos_i .== 1.0) .| (cos_e .== 1.0)

    f_exp(x, y) = exp.(-x .* 2 * y / pi)
    f_exp_2(x, y) = exp.(-x .^ 2 .* y^2 / pi)
    prime_term(cos_x, sin_x, cot_r, cos_psi, sin_psi_div_2_sq, psi, cot_a, cot_b, cot_c, cot_d, index) = ifelse.(index, cos_x .+ sin_x ./ cot_r .* (cos_psi .* f_exp_2(cot_a, cot_r) .+ sin_psi_div_2_sq .* f_exp_2(cot_b, cot_r)) ./ (2 .- f_exp(cot_c, cot_r) .- psi ./ pi .* f_exp(cot_d, cot_r)), 0.0)
    prime_zero_term(cos_x, sin_x, cot_x, cot_r) = cos_x .+ sin_x ./ cot_r .* f_exp_2(cot_x, cot_r) ./ (2 .- f_exp(cot_x, cot_r))

    factor = 1 / sqrt(1 + pi / cot_rough^2)
    f_psi = exp.(-2 .* sin_psi ./ (1 .+ cos_psi))

    cos_i_s0 = factor .* prime_zero_term(cos_i, sin_i, cot_i, cot_rough)
    cos_e_s0 = factor .* prime_zero_term(cos_e, sin_e, cot_e, cot_rough)

    cos_i_s = zeros(size(cos_i))
    cos_i_s += prime_term(cos_i, sin_i, cot_rough, cos_psi, sin_psi_div_2_sq, psi, cot_e, cot_i, cot_e, cot_i, ile)
    cos_i_s += prime_term(cos_i, sin_i, cot_rough, ones(size(cos_psi)), -sin_psi_div_2_sq, psi, cot_i, cot_e, cot_i, cot_e, ige)
    cos_i_s .*= factor
    cos_i_s = ifelse.(mask, cos_i, cos_i_s)

    cos_e_s = zeros(size(cos_e))
    cos_e_s += prime_term(cos_e, sin_e, cot_rough, ones(size(cos_psi)), -sin_psi_div_2_sq, psi, cot_e, cot_i, cot_e, cot_i, ile)
    cos_e_s += prime_term(cos_e, sin_e, cot_rough, cos_psi, sin_psi_div_2_sq, psi, cot_i, cot_e, cot_i, cot_e, ige)
    cos_e_s .*= factor
    cos_e_s = ifelse.(mask, cos_e, cos_e_s)

    s = factor .* (cos_e_s ./ cos_e_s0) .* (cos_i ./ cos_i_s0)
    s ./= ifelse.(
        ile,
        (1 .- f_psi .+ f_psi .* factor .* cos_i ./ cos_i_s0),
        (1 .- f_psi .+ f_psi .* factor .* cos_e ./ cos_e_s0),
    )
    s = ifelse.(mask, 1.0, s)

    return s, cos_i_s, cos_e_s
end
