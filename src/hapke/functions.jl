"""
Calculates the Hapke function for a given set of parameters.

# Arguments:
- `x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The input parameter.
- `w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The weight array.
- `level::Int64=1`: The accuracy level of the Ambartsumian–Chandrasekhar function to calculate.
- `derivative::Bool=false`: Whether to calculate the derivative.

# Returns:
- `Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The calculated Hapke function value.

# Raises:
- ArgumentError: If an invalid level is provided.
"""
function h_function(
    x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    level::Int64=1,
    derivative::Bool=false,
)
    if (level < 0) | (level > 2)
        throw(ArgumentError("level should be 1 or 2!"))
    end

    if level == 1
        h_function_1(x, w)
    elseif level == 2
        h_function_2(x, w)
    end
end

"""
Calculates the Ambartsumian–Chandrasekhar function level 1 for a given set of parameters.

# Arguments:
- `x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The input parameter.
- `w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The weight array.

# Returns:
- `Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The calculated Ambartsumian–Chandrasekhar function level 1 value.

# Note:
Equation 2 in [hapke2002chapter5](@cite).
"""
function h_function_1(
    x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}
)

    gamma = @. sqrt(1 - w)
    return @. (1 + 2 * x) / (1 + 2 * x * gamma)
end

"""
Calculates the Ambartsumian–Chandrasekhar function level 2 for a given set of parameters.

# Arguments
- `x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The input parameter.
- `w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The weight array.

# Returns:
- `Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: The calculated Ambartsumian–Chandrasekhar function level 2 value.

# Note:
Equation 13 in [hapke2002chapter5](@cite).
"""
function h_function_2(
    x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}
)
    gamma = @. sqrt(1 - w)
    r0 = @. (1 - gamma) / (1 + gamma)
    x_log_term = @. log(1 + 1 / x)
    h = @. 1 / (1 - w * x * (r0 + (1 - 2 * r0 * x) / 2 * x_log_term))
    @. h[x==0] = 1.0
    return h
end

"""
Calculates the phase function for the double Henyey-Greenstein model.

# Arguments:
- `cos_g::AbstractArray{<:AbstractFloat}`: The cosine of the scattering angle.
- `b::Float64=0.21`: The asymmetry parameter.
- `c::Float64=0.7`: The backscatter fraction.

# Returns:
- `AbstractArray{<:AbstractFloat}`: The phase function value.

# Note:
Default values for `b` and `c` are for the moon, taken from [warell2004](@cite).
"""
function double_henyey_greenstein(
    cos_g::AbstractArray{<:AbstractFloat},
    b::Float64=0.21,
    c::Float64=0.7,
)
    p_g = (1 + c) / 2 * (1 - b^2) ./ (1 .- 2 * b .* cos_g .+ b^2) .^ 1.5 + (1 - c) / 2 * (1 - b^2) ./ (1 .+ 2 * b .* cos_g .+ b^2) .^ 1.5
    return p_g
end
