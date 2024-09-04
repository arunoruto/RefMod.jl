"""
    coef_a(n::Int64=15)

Calculates the coefficients 'a_n' for the Legendre polynomial series.

# Arguments
- `n::Int64=15`: The number of coefficients to calculate. Default is 15.

# Returns
- `Array{Float64}`: An array of coefficients 'a_n' for the Legendre polynomial series.

# Note
Equation 27 in [hapke2002chapter5](@cite).
"""
function coef_a(n::Int64=15)
    a_n = zeros(n + 1)
    a_n[2] = -0.5
    for i in 4:2:(n+1)
        a_n[i] = (3 / i - 1) * a_n[i-2]
    end
    return a_n
end

"""
    coef_b(b::Float64=0.21, c::Float64=0.7, n::Int=15)

Calculates the coefficients for the Hapke reflectance model Legendre polynomial expansion.

# Arguments
- `b::Float64=0.21`: The single scattering albedo.
- `c::Float64=0.7`: The asymmetry factor.
- `n::Int=15`: The number of coefficients to calculate.

# Returns
- `Array{Float64}`: The calculated coefficients for the Legendre polynomial expansion.

# Note
Equation on page 530 in [hapke2002chapter5](@cite).
Default values for `b` and `c` are for the moon, taken from [warell2004](@cite).
"""
function coef_b(b::Float64=0.21, c::Float64=0.7, n::Int=15)::Array{Float64}
    range = 0:n
    b_n = c .* (2 .* range .+ 1) .* b .^ range
    # Check why this holds...
    # b_n[1] = 1.0
    return b_n
end

"""
    function_p(x::AbstractArray{<:AbstractFloat}, b::Float64=0.21, c::Float64=0.7, n::Int=15)

Calculates the P function using the Hapke reflectance model.

# Arguments
- `x::AbstractArray{<:AbstractFloat}`: The input array.
- `b::Float64=0.21`: The single scattering albedo.
- `c::Float64=0.7`: The asymmetry factor.
- `n::Int=15`: The number of coefficients to calculate.

# Returns
- `Array{Float64}`: The calculated P function.

# Note
Equations 23 and 24 in [hapke2002chapter5](@cite).
Default values for `b` and `c` are for the moon, taken from [warell2004](@cite).
"""
function function_p(
    x::AbstractArray{<:AbstractFloat},
    b::Float64=0.21,
    c::Float64=0.7,
    n::Int=15
)
    a_n = coef_a(n)
    b_n = coef_b(b, c, n)

    previous_2 = ones(size(x))
    previous_1 = x

    res = a_n[1] * b_n[1] .* previous_2
    res += a_n[2] * b_n[2] .* previous_1

    for i in 2:n
        temp = previous_1
        previous_1 = (2 - 1 / i) .* x .* previous_1 - (1 - 1 / i) .* previous_2
        previous_2 = temp
        res += a_n[i+1] * b_n[i+1] .* previous_1
    end

    return 1 .+ res
end

"""
    value_p(b::Float64=0.21, c::Float64=0.7, n::Int=15)

Calculates the value of the P function.

# Arguments
- `b::Float64=0.21`: The single scattering albedo.
- `c::Float64=0.7`: The asymmetry factor.
- `n::Int=15`: The number of coefficients to calculate.

# Returns
- `Float64`: The calculated P value.

# Note
Equations 25 in [hapke2002chapter5](@cite).
Default values for `b` and `c` are for the moon, taken from [warell2004](@cite).
"""
function value_p(
    b::Float64=0.21,
    c::Float64=0.7,
    n::Int=15
)
    a_n = coef_a(n)
    b_n = coef_b(b, c, n)
    # TODO: Change 1+ to 1-... Maybe?
    return 1 + sum(a_n .^ 2 .* b_n)
end
