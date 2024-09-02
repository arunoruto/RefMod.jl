# using AssociatedLegendrePolynomials

function coef_a(n::Int64=15)
    a_n = zeros(n + 1)
    a_n[2] = -0.5
    for i in 4:2:(n+1)
        a_n[i] = (3 / i - 1) * a_n[i-2]
    end
    return a_n
end

function coef_b(b::Float64=0.21, c::Float64=0.7, n::Int=15)::Array{Float64}
    range = 0:n
    b_n = c .* (2 .* range .+ 1) .* b .^ range
    # Check why this holds...
    # b_n[1] = 1.0
    return b_n
end

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

function value_p(
    b::Float64=0.21,
    c::Float64=0.7,
    n::Int=15
)
    a_n = coef_a(n)
    b_n = coef_b(b, c, n)
    # TODO: Change 1+ to 1-... Maybe?
    return 1 + sum(@. a_n^2 * b_n)
end
