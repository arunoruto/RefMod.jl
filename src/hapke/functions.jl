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
# Arguments
- `x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: argument of Chandrasekhar's h function.
- `w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}`: albedo parameter of the function.
"""
function h_function_1(
    x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}
)
    gamma = @. sqrt(1 - w)
    return @. (1 + 2 * x) / (1 + 2 * x * gamma)
    # return (1 .+ 2 .* x) ./ (1 .+ 2 .* x .* gamma)
end

function h_function_2(
    x::Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
    w::Union{AbstractFloat,AbstractArray{<:AbstractFloat}}
)
    gamma = sqrt.(1 .- w)
    r0 = (1 .- gamma) ./ (1 .+ gamma)
    x_log_term = log.(1 .+ 1 ./ x)
    h = 1 ./ (1 .- w .* x .* (r0 .+ (1 .- 2 .* r0 .* x) ./ 2 .* x_log_term))
    h[x.==0] .= 1.0
    return h
end

function double_henyey_greenstein(cos_g, b::Float64=0.21, c::Float64=0.7)
    p_g = (1 + c) / 2 * (1 - b^2) ./ (1 .- 2 * b .* cos_g .+ b^2) .^ 1.5 + (1 - c) / 2 * (1 - b^2) ./ (1 .+ 2 * b .* cos_g .+ b^2) .^ 1.5
end
