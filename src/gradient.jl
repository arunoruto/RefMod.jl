function gradient(F, h = missing)
    G = zeros((ndims(F), size(F)...))
    for i = 1:ndims(F)
        l = size(F, i)
        grad_i = zeros(size(F))
        selectdim(grad_i, i, 2:l-1) .= 0.5 * (selectdim(F, i, 3:l) - selectdim(F, i, 1:l-2))
        selectdim(grad_i, i, 1) .= selectdim(F, i, 2) - selectdim(F, i, 1)
        selectdim(grad_i, i, l) .= selectdim(F, i, l) - selectdim(F, i, l - 1)
        selectdim(G, 1, i) .= grad_i
    end

    G = reverse(G, dims = 1)

    if !ismissing(h)
        if (ndims(h) == 0) || ((ndims(h) == 1) && (length(h) <= 2))
            G ./= h
        else
            throw(ArgumentError("h needs to be a scalar, vector, or a matrix"))
        end
    end

    return (selectdim(G, 1, i) for i = 1:size(G, 1))
end
