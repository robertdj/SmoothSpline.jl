struct SplineRegData
    B::BSpline
    X::Vector{Float64}
    Y::Vector{Float64}
    W::Vector{Float64}

    function SplineRegData(B, x, y, w)
        length_x = length(x)
        length_y = length(y)
        length_w = length(w)

        if !(length_x == length_y == length_w)
            throw(ArgumentError("X, Y and W must have the same length"))
        end

        new(B, x, y, w)
        # new(B, float.(x), float.(y), float.(w))
    end
end


function SplineRegData(x, y)
    w = ones(length(x))

    SplineRegData(x, y, w)
end


function SplineRegData(x, y, w)
    sort_perm_for_x = sortperm(x)
    sorted_x = x[sort_perm_for_x]
    sorted_y = y[sort_perm_for_x]

    # TODO: Should also work for almost equal x's
    # TODO: Not for general w
    unique_x, weights = StatsBase.rle(sorted_x)
    unique_x_length = length(unique_x)
    unique_y = zeros(unique_x_length)

    y_index_count = 0
    for x_idx in 1:unique_x_length
        y_mean_for_same_x = 0.0
        same_x_count = 0
        for l in 1:weights[x_idx]
            y_index_count += 1
            same_x_count += 1
            y_mean_for_same_x += (sorted_y[y_index_count] - y_mean_for_same_x) / same_x_count
        end
        unique_y[x_idx] += y_mean_for_same_x
    end

    internal_knots = compute_internal_knots_from_obs(unique_x)
    knots = extend_internal_knots(internal_knots)
    B = BSpline(knots)

    SplineRegData(B, unique_x, unique_y, weights)
end


function compute_internal_knots_from_obs(x)
    length_x = length(x)
    number_of_knots = compute_number_of_knots(length_x)

    x_indicies_for_knots = floor.(Int64, range(1, stop = length_x, length = number_of_knots))

    return x[x_indicies_for_knots]
end


function extend_internal_knots(knots, p = 3)
    vcat(
        repeat(knots[1:1], p),
        knots,
        repeat(knots[end:end], p)
    )
end


function compute_number_of_knots(n)
    if n < 50
        return n
    end

    a1 = log2(50)
    a2 = log2(100)
    a3 = log2(140)
    a4 = log2(200)

    if	n < 200
        m = 2^(a1 + (a2 - a1)*(n - 50)/150)
    elseif n < 800
        m =  2^(a2 + (a3 - a2)*(n - 200)/600)
    elseif n < 3200
        m =  2^(a3 + (a4 - a3)*(n - 800)/2400)
    else  
        m =  200 + (n - 3200)^0.2
    end

    return trunc(Int, m)
end


function compute_design_matrix(sr)
    number_of_functions = max_function_index(sr.B) + 1
    length_x = length(sr.X)

    p = sr.B.P
    design_matrix0 = OffsetArrays.OffsetMatrix(zeros(length_x, number_of_functions), 0:length_x - 1, 0:number_of_functions - 1)
    X = OffsetArrays.OffsetVector(sr.X, 0:length_x - 1)

    for row = 0:length_x - 1
        u = X[row]
        span_index = find_span(u, sr.B)
        spline_value = basis_funs(span_index, u, sr.B)

        col_index = span_index - p .+ (0:p)
        design_matrix0[row, col_index] = spline_value
    end

    design_matrix = design_matrix0.parent

    return design_matrix
end


function compute_gram_matrix(sr)
    N = max_function_index(sr.B)
    sigma = OffsetArrays.OffsetMatrix(zeros(N + 1, N + 1), 0:N, 0:N)

    p = sr.B.P

    for i = 0:N
        j_upper = min(N, i + p + 1)
        for j = i:j_upper
            entry = compute_gram_matrix_entry(i, j, sr.B)
            sigma[i, j] = entry
            sigma[j, i] = entry
        end
    end

    return sigma.parent
end


"""
Dummy function to use different approximations of `1/3`.
"""
OneThird() = 1/3


function compute_gram_matrix_entry(i, j, B)
    p = B.P
    if abs(i - j) > p
        return 0.0
    end

    min_index, max_index = minmax(i, j)
    M = max_knot_index(B)

    one_third = (@mock OneThird())

    result = 0.0

    k_start = min_index
    k_stop = min(max_index + p + 1, M - 1)
    for k in k_start:k_stop
        yi1 = single_basis_function_deriv(i, B.Knots[k], 2, B)[2]
        yi2 = single_basis_function_deriv(i, B.Knots[k + 1], 2, B)[2]
        Δyi = yi2 - yi1

        yj1 = single_basis_function_deriv(j, B.Knots[k], 2, B)[2]
        yj2 = single_basis_function_deriv(j, B.Knots[k + 1], 2, B)[2]
        Δyj = yj2 - yj1

        Δₖ = B.Knots[k + 1] - B.Knots[k]

        kth_term = yi1 * yj1 + 0.5 * (yi1 * Δyj + yj1 * Δyi) + Δyj * Δyi * one_third
        kth_term *= Δₖ

        result += kth_term
    end

    return result
end
