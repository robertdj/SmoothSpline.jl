macro swap(x, y)
   quote
      local tmp = $(esc(x))
      $(esc(x)) = $(esc(y))
      $(esc(y)) = tmp
    end
end


struct BSpline
    Knots::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    P::Int64
    M::Int64
    N::Int64

    function BSpline(knots, p = 3)
        sorted_knots = sort(knots)
        vals, lens = StatsBase.rle(sorted_knots)

        lens[1] = p + 1
        lens[end] = p + 1
        extended_knots = StatsBase.inverse_rle(vals, lens)

        m = length(extended_knots) - 1
        n = m - p - 1

        reindexed_knots = OffsetArrays.OffsetVector(extended_knots, 0:length(extended_knots) - 1)

        new(reindexed_knots, p, m, n)
    end
end


function max_knot_index(B::BSpline)
    B.M
end


function max_function_index(B::BSpline)
    B.N
end


# searchsortedlast
function find_span(u, B)
    m = max_knot_index(B)
    if u < B.Knots[0] || u > B.Knots[m]
        throw(ArgumentError("u is not in support of splines"))
    end

    n = max_function_index(B)
    if u == B.Knots[n + 1]
        return n
    end

    p = B.P
    low = p
    high = n + 1
    mid = div(low + high, 2)

    while u < B.Knots[mid] || u >= B.Knots[mid + 1]
        if u < B.Knots[mid]
            high = mid
        else
            low = mid
        end

        mid = div(low + high, 2)
    end

    return mid
end


function basis_funs(i, u, B)
    p = B.P
    N = OffsetArrays.OffsetVector(zeros(p + 1), 0:p)
    N[0] = 1.0

    # TODO: Does initialization matter for left and right?
    left = OffsetArrays.OffsetVector(zeros(p + 1), 0:p)
    right = OffsetArrays.OffsetVector(zeros(p + 1), 0:p)

    for j in 1:p
        left[j] = u - B.Knots[i + 1 - j]
        right[j] = B.Knots[i + j] - u

        saved = 0.0
        for r in 0:j - 1
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        end

        N[j] = saved
    end

    return N
end


function basis_funs_deriv(i, u, n, B)
    p = B.P
    ndu = OffsetArrays.OffsetMatrix(zeros(p + 1, p + 1), 0:p, 0:p)
    ndu[0, 0] = 1.0

    left = Vector{Float64}(undef, p)
    right = similar(left)

    for j in 1:p
        left[j] = u - B.Knots[i + 1 - j]
        right[j] = B.Knots[i + j] - u

        saved = 0.0
        for r in 0:j - 1
            ndu[j, r] = right[r + 1] + left[j - r]
            temp = ndu[r, j - 1] / ndu[j, r]

            ndu[r, j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        end

        ndu[j, j] = saved
    end

    ders = OffsetArrays.OffsetMatrix(zeros(n + 1, p + 1), 0:n, 0:p)
    ders[0, :] = ndu[:, p]

    a = OffsetArrays.OffsetMatrix(zeros(2, p + 1), 0:1, 0:p)
    for r in 0:p
        s1 = 0
        s2 = 1
        a[0, 0] = 1.0

        for k in 1:n
            d = 0.0
            rk = r - k
            pk = p - k
            if r >= k
                a[s2, 0] = a[s1, 0] / ndu[pk + 1, rk]
                d = a[s2, 0] * ndu[rk, pk]
            end

            j1 = rk >= -1 ? 1 : -rk
            j2 = r - 1 <= pk ? k - 1 : p - r

            for j in j1:j2
                a[s2, j] = (a[s1, j] - a[s1, j - 1]) / ndu[pk + 1, rk + j]
                d += a[s2, j] * ndu[rk + j, pk]
            end

            if r <= pk 
                a[s2, k] = -a[s1, k - 1] / ndu[pk + 1, r]
                d += a[s2, k] * ndu[r, pk]
            end

            ders[k, r] = d

            # s1, s2 = s2, s1
            @swap(s1, s2)
        end
    end

    r = p
    for k in 1:n
        ders[k, :] *= r
        r *= p - k
    end

    return ders
end


function single_basis_function_deriv(function_idx, u, n, B)
    span_index = find_span(u, B)

    p = B.P
    if function_idx < span_index - p || function_idx > span_index
        return OffsetArrays.OffsetVector(zeros(n + 1), 0:n)
    end

    ders = basis_funs_deriv(span_index, u, n, B)
    function_index_in_ders = function_idx - span_index + p

    return ders[:, function_index_in_ders]
end
