struct SplineRegression
    Coef::Vector{Float64}
    Data::SplineRegData
    xmin::Float64
    xmax::Float64
end

function smooth_spline(x, y, spar)
    spline_data = SplineRegData(x, y)
    coef = regression(spline_data, spar)

    xmin, xmax = extrema(x)

    SplineRegression(coef, spline_data, xmin, xmax)
end


function predict(sr::SplineRegression, x::Float64)
    if x < sr.xmin || x > sr.xmax
        throw(DomainError(x, "Extrapolation not implemented"))
    end

    spline_interval = find_span(x, sr.Data.B)

    coef_indices = spline_interval - 2:spline_interval + sr.Data.B.P - 2
    spline_coef = @view sr.Coef[coef_indices]
    spline_vals = basis_funs(spline_interval, x, sr.Data.B)

    LinearAlgebra.dot(spline_coef, spline_vals)
end


function predict(sr::SplineRegression, x::AbstractVector{Float64})
    [SmoothSpline.predict(sr, x) for x in x]
end


function regression(sr::SplineRegData, spar)
    ridge_normal_matrix, scaled_y = compute_tikhonov_matrix(sr, spar)
    LinearAlgebra.LAPACK.posv!('U', ridge_normal_matrix, scaled_y)

    return scaled_y
end


function tr(A::Matrix{T}, lead, lag) where T
    n = LinearAlgebra.checksquare(A)
    t = zero(T)
    for i in lead:(n - lag)
        t += A[i,i]
    end
    t
end


function compute_tikhonov_matrix(sr::SplineRegData, spar)
    design_matrix = compute_design_matrix(sr)
    weight_matrix = LinearAlgebra.Diagonal(sr.W)
    sigma = compute_gram_matrix(sr)

    normal_matrix = LinearAlgebra.transpose(design_matrix) * weight_matrix * design_matrix

    normal_trace = (@mock tr(normal_matrix, 3, 3))
    sigma_trace = (@mock tr(sigma, 3, 3))
    trace_ratio = normal_trace / sigma_trace
    λ = trace_ratio * 256^(3 * spar - 1)

    ridge_normal_matrix = normal_matrix + λ * sigma
    scaled_y = LinearAlgebra.transpose(design_matrix) * weight_matrix * sr.Y

    return ridge_normal_matrix, scaled_y
end
