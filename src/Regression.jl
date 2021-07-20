function regression(sr::SplineRegData, spar = 0.5)
    # design_matrix = compute_design_matrix(sr)
    # weight_matrix = Diagonal(sr.W)
    # sigma = compute_ridge_term(sr)

    # normal_matrix = transpose(design_matrix) * weight_matrix * design_matrix

    # r = tr(normal_matrix) / tr(sigma)
    # λ = r * 256^(3 * spar - 1)

    # ridge_normal_matrix = normal_matrix + λ * sigma
    # LinearAlgebra.axpy!(λ, sigma, normal_matrix)

    # ridge_normal_matrix, scaled_y = compute_tikhonov_matrix(sr, spar)
    # LinearAlgebra.LAPACK.posv!('U', ridge_normal_matrix, scaled_y)

    spline_model_r = RCall.rcopy(R"smooth.spline(x = $(sr.X), y = $(sr.Y), w = $(sr.W), spar = $spar, keep.stuff = TRUE)")
    cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
    tikhonov_matrix_r = transpose(cholesky_r) * cholesky_r
    scaled_y = spline_model_r[:auxM][:XWy]

    LinearAlgebra.LAPACK.posv!('U', tikhonov_matrix_r, scaled_y)

    return scaled_y
end


function compute_tikhonov_matrix(sr::SplineRegData, spar)
    design_matrix = compute_design_matrix(sr)
    weight_matrix = LinearAlgebra.Diagonal(sr.W)
    sigma = compute_ridge_term(sr)

    normal_matrix = LinearAlgebra.transpose(design_matrix) * weight_matrix * design_matrix

    r = LinearAlgebra.tr(normal_matrix) / LinearAlgebra.tr(sigma)
    λ = r * 256^(3 * spar - 1)

    ridge_normal_matrix = normal_matrix + λ * sigma
    scaled_y = LinearAlgebra.transpose(design_matrix) * weight_matrix * sr.Y

    return ridge_normal_matrix, scaled_y
end

function vector_to_upper_triangular(x)
    nrows = div(length(x), 4)

    X = zeros(nrows, nrows)
    
    # factorization in R with dpbfa.f. Storage format is peculiar; works backwards
    X[1, 1] = x[4]
    X[1:2, 2] = x[7:8]
    X[1:3, 3] = x[10:12]
    for col in 4:nrows
        idx1 = col - 4 .+ (1:4)
        idx2 = (col - 1) * 4 .+ (1:4)
        X[idx1, col] = x[idx2]
    end

    return X
end