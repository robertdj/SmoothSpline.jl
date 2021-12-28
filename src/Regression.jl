function regression2(sr::SplineRegData, spar = 0.5)
    ridge_normal_matrix, scaled_y = compute_tikhonov_matrix(sr, spar)
    LinearAlgebra.LAPACK.posv!('U', ridge_normal_matrix, scaled_y)

    return scaled_y
end


function regression(sr::SplineRegData, spar = 0.5)
    spline_model_r = RCall.rcopy(R"smooth.spline(x = $(sr.X), y = $(sr.Y), w = $(sr.W), spar = $spar, keep.stuff = TRUE)")

    # XWX = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
    # @show XWX[LinearAlgebra.diagind(XWX)]
    # Sigma = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])
    # @show Sigma[LinearAlgebra.diagind(Sigma)]

    # @show spline_model_r[:lambda]

    cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
    tikhonov_matrix_r = transpose(cholesky_r) * cholesky_r
    scaled_y = spline_model_r[:auxM][:XWy]

    LinearAlgebra.LAPACK.posv!('U', tikhonov_matrix_r, scaled_y)

    return scaled_y
end

function mytr(A::Matrix{T}, lead, lag) where T
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

    # @show normal_matrix[LinearAlgebra.diagind(normal_matrix)]
    # @show sigma[LinearAlgebra.diagind(sigma)]

    r = mytr(normal_matrix, 3, 3) / mytr(sigma, 3, 3)
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
