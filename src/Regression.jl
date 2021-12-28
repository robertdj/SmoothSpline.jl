function regression(sr::SplineRegData, spar = 0.5)
    ridge_normal_matrix, scaled_y = compute_tikhonov_matrix(sr, spar)
    LinearAlgebra.LAPACK.posv!('U', ridge_normal_matrix, scaled_y)

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

    r = mytr(normal_matrix, 3, 3) / mytr(sigma, 3, 3)
    λ = r * 256^(3 * spar - 1)

    ridge_normal_matrix = normal_matrix + λ * sigma
    scaled_y = LinearAlgebra.transpose(design_matrix) * weight_matrix * sr.Y

    return ridge_normal_matrix, scaled_y
end
