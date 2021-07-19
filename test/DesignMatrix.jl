using SmoothSpline
using Test

# @testset "Design matrix" begin
#     if Base.find_package("RCall") === nothing
#         # @test_skip
#     end
# end

# import RCall
using RCall
import RDatasets
import LinearAlgebra

mtcars = RDatasets.dataset("datasets", "mtcars")
hp = mtcars[!, :HP]
mpg = mtcars[!, :MPG]

spline_model_r = RCall.rcopy(R"smooth.spline($hp, $mpg, spar = 0.5, keep.stuff = TRUE)")

# spline_data = SplineRegData(hp, mpg)
spline_data = SplineRegData(float.(hp), float.(mpg))


# spline_reg_data = SplineRegData(hp, mpg)
# coef = regression(spline_reg_data, spar = 0.5)

function vector_to_banded_matrix(x)
    nrows = div(length(x), 4)

    X = zeros(nrows, nrows)

    X[LinearAlgebra.diagind(X)] = x[1:nrows]

    X[LinearAlgebra.diagind(X, 1)] = x[nrows + 1:2*nrows - 1]
    X[LinearAlgebra.diagind(X, -1)] = x[nrows + 1:2*nrows - 1]

    X[LinearAlgebra.diagind(X, 2)] = x[2*nrows + 1:3*nrows - 2]
    X[LinearAlgebra.diagind(X, -2)] = x[2*nrows + 1:3*nrows - 2]

    X[LinearAlgebra.diagind(X, 3)] = x[3*nrows + 1:4*nrows - 3]
    X[LinearAlgebra.diagind(X, -3)] = x[3*nrows + 1:4*nrows - 3]

    return X
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


@testset "Design matrix" begin
    design_matrix_julia = SmoothSpline.compute_design_matrix(spline_data)
    W = LinearAlgebra.diagm(spline_data.W)
    weighted_design_matrix = transpose(design_matrix_julia) * W * design_matrix_julia

    weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])

    @test weighted_design_matrix ≈ weighted_design_matrix_r
    @test maximum(abs, weighted_design_matrix - weighted_design_matrix_r) < 1e-14
end


@testset "R's storage of Tikhonov matrix" begin
    weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
    lambda_r = spline_model_r[:lambda]
    gram_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])

    tikhonov_matrix_r = weighted_design_matrix_r + lambda_r * gram_r
    cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
    
    @test transpose(cholesky_r) * cholesky_r ≈ tikhonov_matrix_r
    @test maximum(abs, transpose(cholesky_r) * cholesky_r - tikhonov_matrix_r) < 1e-14
end


@testset "Gram matrix" begin
    sigma_julia = SmoothSpline.compute_ridge_term(spline_data)
    sigma_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])
end


@testset "Cholesky factor" begin
    # lambda_r = spline_model_r[:lambda]
    # rnm = XWX + lambda * sigma
    # c = cholesky(rnm)
    # maximum(abs, c.U - R)
end




# R = zeros(size(sigma))
# # factorization in R with dpbfa.f. Storage format is peculiar; works backwards
# R[1, 1] = ss_r[:auxM][:R][4]
# R[1:2, 2] = ss_r[:auxM][:R][7:8]
# R[1:3, 3] = ss_r[:auxM][:R][10:12]
# for col in 4:24
#     idx1 = col - 4 .+ (1:4)
#     idx2 = (col - 1) * 4 .+ (1:4)
#     R[idx1, col] = ss_r[:auxM][:R][idx2]
# end

# lambda = ss_r[:lambda]

# rnm = XWX + lambda * sigma

# # Very badly conditioned rnm
# cond(rnm)
# c = cholesky(rnm)

# maximum(abs, c.U - R)
# maximum(abs, ridge_normal_matrix - rnm)
# maximum(abs, scaled_y - ss_r[:auxM][:XWy])

# # Like dpbsl.f in R. See also sslvrg.f
# A = deepcopy(rnm)
# B = deepcopy(scaled_y)
# LinearAlgebra.LAPACK.posv!('U', A, B)

# # Alternatively:
# # # Cholesky of A
# # LinearAlgebra.LAPACK.potrf!('U', A)
# # LinearAlgebra.LAPACK.potrs!('U', A, B)

# coef = ss_r[:fit][:coef]
# maximum(abs, coef - B)
