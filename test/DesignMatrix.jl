using SmoothSpline
using Mocking
using Test

using RCall
import QuadGK
import LinearAlgebra

import Random

include("ReshapeROutput.jl")

@testset "Design matrix" begin
    Random.seed!(1)

    x = [0.0; rand(15); 1.0] |> sort
    y = rand(length(x))
    
    spline_model_r = RCall.rcopy(R"smooth.spline($x, $y, spar = 0.5, keep.stuff = TRUE)")
    
    spline_data = SplineRegData(x, y)
    
    @testset "Compare design matrix with R" begin
        design_matrix_julia = SmoothSpline.compute_design_matrix(spline_data)
        W = LinearAlgebra.diagm(spline_data.W)
        weighted_design_matrix = transpose(design_matrix_julia) * W * design_matrix_julia
    
        weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
    
        @test weighted_design_matrix ≈ weighted_design_matrix_r
    end
    
    
    @testset "Compare scaled obervations with R" begin
        _, scaled_observations_julia = compute_tikhonov_matrix(spline_data, 0.5)
    
        scaled_observations_r = spline_model_r[:auxM][:XWy]
    
        @test scaled_observations_julia ≈ scaled_observations_r
    end
    
    
    @testset "R's storage of Tikhonov matrix" begin
        weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
        gram_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])
        tikhonov_matrix_r1 = weighted_design_matrix_r + spline_model_r[:lambda] * gram_r

        cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
        tikhonov_matrix_r2 = transpose(cholesky_r) * cholesky_r
        
        @test tikhonov_matrix_r1 ≈ tikhonov_matrix_r2
    end


    @testset "Check computation of Gram matrix entries" begin
        gram_julia = SmoothSpline.compute_gram_matrix(spline_data)

        gram_quadgk = similar(gram_julia)

        knots = unique(spline_data.B.Knots)
        max_integration_error = 0

        for j in 1:size(gram_julia, 2), i in 1:size(gram_julia, 1)
            gram_quadgk[i, j], integration_error = QuadGK.quadgk(
                u -> SmoothSpline.single_basis_function_deriv(i - 1, u, 2, spline_data.B)[2] * 
                SmoothSpline.single_basis_function_deriv(j - 1, u, 2, spline_data.B)[2],
                knots...
            )

            max_integration_error = max(integration_error, max_integration_error)
        end

        @test gram_quadgk ≈ transpose(gram_quadgk)
        @test gram_julia ≈ transpose(gram_julia)

        @test gram_julia ≈ gram_quadgk
        @test maximum(abs, gram_julia - gram_quadgk) < 10 * max_integration_error
    end

    
    @testset "Compare Gram matrix with R" begin
        gram_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])

        Mocking.activate()
        patch = @patch SmoothSpline.OneThird() = 0.333

        apply(patch) do
            gram_julia = SmoothSpline.compute_gram_matrix(spline_data)

            @test gram_r ≈ gram_julia
        end
    end
end
