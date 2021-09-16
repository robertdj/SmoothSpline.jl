using SmoothSpline
using Test

# import RCall
using RCall
import RDatasets
import QuadGK
import LinearAlgebra

import Random

include("ReshapeROutput.jl")
# include("ExampleSpline.jl")

@testset "Design matrix" begin
    mtcars = RDatasets.dataset("datasets", "mtcars")
    hp = mtcars[!, :HP]
    mpg = mtcars[!, :MPG]

    Random.seed!(1)
    hp = [0.0; 1/3; 0.5; 1.0]
    # hp = [0.0; rand(4); 1.0] |> sort
    mpg = rand(length(hp))
    # hp = 1:4
    # mpg = 1:4
    
    spline_model_r = RCall.rcopy(R"smooth.spline($hp, $mpg, spar = 0.5, keep.stuff = TRUE)")
    
    # spline_data = SplineRegData(hp, mpg)
    spline_data = SplineRegData(float.(hp), float.(mpg))
    
    
    # spline_reg_data = SplineRegData(hp, mpg)
    # coef = regression(spline_reg_data, spar = 0.5)
    
    @testset "Compare design matrix with R" begin
        design_matrix_julia = SmoothSpline.compute_design_matrix(spline_data)
        W = LinearAlgebra.diagm(spline_data.W)
        weighted_design_matrix = transpose(design_matrix_julia) * W * design_matrix_julia
    
        weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
    
        @test weighted_design_matrix ≈ weighted_design_matrix_r
        @test maximum(abs, weighted_design_matrix - weighted_design_matrix_r) < 1e-14
    end
    
    
    @testset "Compare scaled obervations with R" begin
        # design_matrix_julia = SmoothSpline.compute_design_matrix(spline_data)
        # W = LinearAlgebra.diagm(spline_data.W)
        # scaled_observations_julia = transpose(design_matrix_julia) * W * spline_data.Y
    
        _, scaled_observations_julia = compute_tikhonov_matrix(spline_data, 0.5)
    
        scaled_observations_r = spline_model_r[:auxM][:XWy]
    
        @test scaled_observations_julia ≈ scaled_observations_r
        @test maximum(abs, scaled_observations_julia - scaled_observations_r) < 1e-13
    end
    
    
    @testset "R's storage of Tikhonov matrix" begin
        weighted_design_matrix_r = vector_to_banded_matrix(spline_model_r[:auxM][:XWX])
        gram_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])
        tikhonov_matrix_r1 = weighted_design_matrix_r + spline_model_r[:lambda] * gram_r

        cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
        tikhonov_matrix_r2 = transpose(cholesky_r) * cholesky_r
        
        @test tikhonov_matrix_r1 ≈ tikhonov_matrix_r2
        @test maximum(abs, tikhonov_matrix_r1 - tikhonov_matrix_r2) < 1e-10
    end
    

    @testset "Check computation of Gram matrix entries" begin
        gram_julia = SmoothSpline.compute_gram_matrix(spline_data)

        gram_quadgk = similar(gram_julia)

        knots = unique(spline_data.B.Knots)
        max_integration_error = 0

        for i in 1:size(gram_julia, 1)
            for j in 1:size(gram_julia, 2)
                gram_quadgk[i, j], integration_error = QuadGK.quadgk(
                    u -> SmoothSpline.single_basis_function_deriv(i - 1, u, 2, spline_data.B)[2] * 
                    SmoothSpline.single_basis_function_deriv(j - 1, u, 2, spline_data.B)[2],
                    knots...
                )

                max_integration_error = max(integration_error, max_integration_error)
            end
        end

        @test gram_quadgk ≈ transpose(gram_quadgk)
        @test gram_julia ≈ transpose(gram_julia)

        @test max_integration_error <= 1e-8

        @test maximum(abs, gram_julia - gram_quadgk) < 10 * max_integration_error
    end

    
    @testset "Compare Gram matrix with R" begin
        gram_julia = SmoothSpline.compute_gram_matrix(spline_data)
    
        gram_r = vector_to_banded_matrix(spline_model_r[:auxM][:Sigma])

        @show maximum(abs, gram_julia - gram_r)
        # @test maximum(abs, gram_julia - gram_r) < 1e-5
    end
    
    
    @testset "Julia's Tikhonov matrix is better conditioned than R's" begin
        cholesky_r = vector_to_upper_triangular(spline_model_r[:auxM][:R])
        tikhonov_matrix_r = transpose(cholesky_r) * cholesky_r
    
        tikhonov_matrix_julia, _ = compute_tikhonov_matrix(spline_data, 0.5)

        condition_number_r = LinearAlgebra.cond(tikhonov_matrix_r)
        condition_number_julia = LinearAlgebra.cond(tikhonov_matrix_julia)

        @test condition_number_julia <= condition_number_r
    end
end
