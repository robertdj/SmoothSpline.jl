using SmoothSpline
using Mocking
using Test

using RCall
import RDatasets
import LinearAlgebra

using Random

@testset "Regression" begin
    mtcars = RDatasets.dataset("datasets", "mtcars")
    hp = mtcars[!, :HP]
    mpg = mtcars[!, :MPG]

    hp = 1:10
    mpg = 1:10

    # It seems to make a difference that the x's are equidistant
    hp = collect((0:3) / 3)
    mpg = collect(0:3)

    # Random.seed!(1)
    hp = [0.0; rand(5); 1.0] |> sort
    mpg = rand(length(hp))
    
    spline_model_r = RCall.rcopy(R"smooth.spline(x = $hp, y = $mpg, spar = 0.5, keep.stuff = TRUE)")
    
    # spline_data = SplineRegData(hp, mpg)
    spline_data = SplineRegData(float.(hp), float.(mpg))
    
    @testset "Compare regression coefficients with R" begin
        # coefficients_julia2 = SmoothSpline.regression2(spline_data, 0.5)
        coefficients_julia = SmoothSpline.regression(spline_data, 0.5)
        coefficients_r = spline_model_r[:fit][:coef]

        Mocking.activate()
        patch = @patch SmoothSpline.OneThird() = 0.333

        apply(patch) do
            coefficients_julia2 = SmoothSpline.regression2(spline_data, 0.5)

            @test coefficients_julia2 ≈ coefficients_r
            # @show maximum(abs, coefficients_julia2 - coefficients_r)
            # @test maximum(abs, coefficients_julia2 - coefficients_r) < 1e-12
        end


        # @test coefficients_julia ≈ coefficients_r
        # @test maximum(abs, coefficients_julia - coefficients_r) < 1e-12
        # @test maximum(abs, coefficients_julia2 - coefficients_r) < 1e-12
    end
end
