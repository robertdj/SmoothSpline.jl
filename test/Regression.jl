using SmoothSpline
using Mocking
using Test

using RCall
import RDatasets
import LinearAlgebra

using Random

@testset "Regression" begin
    mtcars = RDatasets.dataset("datasets", "mtcars")
    hp = mtcars[!, :HP] .|> float
    mpg = mtcars[!, :MPG] .|> float

    # Random.seed!(1)
    # hp = [0.0; rand(5); 1.0] |> sort
    # mpg = rand(length(hp))
    
    spline_model_r = RCall.rcopy(R"smooth.spline(x = $hp, y = $mpg, spar = 0.5, keep.stuff = TRUE)")
    
    spline_data = SplineRegData(hp, mpg)
    
    @testset "Compare regression coefficients with R" begin
        coefficients_r = spline_model_r[:fit][:coef]

        Mocking.activate()
        patch = @patch SmoothSpline.OneThird() = 0.333

        apply(patch) do
            coefficients_julia = SmoothSpline.regression(spline_data, 0.5)

            @test coefficients_julia ≈ coefficients_r
        end
    end
end
