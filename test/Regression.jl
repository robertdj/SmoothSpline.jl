using SmoothSpline
using Test

using RCall
import RDatasets
import LinearAlgebra

mtcars = RDatasets.dataset("datasets", "mtcars")
hp = mtcars[!, :HP]
mpg = mtcars[!, :MPG]

spline_model_r = RCall.rcopy(R"smooth.spline($hp, $mpg, spar = 0.5, keep.stuff = TRUE)")

# spline_data = SplineRegData(hp, mpg)
spline_data = SplineRegData(float.(hp), float.(mpg))

coefficients_julia = regression(spline_data, 0.5)

coefficients_r = spline_model_r[:fit][:coef]

@testset "Compare regression coefficients with R" begin
    @test coefficients_julia ≈ coefficients_r
    @test maximum(abs, coefficients_julia - coefficients_r) < 1e-12
end
