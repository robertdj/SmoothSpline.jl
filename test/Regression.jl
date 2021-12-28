using SmoothSpline
using Mocking
using Test

using RCall

using Random

@testset "Regression" begin
    Random.seed!(1)

    Mocking.activate()
    patch = @patch SmoothSpline.OneThird() = 0.333

    # smooth.spline scales x values to the closed unit interval before calling the compiled code
    # Test both scaled and unscaled here
    N = rand(3:10)
    x_values = [10 * rand(5), [0.0; rand(5); 1.0]]
    sort!.(x_values)

    @testset "Compare regression coefficients with R" for x in x_values
        y = rand(length(x))
    
        spline_model_r = RCall.rcopy(R"smooth.spline(x = $x, y = $y, spar = 0.5, keep.stuff = TRUE)")
        coefficients_r = spline_model_r[:fit][:coef]

        spline_data = SplineRegData(x, y)
    
        apply(patch) do
            coefficients_julia = SmoothSpline.regression(spline_data, 0.5)

            @test coefficients_julia ≈ coefficients_r
        end
    end
end
