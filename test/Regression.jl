using SmoothSpline
using Mocking
using Test

using RCall

using Random

@testset "Regression" begin
    Random.seed!(1)

    Mocking.activate()
    patch = @patch SmoothSpline.OneThird() = 0.333

    # smooth.spline scales x values to the closed unit interval before calling the compiled code.
    # Test both scaled and unscaled here
    N = rand(5:20)
    x_values = [10 * rand(N), [0.0; rand(N); 1.0]]
    sort!.(x_values)

    spar_values = rand(2)

    @testset "Compare regression coefficients with R" for x in x_values, spar in spar_values
        y = rand(length(x))
    
        spline_model_r = RCall.rcopy(R"smooth.spline(x = $x, y = $y, spar = $spar, keep.stuff = TRUE)")
        coefficients_r = spline_model_r[:fit][:coef]

        spline_data = SplineRegData(x, y)
    
        apply(patch) do
            coefficients_julia = SmoothSpline.regression(spline_data, spar)

            @test coefficients_julia ≈ coefficients_r
        end
    end
end
