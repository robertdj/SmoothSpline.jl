using Mocking
using Random
using SmoothSpline
using Test

using RCall


@testset "Regression" begin
    rng = Random.Xoshiro(1)

    # smooth.spline scales x values to the closed unit interval before calling the compiled code.
    # Results should not depend on such an assumption
    N = rand(rng, 5:20)
    obs_x = 10 * rand(rng, N)
    obs_y = 10 * rand(rng, N)

    spar = rand(rng)
    spline_model_r = RCall.rcopy(R"spline_model <- smooth.spline(x = $obs_x, y = $obs_y, spar = $spar)")

    Mocking.activate()
    patch = @patch SmoothSpline.OneThird() = 0.333

    spline_model_julia = apply(patch) do
        SmoothSpline.smooth_spline(obs_x, obs_y, spar)
    end

    @testset "Compare regression coefficients with R" begin
        coefficients_r = spline_model_r[:fit][:coef]

        @test spline_model_julia.Coef ≈ coefficients_r
    end

    @testset "Compare predictions with R" begin
        x = range(extrema(obs_x)...; length = 100)

        predictions_r = RCall.rcopy(R"predict(spline_model, $x)")
        
        predictions_julia = SmoothSpline.predict(spline_model_julia, x)

        @test predictions_julia ≈ predictions_r[:y]
    end
end
