import Random
using SmoothSpline
using Test

include("ExampleSpline.jl")

@testset "Splines" begin
    knots = [0.0, 1, 2, 3, 4, 4, 5]
    B = SmoothSpline.BSpline(knots, 2)

    # Example 2.2 on page 52 in "The NURBS Book"
    us = [0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0]
    indices = [2, 2, 3, 4, 5, 7, 7]
    @testset "Find interval for $u" for (u, index) in zip(us, indices)
        idx = SmoothSpline.find_span(u, B)

        @test idx == index
    end


    @testset "Error looking outside knots" for u in [-1, 6]
        @test_throws ArgumentError SmoothSpline.find_span(u, B)
    end


    # Property 2.4 on page 57 in "The NURBS Book"
    @testset "Splines are a division of unity" begin
        u = 5 * rand()
        idx = SmoothSpline.find_span(u, B)
        y = SmoothSpline.basis_funs(idx, u, B)

        @test sum(y) ≈ 1
    end


    # Example 2.2 on page 54 in "The NURBS Book"
    @testset "Specific spline values in $u" for u in 0:0.5:5
        rng = Random.Xoshiro(1)

        u += 0.05 * randn(rng)
        u = min(max(0, u), 5)

        knot_idx = SmoothSpline.find_span(u, B)
        y = SmoothSpline.basis_funs_deriv(knot_idx, u, 2, B)
        expected_y = example_spline.(u, knot_idx .- (B.P:-1:0))
        expected_y_diff1 = example_spline_diff1.(u, knot_idx .- (B.P:-1:0))
        expected_y_diff2 = example_spline_diff2.(u, knot_idx .- (B.P:-1:0))

        @test y[0, :].parent ≈ expected_y
        @test y[1, :].parent ≈ expected_y_diff1
        @test y[2, :].parent ≈ expected_y_diff2
    end
end
