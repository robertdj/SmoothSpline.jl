module SmoothSpline

using Mocking
using RCall
using RecipesBase

import LinearAlgebra
import OffsetArrays
import StatsBase

export
    BSpline,
    SplineRegData,

    basis_funs,
    basis_funs_deriv,
    find_span,
    regression,
    compute_tikhonov_matrix,
    predict,
    smooth_spline

include("Splines.jl")
include("DesignMatrix.jl")
include("Regression.jl")
include("Plot.jl")

end # module
