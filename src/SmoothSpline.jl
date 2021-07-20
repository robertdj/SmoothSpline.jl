module SmoothSpline

using RCall

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
    compute_tikhonov_matrix

include("Splines.jl")
include("DesignMatrix.jl")
include("Regression.jl")

end # module
