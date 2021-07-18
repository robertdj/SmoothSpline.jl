module SmoothSpline

import LinearAlgebra
import OffsetArrays
import StatsBase

export
    BSpline,
    SplineRegData,

    basis_funs,
    basis_funs_deriv,
    find_span,
    regression

include("Splines.jl")
include("DesignMatrix.jl")
include("Regression.jl")

end # module
