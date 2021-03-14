module SmoothSpline

import LinearAlgebra
import OffsetArrays
import StatsBase

export
    find_span,
    basis_funs,
    basis_funs_deriv,
    BSpline

include("Splines.jl")

end # module
