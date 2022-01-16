module SmoothSpline

using Mocking

import LinearAlgebra
import OffsetArrays
import StatsBase

export
    smooth_spline

include("Splines.jl")
include("DesignMatrix.jl")
include("Regression.jl")

end # module
