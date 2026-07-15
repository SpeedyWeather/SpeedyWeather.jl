using SpeedyTransforms
using RingGrids, LowerTriangularArrays
using Test

include("spectral_transform.jl")
include("dispatch.jl")
include("spectral_gradients.jl")
include("power_spectrum.jl")
include("resolutions.jl")
include("array_utils.jl")
include("matrix_transform.jl")

# this will load Enzyme, EnzymeTestUtils
include("spectral_transform_ad_rules.jl")
