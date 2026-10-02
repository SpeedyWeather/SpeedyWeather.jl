using SpeedyTransforms
using RingGrids
using LowerTriangularArrays
using JET
using Test

include("spectral_transform.jl")
include("type_name_length.jl")
include("dispatch.jl")
include("spectral_gradients.jl")
include("power_spectrum.jl")
include("resolutions.jl")
include("array_utils.jl")
include("matrix_transform.jl")

# this will load Enzyme, EnzymeTestUtils
include("spectral_transform_ad_rules.jl")
