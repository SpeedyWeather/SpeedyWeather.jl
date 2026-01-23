using SpeedyTransforms
using RingGrids, LowerTriangularArrays
using Test

include("spectral_transform.jl")
include("spectral_gradients.jl")
include("power_spectrum.jl")
include("resolutions.jl")

# this will load Enzyme, EnzymeTestUtils
if VERSION < v"1.12"
    include("spectral_transform_ad_rules.jl")
end
