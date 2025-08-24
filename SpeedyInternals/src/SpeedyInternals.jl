"""
    SpeedyInternals

Common (internal) functions for SpeedyWeather.jl and related packages 
like RingGrids.jl, LowerTriangularArrays.jl and SpeedyTransforms.jl
"""
module SpeedyInternals

    # ARCHITECTURES (Device handling)
    include("Architectures/Architectures.jl")

    # UTILS (Kernel launching and various utilities)
    include("Utils/Utils.jl")

end
