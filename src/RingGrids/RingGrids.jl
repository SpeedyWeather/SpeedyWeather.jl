#TODO turn into module?

include("grids_general.jl") # def AbstractGrid
include("show.jl")

include("full_grids.jl")    # def AbstractFullGrid, Full*Grid, * = Gaussian, Clenshaw, HEALPix, OctaHEALPix
include("octahedral.jl")    # def AbstractOctahedralGrid, OctahedralGaussianGrid, OctahedralClenshawGrid
include("healpix.jl")       # def AbstractHEALPixGrid, HEALPixGrid
include("octahealpix.jl")      # def AbstractOctaHEALPixGrid, OctaHEALPixGrid

include("quadrature_weights.jl")    # quadrature weights and solid angles for grids

include("interpolation.jl")