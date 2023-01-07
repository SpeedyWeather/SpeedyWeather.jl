#TODO turn into module?

include("grids_general.jl") # def AbstractGrid

include("full_grids.jl")    # def AbstractFullGrid, Full*Grid, * = Gaussian, Clenshaw, HEALPix, HEALPix4
include("octahedral.jl")    # def AbstractOctahedralGrid, OctahedralGaussianGrid, OctahedralClenshawGrid
include("healpix.jl")       # def AbstractHEALPixGrid, HEALPixGrid
include("healpix4.jl")      # def AbstractHEALPix4Grid, HEALPix4Grid

include("quadrature_weights.jl")    # quadrature weights and solid angles for grids

include("interpolation.jl")
