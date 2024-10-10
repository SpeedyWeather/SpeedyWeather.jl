module RingGrids

# DOCUMENTATION AND VISUALISATION
using  DocStringExtensions
import UnicodePlots

# NUMERICS
import Statistics: mean
import FastGaussQuadrature
import LinearAlgebra

# GPU
import Adapt
import GPUArrays

# ABSTRACT GRIDS (2D) AND GRIDARRAYS (3D+)
export  AbstractGridArray,
        AbstractFullGridArray,
        AbstractReducedGridArray

export  AbstractGrid,
        AbstractFullGrid,
        AbstractReducedGrid

# CONCRETE GRIDS (2D) AND GRIDARRAYS (3D+)
export  FullGaussianArray,
        FullClenshawArray,
        FullHEALPixArray,
        FullOctaHEALPixArray

export  FullGaussianGrid,
        FullClenshawGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid

export  OctahedralGaussianArray,
        OctahedralClenshawArray,
        HEALPixArray,
        OctaHEALPixArray

export  OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        HEALPixGrid,
        OctaHEALPixGrid

# SIZE
export  grids_match,
        get_nlat,
        get_nlat_half,
        get_npoints2D

# COORDINATES
export  get_latdlonds,
        get_lat,
        get_colat,
        get_latd,
        get_lond,
        get_nlons,
        get_nlon_max

# INTEGRATION
export  get_quadrature_weights,
        get_solid_angles

# ITERATORS
export  eachgrid,
        eachring,
        whichring,
        eachgridpoint,
        each_index_in_ring,
        each_index_in_ring!

# SCALING
export  scale_coslat!,
        scale_coslat²!,
        scale_coslat⁻¹!,
        scale_coslat⁻²!,
        scale_coslat,
        scale_coslat²,
        scale_coslat⁻¹,
        scale_coslat⁻²

# INTERPOLATION
export  AbstractInterpolator,
        GridGeometry,
        AbstractLocator,
        AnvilLocator,
        AnvilInterpolator,
        DEFAULT_INTERPOLATOR

export  interpolate,
        interpolate!,
        update_locator,
        update_locator!

include("utility_functions.jl")

# GENERAL
include("general.jl")
include("full_grids.jl")
include("reduced_grids.jl")
include("scaling.jl")

# FULL GRIDS
include("grids/full_gaussian.jl")
include("grids/full_clenshaw.jl")
include("grids/full_healpix.jl")
include("grids/full_octahealpix.jl")

# REDUCED GRIDS
include("grids/octahedral_gaussian.jl")
include("grids/octahedral_clenshaw.jl")
include("grids/healpix.jl")
include("grids/octahealpix.jl")

# INTEGRATION AND INTERPOLATION
include("quadrature_weights.jl")
include("interpolation.jl")

# OUTPUT
include("show.jl")

end