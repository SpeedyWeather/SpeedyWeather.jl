module RingGrids

# DOCUMENTATION AND VISUALISATION
using  DocStringExtensions
import UnicodePlots

# NUMERICS
import Statistics: mean
import FastGaussQuadrature
import LinearAlgebra

# GRIDS
export  AbstractGrid,
        AbstractFullGrid,
        AbstractOctahedralGrid,
        AbstractHEALPixGrid,
        AbstractOctaHEALPixGrid

export  FullGaussianGrid,
        FullClenshawGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid

export  OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        HEALPixGrid,
        OctaHEALPixGrid

# GRID FUNCTIONS
export  grids_match,
        get_truncation,
        get_resolution,
        get_nlat,
        get_nlat_half,
        get_npoints,
        get_latdlonds,
        get_lat,
        get_colat,
        get_latd,
        get_lond,
        each_index_in_ring,
        each_index_in_ring!,
        eachgridpoint,
        eachring,
        whichring,
        get_nlons,
        get_nlon_max,
        get_quadrature_weights,
        get_solid_angles

# SCALING
export  scale_coslat!,
        scale_coslat²!,
        scale_coslat⁻¹!,
        scale_coslat⁻²!

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

# export  plot

include("utility_functions.jl")

include("grids_general.jl")
include("full_grids.jl")
include("octahedral.jl")
include("healpix.jl")
include("octahealpix.jl")
include("quadrature_weights.jl")
include("interpolation.jl")
include("show.jl")
include("similar.jl")
include("scaling.jl")

end