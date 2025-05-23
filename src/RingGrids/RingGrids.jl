module RingGrids

# DOCUMENTATION
using  DocStringExtensions
import Printf

# NUMERICS
import Statistics: Statistics, mean
import FastGaussQuadrature
import LinearAlgebra
export rotate, rotate!

import Adapt
import GPUArrays
import KernelAbstractions

# ABSTRACT GRIDS
export  AbstractGrid,
        AbstractFullGrid,
        AbstractReducedGrid

# CONCRETE GRIDS
export  FullGaussianGrid,
        FullClenshawGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid

export  OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        HEALPixGrid,
        OctaHEALPixGrid,
        OctaminimalGaussianGrid

# FIELDS (Data on grids)
export  AbstractField, AbstractField2D, AbstractField3D,
        Field, Field2D, Field3D

export  FullGaussianField,
        FullClenshawField,
        FullHEALPixField,
        FullOctaHEALPixField,
        OctahedralGaussianField,
        OctahedralClenshawField,
        HEALPixField,
        OctaHEALPixField,
        OctaminimalGaussianField

# SIZE
export  grids_match,
        get_nlat,
        get_nlat_half,
        get_npoints,
        get_npoints2D

# COORDINATES
export  get_londlatds,
        get_lonlats,
        get_loncolats,
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
export  eachlayer,
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

# STATISTICS
export zonal_mean

include("utility_functions.jl")
include("architecture.jl")

# GENERAL
include("abstract_types.jl")
include("field.jl")
include("grid.jl")
include("scaling.jl")
include("geodesics.jl")
include("reverse.jl")
include("rotate.jl")

# FULL GRIDS
include("grids/full_grids.jl")
include("grids/full_gaussian.jl")
include("grids/full_clenshaw.jl")
include("grids/full_healpix.jl")
include("grids/full_octahealpix.jl")

# REDUCED GRIDS
include("grids/reduced_grids.jl")
include("grids/octahedral_gaussian.jl")
include("grids/octahedral_clenshaw.jl")
include("grids/healpix.jl")
include("grids/octahealpix.jl")
include("grids/octaminimal_gaussian.jl")

# INTEGRATION AND INTERPOLATION
include("quadrature_weights.jl")
include("interpolation.jl")
include("vertices.jl")
include("statistics.jl")

end