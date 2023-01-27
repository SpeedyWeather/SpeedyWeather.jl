module RingGrids

    import Statistics: mean

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
            get_npoints,
            get_latdlonds,
            get_latd,
            each_index_in_ring,
            each_index_in_ring!,
            eachgridpoint,
            eachring,
            get_nlons,
            get_nlon_max,
            get_quadrature_weights,
            get_solid_angles

    # INTERPOLATION
    export  AbstractInterpolator,
            GridGeometry,
            AbstractLocator,
            AnvilLocator,
            AnvilInterpolator,
            DEFAULT_INTERPOLATOR

    export  interpolate,
            interpolate!,
            update_locator

    include("grids_general.jl")
    include("show.jl")
    include("full_grids.jl")
    include("octahedral.jl")
    include("healpix.jl")
    include("octahealpix.jl")
    include("quadrature_weights.jl")
end