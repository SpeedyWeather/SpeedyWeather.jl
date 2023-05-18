Base.@kwdef struct SpectralGrid
    NF::Type{<:AbstractFloat} = Float32
    
    # HORIZONTAL
    trunc::Int = 31             # max degree of spherical harmonics
    Grid::Type{<:AbstractGrid} = OctahedralGaussianGrid # horizontal grid
    dealiasing::Float64 = 2     # dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid
    radius::Float64 = 6.371e6   # radius of the sphere [m]

    # SIZE OF GRID from trunc, Grid, dealiasing:
    # nlat_half is the number of latitude rings on one hemisphere (Equator incl)
    nlat_half::Int = SpeedyTransforms.get_nlat_half(trunc,dealiasing)
    npoints::Int = RingGrids.get_npoints(Grid,nlat_half)           # total number of grid points

    # VERTICAL
    nlev::Int = 8
    vertical_coordinates::VerticalCoordinates = SigmaCoordinates(;nlev)

    SpectralGrid(NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,nlev,vertical_coordinates) = nlev == vertical_coordinates.nlev ?
        new(NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,nlev,vertical_coordinates) :
        error("nlev does not match. $nlev vs $(vertical_coordinates.nlev)")
end

SpectralGrid(NF::Type{<:AbstractFloat};kwargs...) = SpectralGrid(NF=NF;kwargs...)
SpectralGrid(Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(Grid=Grid;kwargs...)
SpectralGrid(NF::Type{<:AbstractFloat},Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(Grid=Grid;kwargs...)

function Base.show(io::IO,SG::SpectralGrid)
    (;NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,nlev,vertical_coordinates) = SG
    truncation = if dealiasing < 2 "linear" elseif dealiasing < 3 "quadratic" else "cubic" end
    res = sqrt(4π*radius^2/npoints)/1000  # in [km]
    println(io,"Spectral:   T$trunc LowerTriangularMatrix{Complex{$NF}}, radius = $radius m")
    println(io,"Grid:       $npoints-element, $(get_nlat(Grid,nlat_half))-ring $Grid{$NF} ($truncation)")
    println(io,"Resolution: $(@sprintf("%.3g",res))km (average)")
      print(io,"Vertical:   $nlev-level $(typeof(vertical_coordinates))")
end

# TODO adjust nlev when vertical_coordinates is provided
# SpectralGrid(vertical_coordinates::VerticalCoordinates;kwargs...) = SpectralGrid(nlev=vertical_coordinates.nlev;vertical_coordinates,kwargs...)

"""
    G = Geometry(::SpectralGrid)

Construct Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct Geometry{NF<:AbstractFloat} <: AbstractGeometry{NF}     # NF: Number format

    "SpectralGrid that defines spectral and grid resolution"
    spectral_grid::SpectralGrid

    "grid of the dynamical core"
    Grid::Type{<:AbstractGrid} = spectral_grid.Grid

    "resolution parameter nlat_half of Grid, # of latitudes on one hemisphere (incl Equator)"
    nlat_half::Int = spectral_grid.nlat_half      


    # GRID-POINT SPACE
    "maximum number of longitudes (at/around Equator)"
    nlon_max::Int = get_nlon_max(Grid,nlat_half)

    "=nlon_max, same (used for compatibility), TODO: still needed?"
    nlon::Int = nlon_max

    "number of latitude rings"
    nlat::Int = get_nlat(Grid,nlat_half)

    "number of vertical levels"
    nlev::Int = spectral_grid.nlev

    "total number of grid points"
    npoints::Int = spectral_grid.npoints

    "Planet's radius [m]"
    radius::Float64 = spectral_grid.radius    


    # ARRAYS OF LANGITUDES/LONGITUDES
    "array of latitudes in degrees (90˚...-90˚)"
    latd::Vector{Float64} = get_latd(Grid,nlat_half)

    "array of longitudes in degrees (0...360˚), empty for non-full grids"
    lond::Vector{Float64} = get_lond(Grid,nlat_half)

    "longitude (-180˚...180˚) for each grid point in ring order"
    londs::Vector{NF} = get_latdlonds(Grid,nlat_half)[2]
    
    "latitude (-90˚...˚90) for each grid point in ring order"
    latds::Vector{NF} = get_latdlonds(Grid,nlat_half)[1]

    "sin of latitudes"
    sinlat::Vector{NF} = sind.(latd)
    
    "cos of latitudes"
    coslat::Vector{NF} = cosd.(latd)
    
    "= 1/cos(lat)"
    coslat⁻¹::Vector{NF} = 1 ./ coslat

    "= cos²(lat)"
    coslat²::Vector{NF} = coslat.^2

    "# = 1/cos²(lat)"
    coslat⁻²::Vector{NF} = 1 ./ coslat²            


    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::Vector{NF} = spectral_grid.vertical_coordinates.σ_half

    "σ at full levels, σₖ"
    σ_levels_full::Vector{NF} = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])  
    
    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::Vector{NF} = σ_levels_half[2:end] - σ_levels_half[1:end-1]      

    "log of σ at full levels, include surface (σ=1) as last element"
    ln_σ_levels_full::Vector{NF} = log.(vcat(σ_levels_full,1))
end

function Geometry(spectral_grid::SpectralGrid)
    return Geometry{spectral_grid.NF}(;spectral_grid)
end

function Base.show(io::IO,G::Geometry)
    println(io,"$(typeof(G))(")
    print(io,"$(G.spectral_grid))")
    # Base.show(io,G.σ_levels_half)
end

"""
    S = SpectralTransform(::SpectralGrid)

Generator function for a SpectralTransform struct pulling in parameters from a SpectralGrid struct."""
function SpeedyTransforms.SpectralTransform(spectral_grid::SpectralGrid;
                                            recompute_legendre::Bool = false,
                                            kwargs...)
    (;NF, Grid, trunc, dealiasing) = spectral_grid
    return SpectralTransform(NF,Grid,trunc;recompute_legendre,dealiasing,kwargs...)
end

