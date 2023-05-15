Base.@kwdef struct SpectralGrid
    NF::Type{<:AbstractFloat} = Float32
    
    # HORIZONTAL
    trunc::Int = 31             # max degree of spherical harmonics
    Grid::Type{<:AbstractGrid} = OctahedralGaussianGrid # horizontal grid
    dealiasing::Float64 = 2     # dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid

    # SIZE OF GRID from trunc, Grid, dealiasing:
    # nlat_half is the number of latitude rings on one hemisphere (Equator incl)
    nlat_half = SpeedyTransforms.get_nlat_half(trunc,dealiasing)
    npoints = RingGrids.get_npoints(Grid,nlat_half)           # total number of grid points

    # VERTICAL
    nlev::Int = 8
    vertical_coordinates::VerticalCoordinates = SigmaCoordinates(;nlev)

    SpectralGrid(NF,trunc,Grid,dealiasing,nlev,vertical_coordinates) = nlev == vertical_coordinates.nlev ?
        new(NF,trunc,Grid,dealiasing,nlev,vertical_coordinates) :
        error("nlev does not match. $nlev vs $(vertical_coordinates.nlev)")
end

SpectralGrid(NF::Type{<:AbstractFloat};kwargs...) = SpectralGrid(NF=NF;kwargs...)
SpectralGrid(Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(Grid=Grid;kwargs...)
SpectralGrid(NF::Type{<:AbstractFloat},Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(Grid=Grid;kwargs...)

# TODO adjust nlev when vertical_coordinates is provided
# SpectralGrid(vertical_coordinates::VerticalCoordinates;kwargs...) = SpectralGrid(nlev=vertical_coordinates.nlev;vertical_coordinates,kwargs...)