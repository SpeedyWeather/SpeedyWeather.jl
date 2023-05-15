Base.@kwdef struct SpectralGrid
    NF::Type{<:AbstractFloat} = Float32
    
    # HORIZONTAL
    trunc::Int = 31
    Grid::Type{<:AbstractGrid} = OctahedralGaussianGrid
    dealiasing::Float64 = 2

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