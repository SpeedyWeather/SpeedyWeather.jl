"""
Boundaries holds the boundary arrays in grid-point space

ϕ0::Array{T,2}              # Surface geopotential (i.e. orography) [m^2/s^2]
ϕ0trunc::Array{T,2}         # Spectrally truncated surface geopotential [m^2/s^2]
land_sea_mask::Array{T,2}   # land-sea mask
albedo::Array{T,2}          # Annual mean surface albedo
"""
struct Boundaries{T<:AbstractFloat}
    ϕ0::Array{T,2}              # Surface geopotential (i.e. orography) [m^2/s^2]
    ϕ0trunc::Array{T,2}         # Spectrally truncated surface geopotential [m^2/s^2]
    land_sea_mask::Array{T,2}   # land-sea mask
    albedo::Array{T,2}          # Annual mean surface albedo
end

""" Generator function for a Boundaries struct. Loads the boundary conditions,
orography, land-sea mask and albedo from an netCDF file and stores the in a
Boundaries-struct."""
function Boundaries{T}( P::Params,
                        G::GeoSpectral{T}) where {T<:AbstractFloat}

    @unpack boundary_path, boundary_file, g = P

    # LOAD NETCDF FILE
    if boundary_path == ""  # default: take the surface.nc file in data
        path = joinpath(@__DIR__,"../data",boundary_file)
    else                    
        path = joinpath(boundary_path,boundary_file)
    end
    nc = NetCDF.open(path)
    orog = nc.vars["orog"][:,end:-1:1]  # latitude is North to South in file
    lsm = nc.vars["lsm"][:,end:-1:1]
    alb = nc.vars["alb"][:,end:-1:1]

    # GEOPOTENTIAL, truncate in spectral space
    ϕ0 = g*orog
    ϕ0trunc = spectral_truncation(ϕ0,G)

    Boundaries{T}(ϕ0,ϕ0trunc,lsm,alb)
end
