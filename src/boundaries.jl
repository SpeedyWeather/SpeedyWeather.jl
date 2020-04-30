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

""" Generator function for a Boundaries struct."""
function Boundaries{T}( constants::Constants,
                        geometry::Geometry,
                        spectral_trans::SpectralTrans) where {T<:AbstractFloat}

    @unpack g = constants

    # LOAD NETCDF FILE
    path = "data"   #TODO pass on through parameters
    nc = NetCDF.open(joinpath(path,"surface.nc"))
    orog = nc.vars["orog"][:,end:-1:1]      # latitude is North to South in file
    lsm = nc.vars["lsm"][:,end:-1:1]
    alb = nc.vars["alb"][:,end:-1:1]

    # GEOPOTENTIAL, truncate in spectral space
    ϕ0 = g*orog
    ϕ0trunc = spectral_truncation(ϕ0,spectral_trans,geometry)

    # not necessary for orog, lsm, alb:
    # Fix undefined values
    # raw_input[raw_input .<= -999.0] .= 0.0
    # raw_input

    Boundaries{T}(ϕ0,ϕ0trunc,lsm,alb)
end
