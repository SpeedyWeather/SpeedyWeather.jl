"""
Struct that holds the boundary arrays in grid-point space

    ϕ0::Array{NF,2}              # surface geopotential [m^2/s^2]
    ϕ0trunc::Array{NF,2}         # spectrally truncated surface geopotential [m^2/s^2]
    land_sea_mask::Array{NF,2}   # land-sea mask
    albedo::Array{NF,2}          # annual mean surface albedo
"""
struct Boundaries{NF<:AbstractFloat}    # number format NF
    ϕ0::Array{NF,2}              # surface geopotential (i.e. orography) [m^2/s^2]
    ϕ0trunc::Array{NF,2}         # spectrally truncated surface geopotential [m^2/s^2]
    land_sea_mask::Array{NF,2}   # land-sea mask
    albedo::Array{NF,2}          # annual mean surface albedo
end

""" Generator function for a Boundaries struct. Loads the boundary conditions,
orography, land-sea mask and albedo from an netCDF file and stores the in a
Boundaries-struct."""
function Boundaries{NF}(P::Params,
                        G::GeoSpectral{NF}) where {NF<:AbstractFloat}

    @unpack boundary_path, boundary_file, g = P

    # LOAD NETCDF FILE
    if boundary_path == ""
        path = joinpath(@__DIR__,"../input_data",boundary_file)
    else
        path = joinpath(boundary_path,boundary_file)
    end
    ncfile = NetCDF.open(path)
    orog = ncfile.vars["orog"][:,end:-1:1]  # latitude is North to South in file
    lsm = ncfile.vars["lsm"][:,end:-1:1]
    alb = ncfile.vars["alb"][:,end:-1:1]

    # GEOPOTENTIAL, truncate in spectral space
    ϕ0 = g*orog
    ϕ0trunc = spectral_truncation(ϕ0,G)

    Boundaries{NF}(ϕ0,ϕ0trunc,lsm,alb)
end
