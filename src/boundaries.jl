"""
Struct that holds the boundary arrays in grid-point space

    geopot_surf     ::Array{Complex{NF},2}  # spectral surface geopotential (orography * gravity) [m^2/s^2]
    landsea_mask    ::Array{NF,2}           # land-sea mask, grid-point
    albedo          ::Array{NF,2}           # annual mean surface albedo, grid-point
"""
struct Boundaries{NF<:AbstractFloat}        # number format NF
    geopot_surf     ::Array{Complex{NF},2}  # surface geopotential (i.e. orography) [m^2/s^2]
    landsea_mask    ::Array{NF,2}           # land-sea mask
    albedo          ::Array{NF,2}           # annual mean surface albedo
end

""" Generator function for a Boundaries struct. Loads the boundary conditions,
orography, land-sea mask and albedo from an netCDF file and stores the in a
`Boundaries` struct."""
function Boundaries(P::Parameters,
                    G::GeoSpectral)

    @unpack boundary_path, boundary_file, gravity = P

    # LOAD NETCDF FILE
    if boundary_path == ""
        path = joinpath(@__DIR__,"../input_data",boundary_file)
    else
        path = joinpath(boundary_path,boundary_file)
    end
    ncfile = NetCDF.open(path)
    
    # READ, latitude is North to South in file
    orography = ncfile.vars["orog"][:,end:-1:1]         # height [m]
    landsea_mask = ncfile.vars["lsm"][:,end:-1:1]       # fraction of land [0-1]
    albedo = ncfile.vars["alb"][:,end:-1:1]             # annual-mean albedo [0-1]

    # GEOPOTENTIAL
    # transform to spectral space and truncate
    P_inputdata = Parameters(trunc=31,nlat=48,nlon=96)  # orography only avail on default grid
    G_inputdata = GeoSpectral(P_inputdata)              # spectral transform for default grid    
    geopot_surf_lowres = spectral(gravity*orography,G_inputdata)
    spectral_truncation!(geopot_surf_lowres,G_inputdata)

    # interpolate orography by adding 0 coefficients to high resolution spectral geopotential
    geopot_surf = zeros(Complex{NF},G.spectral.mx,G.spectral.nx)    # highres with mx,nx from G
    @unpack mx,nx = G_inputdata.spectral                            # lowres mx,nx matching input data
    geopot_surf[1:mx,1:nx] .= geopot_surf_lowres

    Boundaries{P.NF}(geopot_surf,landsea_mask,albedo)
end
