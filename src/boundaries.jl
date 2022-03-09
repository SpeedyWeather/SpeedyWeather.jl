"""
Struct that holds the boundary arrays in grid-point space
    geopot_surf     ::Array{Complex{NF},2}  # spectral surface geopotential (orography * gravity) [m^2/s^2]
    landsea_mask    ::Array{NF,2}           # land-sea mask, grid-point
    albedo          ::Array{NF,2}           # annual mean surface albedo, grid-point
"""
struct Boundaries{NF<:AbstractFloat}        # number format NF
    geopot_surf     ::Array{Complex{NF},2}  # spectral surface geopotential (i.e. orography) [m^2/s^2]
    # landsea_mask    ::Array{NF,2}           # land-sea mask
    # albedo          ::Array{NF,2}           # annual mean surface albedo
end

""" Generator function for a Boundaries struct. Loads the boundary conditions,
orography, land-sea mask and albedo from an netCDF file and stores the in a
`Boundaries` struct."""
function Boundaries(P::Parameters,
                    G::GeoSpectral)

    @unpack orography_path, orography_file, gravity, NF = P

    # LOAD NETCDF FILE
    if orography_path == ""
        path = joinpath(@__DIR__,"../input_data",orography_file)
    else
        path = joinpath(orography_path,orography_file)
    end
    ncfile = NetCDF.open(path)
    
    # READ, TODO check which latitude ordering is required, it's North to South in file
    orography = ncfile.vars["orog"][:,:]        # height [m]
    # landsea_mask = ncfile.vars["lsm"][:,:]    # fraction of land [0-1]
    # albedo = ncfile.vars["alb"][:,:]          # annual-mean albedo [0-1]

    # GEOPOTENTIAL, transform to spectral space and truncate accordingly
    nlon_orography,nlat_orography = size(orography)
    P_inputdata = Parameters(trunc=P.trunc,
                                nlon=nlon_orography,
                                nlat=nlat_orography)
    G_inputdata = GeoSpectral(P_inputdata)              # spectral transform for default grid    
    geopot_surf = spectral(gravity*orography,G_inputdata)
    spectral_truncation!(geopot_surf,G_inputdata)

    # # interpolate orography by adding 0 coefficients to high resolution spectral geopotential
    # geopot_surf = zeros(Complex{NF},G.spectral.mx,G.spectral.nx)    # highres with mx,nx from G
    # @unpack mx,nx = G_inputdata.spectral                            # lowres mx,nx matching input data
    # geopot_surf[1:mx,1:nx] .= geopot_surf_lowres

    Boundaries{NF}(geopot_surf) #,landsea_mask,albedo)
end