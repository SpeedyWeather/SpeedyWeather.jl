"""
Struct holding quantities calculated from the physical parameterisations. All quantities
are in grid-point space.
"""
struct ParametrizationVariables{NF<:AbstractFloat}
    sat_vap_pressure   ::Array{NF,3}   # Saturation vapour pressure
    sat_spec_humidity  ::Array{NF,3}   # Saturation specific humidity
    cloud_top          ::Array{Int,2}  # Cloud-top
    precip_large_scale ::Array{NF,2}   # Large-scale precipitation
    humid_tend_lsc     ::Array{NF,3}   # Humidity tendencies due to large-scale condensation
    temp_tend_lsc      ::Array{NF,3}   # Temperature tendencies due to large-scale condensation
end

"""
Generator function for the ParametrizationVariables struct. Initialises with zeros.
"""
function ParametrizationVariables(G::GeoSpectral{NF}) where NF
    @unpack nlon, nlat, nlev = G.geometry

    sat_vap_pressure   = zeros(NF,nlon,nlat,nlev)  # Saturation vapour pressure
    sat_spec_humidity  = zeros(NF,nlon,nlat,nlev)  # Saturation specific humidity
    cloud_top          = zeros(Int,nlon,nlat)      # Cloud-top
    precip_large_scale = zeros(NF,nlon,nlat)       # Large-scale precipitation
    humid_tend_lsc     = zeros(NF,nlon,nlat,nlev)  # Humidity tendencies due to large-scale condensation
    temp_tend_lsc      = zeros(NF,nlon,nlat,nlev)  # Temperature tendencies due to large-scale condensation

    return ParametrizationVariables(sat_vap_pressure,
                                    sat_spec_humidity,
                                    cloud_top,
                                    precip_large_scale,
                                    humid_tend_lsc,
                                    temp_tend_lsc,
                                    )
end