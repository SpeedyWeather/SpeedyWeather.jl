
"""Struct holding the tendencies of the prognostic spectral variables plus some additional tendencies used in their calculation"""
struct Tendencies{NF<:AbstractFloat}
    vor_tend        ::Array{Complex{NF},3}      # Vorticity of horizontal wind field
    div_tend        ::Array{Complex{NF},3}      # Divergence of horizontal wind field
    temp_tend       ::Array{Complex{NF},3}      # Absolute temperature [K]
    pres_surf_tend  ::Array{Complex{NF},2}      # Log of surface pressure [log(Pa)]
    humid_tend      ::Array{Complex{NF},3}      # Specific humidity [g/kg]
    u_tend          ::Array{Complex{NF},3}      # zonal velocity
    v_tend          ::Array{Complex{NF},3}      # meridonal velocity
end

"""
Defines the Tendencies struct.
"""
function Tendencies{NF}(G::GeoSpectral{NF}) where NF


    @unpack nlev = G.geometry
    @unpack mx, nx = G.spectral

    # conversion to type NF later when creating a PrognosticVariables struct
    vor_tend         = zeros(Complex{Float64},mx,nx,nlev)  # vorticity
    div_tend         = zeros(Complex{Float64},mx,nx,nlev)  # divergence
    temp_tend        = zeros(Complex{Float64},mx,nx,nlev)  # absolute Temperature
    pres_surf_tend   = zeros(Complex{Float64},mx,nx)       # logarithm of surface pressure
    humid_tend       = zeros(Complex{Float64},mx,nx,nlev)  # specific humidity
    u_tend           = zeros(Complex{Float64},mx,nx,nlev)  # zonal velocity
    v_tend           = zeros(Complex{Float64},mx,nx,nlev)  # meridonal velocity

    return Tendencies{NF}(vor_tend,div_tend,temp_tend,pres_surf_tend,humid_tend,u_tend,v_tend)

end


"""Struct holding the core prognostic spectral variables in grid point space, plus some additional grid quantities"""
struct GridVariables{NF<:AbstractFloat}
    vor_grid           ::Array{NF,3}  # Gridpoint field of vorticity
    div_grid           ::Array{NF,3}  # Gridpoint field of divergence
    temp_grid          ::Array{NF,3}  # Gridpoint field of absolute temperature [K]
    pres_surf_grid     ::Array{NF,2}  # Gridpoint field of surface pressure logarithm [log(Pa)]
    humid_grid         ::Array{NF,3}  # Gridpoint field of specific_humidity
    geopot_grid        ::Array{NF,3}  # Gridpoint field of geopotential
    tr_grid            ::Array{NF,3}  # Gridpoint field of tracers
    u_grid             ::Array{NF,3}  # Gridpoint field of zonal velocity
    v_grid             ::Array{NF,3}  # Gridpoint field of meridional velocity
    temp_grid_anomaly  ::Array{NF,3}  # Gridpoint field of absolute temperature anomaly [K]
end


"""
Defines the GridVariables struct.
"""
function GridVariables{NF}(G::GeoSpectral{NF}) where NF


    @unpack nlev,nlon,nlat = G.geometry

    # conversion to type NF later when creating a PrognosticVariables struct
    vor_grid           = zeros(Complex{Float64},nlon,nlat,nlev)  # vorticity
    div_grid           = zeros(Complex{Float64},nlon,nlat,nlev)  # divergence
    temp_grid          = zeros(Complex{Float64},nlon,nlat,nlev)  # absolute Temperature
    pres_surf_grid     = zeros(Complex{Float64},nlon,nlat)       # logarithm of surface pressure
    humid_grid        = zeros(Complex{Float64},nlon,nlat,nlev)  # specific humidity
    geopot_grid       = zeros(Complex{Float64},nlon,nlat,nlev)  # geopotential
    tr_grid           = zeros(Complex{Float64},nlon,nlat,nlev)  # tracers
    u_grid             = zeros(Complex{Float64},nlon,nlat,nlev)  # zonal velocity
    v_grid             = zeros(Complex{Float64},nlon,nlat,nlev)  # meridonal velocity
    temp_grid_anomaly  = zeros(Complex{Float64},nlon,nlat,nlev)  # absolute temperature anolamy

    return GridVariables{NF}(vor_grid,div_grid,temp_grid,pres_surf_grid,humid_grid,geopot_grid,tr_grid,u_grid,v_grid,temp_grid_anomaly)

end



"""Struct holding the diagnostic variables."""
struct DiagnosticVariables{NF<:AbstractFloat}
    tendencies     ::Tendencies{NF}
    grid_variables ::GridVariables{NF}
    #miscvars    ::MiscellaneousVariables{NF}
end

"""Generator function for Diagnostic Variables """
function DiagnosticVariables{NF}(G::GeoSpectral) where {NF<:AbstractFloat}
    tendencies_struct = Tendencies{NF}(G)
    grid_variables_struct = GridVariables{NF}(G)
    return DiagnosticVariables{NF}(tendencies_struct,grid_variables_struct)
end







# struct MiscellaneousVariables{NF<:AbstractFloat}
#     u_mean ::Array{NF,2}  # Mean gridpoint zonal velocity over all levels
#     v_mean ::Array{NF,2}  # Mean gridpoint meridional velocity over all levels
#     d_mean ::Array{NF,2}  # Mean gridpoint divergence over all levels

#     dumc ::Array{NF,3}  # Array for holding x/y gradients of surface pressure

#     px ::Array{NF,2}  # X Grad of pressure in grid point space 
#     py ::Array{NF,2}  # Y Grad of pressure in grid point space 


#     sigma_tend ::Array{NF,3} #vertical velocity in sigma coords
#     sigma_m    ::Array{NF,3}
#     puv :: Array{NF,3}

#     #Arbitrary array just to hold some temporary values. Used in zonal/meridonal_wind_tendency!(). 
#     #Named `temp` in Paxton/Chantry, but not a temperature quantity. `temp` for temporary?
#     #Better name needed
#     arbitrary_array ::Array{NF,3} 


#     d_meanc ::Array{NF,2} #used in get_spectral_tendencies  
#     sigma_tend_c ::Array{NF,3} 
#     dumk ::Array{NF,3} 
     
#  end

