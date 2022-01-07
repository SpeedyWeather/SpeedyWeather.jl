
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



"""Struct holding the core prognostic spectral variables in grid point space, plus some additional grid quantities"""
struct GridVariables{NF<:AbstractFloat}
    vor_grid       ::Array{NF,3}  # Gridpoint field of vorticity
    div_grid       ::Array{NF,3}  # Gridpoint field of divergence
    temp_grid      ::Array{NF,3}  # Gridpoint field of absolute temperature [K]
    pres_surf_grid ::Array{NF,2}  # Gridpoint field of surface pressure logarithm [log(Pa)]
    humid_grid     ::Array{NF,3}  # Gridpoint field of specific_humidity
    
    geopot_grid ::Array{NF,3}  # Gridpoint field of geopotential
    tr_grid     ::Array{NF,3}  # Gridpoint field of tracers
    u_grid      ::Array{NF,3}  # Gridpoint field of zonal velocity
    v_grid      ::Array{NF,3}  # Gridpoint field of meridional velocity
    temp_grid_anomaly   ::Array{NF,3}  # Gridpoint field of absolute temperature anomaly [K]

    
end


struct MiscellaneousVariables{NF<:AbstractFloat}
    u_mean ::Array{NF,2}  # Mean gridpoint zonal velocity over all levels
    v_mean ::Array{NF,2}  # Mean gridpoint meridional velocity over all levels
    d_mean ::Array{NF,2}  # Mean gridpoint divergence over all levels

    dumc ::Array{NF,3}  # Array for holding x/y gradients of surface pressure

    px ::Array{NF,2}  # X Grad of pressure in grid point space 
    py ::Array{NF,2}  # Y Grad of pressure in grid point space 


    sigma_tend ::Array{NF,3} #vertical velocity in sigma coords
    sigma_m    ::Array{NF,3}
    puv :: Array{NF,3}

    #Arbitrary array just to hold some temporary values. Used in zonal/meridonal_wind_tendency!(). 
    #Named `temp` in Paxton/Chantry, but not a temperature quantity. `temp` for temporary?
    #Better name needed
    arbitrary_array ::Array{NF,3} 


    d_meanc ::Array{NF,2} #used in get_spectral_tendencies  
    sigma_tend_c ::Array{NF,3} 
    dumk ::Array{NF,3} 
     
 end

"""Struct holding the diagnostic variables."""
struct DiagnosticVariables{NF<:AbstractFloat}
    tendencies  ::Tendencies{NF}
    gridvars    ::GridVariables{NF}
    miscvars    ::MiscellaneousVariables{NF}
end