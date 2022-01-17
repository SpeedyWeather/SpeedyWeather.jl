struct parametrization_tendencies{NF<:AbstractFloat

#Just a placeholder

temp_convection ::Array{Complex{NF},3} #temperature tendency due to convection, grid point space
humid_convection ::Array{Complex{NF},3} #humidity tendency due to convection, grid point space

end

"""
Compute physical parametrization tendencies
"""
function parametrization_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                     Diag::PrognosticVariables{NF},# Diagnostic variables
                                     C::Constants{NF}
                                     ) where {T<:AbstractFloat}


    #Unpack relevant variables
    @unpack temp_tend, humid_tend,u_tend,v_tend = Diag.Tendencies
    @unpack pres_surf_grid,humid_grid,temp_grid,geopot_grid = Diag.gridvars

    nlat,nlon,nlev = size(humid_grid)


    # 1. Compute the thermodynamic variables

    # 1.1 Declare some surface pressure quantities
    pres_surf_grid_exp = exp(pres_surf_grid)
    pres_surf_grid_exp_inv = 1/pres_surf_grid_exp

    #1.2 Remove negative humidity grid values
    #TK: don't really want to iterate in this way. Is there a nice inbuilt function to do the same thing?
    for k in 1:nlat
        for l in 1:nlon
            for m in 1:nlev
                humid_grid[k,l,m] = max(humid_grid[k,l,m],0)
            end
        end
    end

    #1.3 Static energy as an anomaly from temperature and geopotential
    static_energy_grid = temp_grid + geopot_grid #Following Paxton/Chantry with a normalised cp


    #1.4...
    for k in 1:nlev
        compute_relative_and_saturation_humidity!(Diag,G) #defined in humidity.jl
    end


    # =========================================================================
    # 2.Precipitation (deep convection and large scale condensation)
    # =========================================================================


    # =========================================================================
    # 3.Radiation and surface fluxes 
    # =========================================================================


    # =========================================================================
    # 4.Planetary boundary layer interactions with lower troposphere 
    # =========================================================================


    #Sum of all tendencies
    u_tend = 0.0
    v_tend = 0.0
    temp_tend = temp_convection   # +... other terms
    humid_tend = humid_convection # +...

    #SSPT noise

   #Additive random noise



end




