
"""
Humidity struct containing relevant parameters and arrays
"""
struct compute_relative_and_saturation_humidity{NF<:AbstractFloat}     

#Set of constants used in the saturation humidity calculations
#Where do these come from? What do they represent?

e0 ::Float64= 6.108e-3
c1 ::Float64= 17.269
c2 ::Float64= 21.875
t0 ::Float64= 273.16 #Paxton/Chantry have a zero_c term here? 
t1 ::Float64= 35.86 
t2 ::Float64= 7.66 

end



"""
Compute the saturation specific humidity and relative hum. from specific hum. 
"""
function compute_relative_and_saturation_humidity!(Diag::PrognosticVariables{NF}
                                                   G::GeoSpectral{NF},             # Geometry and spectral struct 
                                                   ) where {NF<:AbstractFloat}



@unpack temp_grid,humid_grid_saturation,humid_grid_relative = Diag.gridvars
@unpack σ_levels_full = G.geometry

nlat,nlon,nlev = size(humid_grid_saturation)


σ_levels_full

#Iterate over every point. Can we vectorize this?
for k in 1:nlat
    for l in 1:nlon
        for m in 1:nlev

            if temp_grid[k,l,m] >=0
               #Saturation relative to liquid water
               humid_grid_saturation[k,l,m] = e0*exp(c1 * (temp_grid[k,l,m] - t0)/(temp_grid[k,l,m] - t1))
            else
               #Saturation relative to ice
               humid_grid_saturation[k,l,m] = e0*exp(c2 * (temp_grid[k,l,m] - t0)/(temp_grid[k,l,m] - t2))
            end if 


        end
    end
end


for m in 1:nlev
    if σ_levels_full[m] <= 0
        humid_grid_saturation[:,:,m]  = 622.0_*humid_grid_saturation[:,:,m] / (pres_surf_grid_exp[1,1,m] - 0.378*humid_grid_saturation[:,:,m])

    else
        humid_grid_saturation  = 622.0_*humid_grid_saturation / (σ_levels_full[m]*pres_surf_grid_exp - 0.378*humid_grid_saturation)

    end if

end

end