"""
Compute the saturation specific humidity and relative hum. from specific hum. 
"""
function relative_and_saturation_humidity!(Diag::PrognosticVariables{NF}
                                           G::GeoSpectral{NF},             # Geometry and spectral struct
                                           C::Constants{NF} 
                                          ) where {NF<:AbstractFloat}



@unpack temp_grid,humid_grid_saturation,humid_grid_relative = Diag.gridvars
@unpack σ_levels_full = G.geometry
@unpack e0,c1,c2,t0,t1,t2,c3,c4 = C.humidity_constants
nlat,nlon,nlev = size(humid_grid_saturation)


#Iterate over every point. Can we vectorize this?
for k in 1:nlat
    for l in 1:nlon
        for m in 1:nlev


            c,t = temp_grid[k,l,m] >= zero(NF) ? c1,t1 : c2,t2
            humid_grid_saturation[k,l,m] = e0*exp(c * (temp_grid[k,l,m] - t0)/(temp_grid[k,l,m] - t))

            # if temp_grid[k,l,m] >=0
            #    #Saturation relative to liquid water
            #    humid_grid_saturation[k,l,m] = e0*exp(c1 * (temp_grid[k,l,m] - t0)/(temp_grid[k,l,m] - t1))
            # else
            #    #Saturation relative to ice
            #    humid_grid_saturation[k,l,m] = e0*exp(c2 * (temp_grid[k,l,m] - t0)/(temp_grid[k,l,m] - t2))
            # end if 


        end
    end
end


for m in 1:nlev
    if σ_levels_full[m] <= 0
        humid_grid_saturation[:,:,m]  = c3*humid_grid_saturation[:,:,m] / (pres_surf_grid_exp[1,1,m] - c4*humid_grid_saturation[:,:,m])

    else
        humid_grid_saturation  = c3*humid_grid_saturation / (σ_levels_full[m]*pres_surf_grid_exp - c4*humid_grid_saturation)

    end if

end

end