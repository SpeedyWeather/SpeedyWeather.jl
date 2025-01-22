# this is mainly there for testing purposes, and a bit more manual then above 
function flatten(prog::PrognosticVariables{NF, ArrayType, NSTEPS, SpectralVariable2D, SpectralVariable3D, GridVariable2D}) where {NF,ArrayType, NSTEPS, SpectralVariable2D, SpectralVariable3D, GridVariable2D}

    (; trunc, nlayers, nlat_half) = prog
    nvars = 5 + length(progn.tracers)

    # do it as a LTA 
    prog_array = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers, NSTEPS, nvars)

    for istep in 1:NSTEPS
        prog_array[:, :, istep, 1] = prog.vor[istep]
        prog_array[:, :, istep, 2] = prog.div[istep]
        prog_array[:, :, istep, 3] = prog.temp[istep]
        prog_array[:, :, istep, 4] = prog.humid[istep]
        prog_array[:, 1, istep, 5] = prog.pres[istep]

        for (i_key, (key, values)) in enumerate(prog.tracers)
            prog_array[:,:, istep, 5+i_key] = values
        end 
    end 

    nvars_grid = 6 
    land_ocean_array = zeros(GridVariable2D, nlat_half, nvars_grid)
    land_ocean_array[:, 1] = prog.ocean.sea_surface_temperature
    land_ocean_array[:, 2] = prog.ocean.sea_ice_concentration
    land_ocean_array[:, 3] = prog.land.land_surface_temperature
    land_ocean_array[:, 4] = prog.land.snow_depth
    land_ocean_array[:, 5] = prog.land.soil_moisture_layer1
    land_ocean_array[:, 6] = prog.land.soil_moisture_layer2

    return prog_array, land_ocean_array
end 
