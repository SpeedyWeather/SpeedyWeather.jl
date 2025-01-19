module SpeedyWeatherFiniteDifferencesExt 

using SpeedyWeather 
import FiniteDifferences
import FiniteDifferences: to_vec 

# FiniteDifferences needs to be able to convert data structures to Vectors and back 
# This doesn't work out of the box with our data types, so we'll define those 
# conversions here.
function FiniteDifferences.to_vec(x::Grid) where Grid <: AbstractGridArray
    x_vec, from_vec = FiniteDifferences.to_vec(Array(x))

    function GridArray_from_vec(x_vec)
        return Grid(reshape(from_vec(x_vec), size(x)), x.nlat_half)
    end 

    return x_vec, GridArray_from_vec
end 

function FiniteDifferences.to_vec(x::LTA) where LTA <: LowerTriangularArray
    x_vec, from_vec = FiniteDifferences.to_vec(x.data)

    function LowerTriangularArray_from_vec(x_vec)
        return LowerTriangularArray(reshape(from_vec(x_vec), size(x)), x.m, x.n)
    end 

    return x_vec, LowerTriangularArray_from_vec
end

###
### PrognosticVariables
###

function flatten(prog::PrognosticVariables{NF, ArrayType, NSTEPS, SpectralVariable2D, SpectralVariable3D, GridVariable2D}) where {NF,ArrayType, NSTEPS, SpectralVariable2D, SpectralVariable3D, GridVariable2D}

    (; trunc, nlayers, nlat_half) = prog
    nvars = 5 + length(prog.tracers)

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

function SpeedyWeather.set!(prog::PrognosticVariables, L::LowerTriangularArray, G::AbstractGridArray)
    
    # the geometry isn't actually used for the set! we do 
    geometry = Geometry(SpectralGrid(trunc = prog.trunc, nlayers=prog.nlayers, Grid=typeof(prog.land.land_surface_temperature)))

    for istep in axes(L, 3)
        set!(prog, geometry; vor = L[:,:,istep,1], lf=istep)
        set!(prog, geometry; div = L[:,:,istep,2], lf=istep)
        set!(prog, geometry; temp = L[:,:,istep,3], lf=istep)
        set!(prog, geometry; humid = L[:,:,istep,4], lf=istep)
        set!(prog, geometry; pres = L[:,1,istep,5], lf=istep)

        for (i_key, (key, values)) in enumerate(prog.tracers)
            values .= L[:,:, istep, 5+i_key] 
        end 
    end 

    set!(prog, geometry; sea_surface_temperature = G[:,1])
    set!(prog, geometry; sea_ice_concentration = G[:,2])
    set!(prog, geometry; land_surface_temperature = G[:,3])
    set!(prog, geometry; snow_depth = G[:,4])
    set!(prog, geometry; soil_moisture_layer1 = G[:,5])
    set!(prog, geometry; soil_moisture_layer2 = G[:,6])

    return nothing
end 

function FiniteDifferences.to_vec(prog::PrognosticVariables)

    flattened_prog_lta, flattened_prog_grid = flatten(prog)

    flattened_prog_lta, lta_from_vec = FiniteDifferences.to_vec(flattened_prog_lta)
    flattened_prog_grid, grid_from_vec = FiniteDifferences.to_vec(flattened_prog_grid)

    N_lta = length(flattened_prog_lta)

    function PrognosticVariables_from_vec(x_vec)
        prog_new = zero(prog)
        set!(prog_new, lta_from_vec(x_vec[1:N_lta]), grid_from_vec(x_vec[N_lta+1:end])) 
        return prog_new 
    end 

    return cat(flattened_prog_lta, flattened_prog_grid, dims=1), PrognosticVariables_from_vec
end 

end 