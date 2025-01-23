# High-level tests whether time stepping in all models work

@testset "Differentiability: Timestepping" begin 
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    
    model_types = [ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel]
    model_type = PrimitiveDryModel
    for model_type in model_types 
        
        nlayer = model_type == ShallowWaterModel ? 1 : 3 

        spectral_grid = SpectralGrid(trunc=15, nlayers=nlayer)          # define resolution
        model = PrimitiveDryModel(; spectral_grid)   # construct model
        simulation = initialize!(model)  
        initialize!(simulation)

        (; prognostic_variables, diagnostic_variables, model) = simulation
        (; Δt, Δt_millisec) = model.time_stepping
        dt = 2Δt

        progn = prognostic_variables
        diagn = diagnostic_variables

        diagn_copy = deepcopy(diagn)
        progn_copy = deepcopy(progn)

        d_progn = zero(progn)
        d_diag = make_zero(diagn)
        d_model = make_zero(model)

        progn_new = zero(progn)
        dprogn_new = one(progn) # seed 

        # test that we can differentiate wrt an IC 
        autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, d_model))

        # nonzero gradient
        @test sum(to_vec(d_progn)[1]) != 0

        # FD comparison 
        dprogn_2 = one(progn) # seed 

        fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> timestep_oop(x, diagn_copy, dt, model), dprogn_2, progn_copy)

        @test isapprox(to_vec(fd_jvp[1])[1], to_vec(d_progn)[1])


        column = diagn.column 
        dcolumn = make_zero(column)
        heat_flux = model.surface_heat_flux 
        
        autodiff(Reverse, surface_heat_flux!, Const, Duplicated(column, dcolumn), Const(heat_flux), Const(model))

    end 

end


function surface_heat_flux!(   
    column::ColumnVariables{NF},
    heat_flux::SurfaceHeatFlux,
    model::PrimitiveEquation,
) where NF
    cₚ = model.atmosphere.heat_capacity
    (; heat_exchange_land, heat_exchange_sea) = heat_flux

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    T_skin_sea = column.skin_temperature_sea
    T_skin_land = column.skin_temperature_land
    T = column.surface_temp
    land_fraction = column.land_fraction

    # drag coefficient
    # the convert might seem redundant but without it Enzyme struggles with its type and activity 
    drag_sea, drag_land = heat_flux.use_boundary_layer_drag ?
                        (column.boundary_layer_drag, column.boundary_layer_drag) : 
                        (heat_exchange_sea, heat_exchange_land)

    # SPEEDY documentation Eq. 54 and 56, land/sea fraction included
    flux_land = ρ*drag_land*V₀*cₚ*(T_skin_land - T)*land_fraction
    flux_sea  = ρ*drag_sea*V₀*cₚ*(T_skin_sea  - T)*(one(NF)-land_fraction)

    # mix fluxes for fractional land-sea mask
    land_available = isfinite(T_skin_land)
    sea_available = isfinite(T_skin_sea)

    # Only flux from land/sea if available (not NaN) otherwise zero flux
    column.flux_temp_upward[end] += land_available ? flux_land : zero(NF)
    column.flux_temp_upward[end] += sea_available ? flux_sea : zero(NF)

    return nothing
end