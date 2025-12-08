using Adapt


@testset "Custom Parameterization" begin

    @kwdef struct SimpleAlbedo{NF<:Number} <: SpeedyWeather.AbstractAlbedo
        land_albedo::NF = 0.35
        seaice_albedo::NF = 0.6
        ocean_albedo::NF = 0.06
    end 

    Adapt.@adapt_structure SimpleAlbedo
    SimpleAlbedo(SG::SpectralGrid; kwargs...) = SimpleAlbedo{SG.NF}(; kwargs...)
    SpeedyWeather.initialize!(::SimpleAlbedo, model::PrimitiveEquation) = nothing 
    SpeedyWeather.variables(::SimpleAlbedo) = (
        DiagnosticVariable(name=:my_albedo, dims=Grid2D(), desc="Albedo", units="1"),
    )

    # note that for albedos the diagn::DiagnosticVariables is needed to indicate that this
    # is evaluated on the land-fraction averaged albedo in diagn. physics.my_albedo not independently for land and ocean
    # like the OceanLandAlbedo does
    # note that this conflicts with the shortwave radiation scheme though as this requires either a composite albedo like OceanLandAlbedo
    # or the `::DiagnosticVariables` removed where then `diagn` = `diagn.physics.land` and `diagn.physics.ocean` one after another
    function SpeedyWeather.parameterization!(ij, diagn::DiagnosticVariables, progn, albedo::SimpleAlbedo, model_parameters) 
        
        (; land_sea_mask) = model_parameters
        (; sea_ice_concentration) = progn.ocean
        (; land_albedo, seaice_albedo, ocean_albedo) = albedo
        
        if land_sea_mask.mask[ij] > 0.95 # if mostly land 
            diagn.physics.my_albedo[ij] = land_albedo
        else # if ocean
            diagn.physics.my_albedo[ij] = ocean_albedo + sea_ice_concentration[ij] * (seaice_albedo - ocean_albedo)
        end
    end 

    # replace the existing albedo with our custom albedo
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    albedo = SimpleAlbedo(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; albedo=albedo)
    simulation = initialize!(model) 
    
    # check if albedo in diagn 
    @test haskey(simulation.diagnostic_variables.physics, :my_albedo)
    
    run!(simulation, steps=3)       # run only a few time steps
    progn, diagn, model = SpeedyWeather.unpack(simulation)

    (; my_albedo) = diagn.physics
    (; mask) = model.land_sea_mask

    for ij in eachindex(my_albedo, mask)
        if mask[ij] > 0.95
            @test my_albedo[ij] == albedo.land_albedo
        else
            @test albedo.ocean_albedo <= my_albedo[ij] <= albedo.seaice_albedo
        end
    end

    # now do the same but with manually adding the albedo to the list of parameterizations
    model = PrimitiveWetModel(spectral_grid; custom_parameterization = SimpleAlbedo(spectral_grid), parameterizations=(:convection, :large_scale_condensation, :custom_parameterization, 
                                                            :surface_condition, :surface_momentum_flux, 
                                                            :surface_heat_flux, :surface_humidity_flux, 
                                                            :stochastic_physics))

    simulation = initialize!(model) 
    run!(simulation, steps=3)    # run only a few time steps
    progn, diagn, model = SpeedyWeather.unpack(simulation)

    (; my_albedo) = diagn.physics
    (; mask) = model.land_sea_mask

    for ij in eachindex(my_albedo, mask)
        if mask[ij] > 0.95
            @test my_albedo[ij] == albedo.land_albedo
        else
            @test albedo.ocean_albedo <= my_albedo[ij] <= albedo.seaice_albedo
        end
    end

end
