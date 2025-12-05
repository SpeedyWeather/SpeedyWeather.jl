using Adapt


@testset "Custom Parameterization" begin

    @kwdef struct SimpleAlbedo{NF<:Number} <: SpeedyWeather.AbstractParameterization
        land_albedo::NF=0.35
        seaice_albedo::NF=0.6
        ocean_albedo::NF=0.06
    end 

    Adapt.@adapt_structure SimpleAlbedo

    SimpleAlbedo(SG::SpectralGrid; kwargs...) = SimpleAlbedo(; kwargs...)

    SpeedyWeather.initialize!(::SimpleAlbedo, model::PrimitiveEquation) = nothing 

    SpeedyWeather.variables(::SimpleAlbedo) = (
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo", units="1"),
    )

    function SpeedyWeather.parameterization!(ij, diagn::DiagnosticVariables, progn::PrognosticVariables, albedo::SimpleAlbedo, model_parameters) 
        
        (; land_sea_mask) = model_parameters
        (; sea_ice_concentration) = progn.ocean
        (; land_albedo, seaice_albedo, ocean_albedo) = albedo
        
        if land_sea_mask.mask[ij] > 0.95 # if mostly land 
            diagn.physics.albedo[ij] = land_albedo
        else # if ocean
            diagn.physics.albedo[ij] = ocean_albedo + sea_ice_concentration[ij] * (seaice_albedo - ocean_albedo)
        end
    end 

    # replace the existing albedo with our custom albedo
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    albedo = SimpleAlbedo()
    model = PrimitiveWetModel(spectral_grid; albedo=albedo)
    simulation = initialize!(model) 
    run!(simulation, period=Day(5)) # spin up the model a little
    progn, diagn, model = SpeedyWeather.unpack(simulation)

    # check if albedo in diagn 
    @test !isnothing(diagn.physics.albedo)

    # check if albedo is set correctly 
    #heatmap(diagn.physics.albedo)

    @test all(diagn.physics.albedo[model.land_sea_mask.mask .> 0.95] .≈ albedo.land_albedo)


    # now do the same but with manually adding the albedo to the list of parameterizations
    model = PrimitiveWetModel(spectral_grid; custom_parameterization = SimpleAlbedo(), parameterizations=(:convection, :large_scale_condensation, :custom_parameterization, 
                                                            :surface_condition, :surface_momentum_flux, 
                                                            :surface_heat_flux, :surface_humidity_flux, 
                                                            :stochastic_physics))

    simulation = initialize!(model) 
    run!(simulation, period=Day(5)) # spin up the model a little
    progn, diagn, model = SpeedyWeather.unpack(simulation)

    # check if albedo is set correctly 
    #heatmap(diagn.physics.albedo)

    @test all(diagn.physics.albedo[model.land_sea_mask.mask .> 0.95] .≈ albedo.land_albedo)

end
