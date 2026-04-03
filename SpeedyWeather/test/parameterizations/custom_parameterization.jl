@testset "Custom Parameterization" begin

    @kwdef struct SimpleAlbedo{NF <: Number} <: SpeedyWeather.AbstractAlbedo
        land_albedo::NF = 0.35
        seaice_albedo::NF = 0.6
        ocean_albedo::NF = 0.06
    end

    SimpleAlbedo(SG::SpectralGrid; kwargs...) = SimpleAlbedo{SG.NF}(; kwargs...)
    SpeedyWeather.initialize!(::SimpleAlbedo, model::PrimitiveEquation) = nothing
    SpeedyWeather.variables(::SimpleAlbedo) = (
        SpeedyWeather.ParameterizationVariable(:my_albedo, SpeedyWeather.Grid2D(), desc = "Albedo", units = "1"),
    )

    # note that for albedos should actually define `albedo!(ij, vars, albedo, model)`
    # as they are applied to ocean/land separately but ignore this here and implement this like any other parameterization
    # in practice this conflicts with the shortwave radiation scheme that also doesn't know about the `my_albedo` variable
    function SpeedyWeather.parameterization!(ij, vars::Variables, albedo::SimpleAlbedo, model)

        (; land_sea_mask) = model
        (; sea_ice_concentration) = vars.prognostic.ocean
        (; land_albedo, seaice_albedo, ocean_albedo) = albedo

        if land_sea_mask.mask[ij] > 0.95 # if mostly land
            vars.parameterizations.my_albedo[ij] = land_albedo
        else # if ocean
            vars.parameterizations.my_albedo[ij] = ocean_albedo + sea_ice_concentration[ij] * (seaice_albedo - ocean_albedo)
        end
    end

    # replace the existing albedo with our custom albedo
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    albedo = SimpleAlbedo(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; albedo = albedo)
    simulation = initialize!(model)

    # check if albedo in vars.parameterizations
    @test haskey(simulation.variables.parameterizations, :my_albedo)

    run!(simulation, steps = 3)       # run only a few time steps
    vars, model = SpeedyWeather.unpack(simulation)

    (; my_albedo) = vars.parameterizations
    (; mask) = model.land_sea_mask

    for ij in eachindex(my_albedo, mask)
        if mask[ij] > 0.95
            @test my_albedo[ij] == albedo.land_albedo
        else
            @test albedo.ocean_albedo <= my_albedo[ij] <= albedo.seaice_albedo
        end
    end

    # now do the same but with manually adding the albedo to the list of parameterizations
    model = PrimitiveWetModel(
        spectral_grid;
        custom_parameterization = SimpleAlbedo(spectral_grid),
        parameterizations = (
            :convection, :large_scale_condensation, :custom_parameterization,
            :surface_condition, :surface_momentum_flux,
            :surface_heat_flux, :surface_humidity_flux,
            :stochastic_physics,
        )
    )

    simulation = initialize!(model)
    run!(simulation, steps = 3)    # run only a few time steps
    vars, model = SpeedyWeather.unpack(simulation)

    (; my_albedo) = vars.parameterizations
    (; mask) = model.land_sea_mask

    for ij in eachindex(my_albedo, mask)
        if mask[ij] > 0.95
            @test my_albedo[ij] == albedo.land_albedo
        else
            @test albedo.ocean_albedo <= my_albedo[ij] <= albedo.seaice_albedo
        end
    end
end
