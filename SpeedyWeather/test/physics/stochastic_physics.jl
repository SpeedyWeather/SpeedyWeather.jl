@testset "Surface perturbation wet/dry model" begin
    @kwdef struct MySurfacePerturbation <: SpeedyWeather.AbstractSurfacePerturbation
        noise_scale_temp::Float32 = 1            # standard deviation of perturbation, in Kelvin
    end
    
    # implement the actual perturbation as a functor 
    function (SP::MySurfacePerturbation)(
        column::ColumnVariables,
        model::PrimitiveEquation,
    )
        (; temp, humid, nlayers) = column
        T, q = temp[nlayers], humid[nlayers]    # unperturbed lowermost layer variables
    
        # perturb temperature additive only
        r = column.random_value
        T += SP.noise_scale_temp*r
    
        return T, q
    end

    spectral_grid = SpectralGrid()
    my_surface_perturbation = MySurfacePerturbation(noise_scale_temp=1)

    # construct convection with that perturbation
    convection = SimplifiedBettsMiller(spectral_grid, surface_temp_humid=my_surface_perturbation)

    # don't forget the random process otherwise your `column.random_value` is always 0!
    random_process = SpectralAR1Process(spectral_grid, seed=1)

    # construct model
    model = PrimitiveWetModel(spectral_grid; random_process, convection)
    simulation = initialize!(model)

    run!(simulation, period=Day(5))

    k = spectral_grid.nlayers
    surface_humid = simulation.diagnostic_variables.grid.humid_grid[:, k]
    
    # reconstruct model with different seed
    random_process = SpectralAR1Process(spectral_grid, seed=2)
    model = PrimitiveWetModel(spectral_grid; random_process, convection)
    # reinitialize, rerun check whether it's different!
    simulation = initialize!(model)     # choose a different random seed for random_process
    run!(simulation, period=Day(5))
    surface_humid2 = simulation.diagnostic_variables.grid.humid_grid[:, k]

    for ij in eachindex(surface_humid, surface_humid2)
        if surface_humid[ij] != 0
            @test surface_humid[ij] != surface_humid2[ij]
        end
    end

    # for the dry model too
    convection = DryBettsMiller(spectral_grid, surface_temp=my_surface_perturbation)
    model = PrimitiveDryModel(spectral_grid; random_process, convection)
    simulation = initialize!(model)
    run!(simulation, period=Day(5))

    k = spectral_grid.nlayers
    surface_temp = simulation.diagnostic_variables.grid.temp_grid[:, k]
    
    # reinitialize, rerun check whether it's different!
    simulation = initialize!(model)     # choose a different random seed for random_process
    run!(simulation, period=Day(5))
    surface_temp2 = simulation.diagnostic_variables.grid.temp_grid[:, k]

    # count the number of grid points that are different in temperature
    N = length(surface_temp)
    local n::Int64 = 0
    for ij in eachindex(surface_temp, surface_temp2)
        n += surface_temp[ij] != surface_temp2[ij]
    end

    # 75% of grid points must be different
    @test n/N > 3/4
end