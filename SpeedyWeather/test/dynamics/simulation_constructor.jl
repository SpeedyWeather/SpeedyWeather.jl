@testset "Simulation(model) equals initialize!(model)" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)

    # Simulation(model) is an alias for initialize!(model)
    model = BarotropicModel(spectral_grid)
    simulation = Simulation(model)
    @test simulation isa SpeedyWeather.Simulation
    @test simulation.model === model

    # keyword arguments are forwarded, e.g. the initial time
    time = DateTime(2020, 5, 1)
    simulation = Simulation(BarotropicModel(spectral_grid), time = time)
    @test simulation.variables.prognostic.clock.time == time

    # both entry points give the same result
    model1 = ShallowWaterModel(spectral_grid)
    model2 = ShallowWaterModel(spectral_grid)
    simulation1 = initialize!(model1, time = time)
    simulation2 = Simulation(model2, time = time)
    @test simulation1.variables.prognostic.vorticity == simulation2.variables.prognostic.vorticity
    @test simulation1.variables.prognostic.clock.time == simulation2.variables.prognostic.clock.time
end

@testset "Printing Simulation types" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    simulation = Simulation(BarotropicModel(spectral_grid))

    # concrete simulation types are printed shortened
    @test sprint(show, typeof(simulation)) == "Simulation{Variables{…}, BarotropicModel{…}}"

    # partially specified types and types with unbound type variables (as in method signatures,
    # stack traces and docstrings) are printed as usual and don't throw. Print as seen from within
    # SpeedyWeather, otherwise names are module-qualified depending on what's loaded into Main
    show_unqualified(x) = sprint(show, x; context = :module => SpeedyWeather)
    @test show_unqualified(SpeedyWeather.Simulation) == "Simulation"
    @test show_unqualified(SpeedyWeather.Simulation{Int}) == "Simulation{Int64}"
    simulation_typevars = Base.unwrap_unionall(SpeedyWeather.Simulation)
    @test show_unqualified(simulation_typevars) == "Simulation{V, M<:AbstractModel}"
    @test show_unqualified(Tuple{Type{simulation_typevars}}) == "Tuple{Type{Simulation{V, M<:AbstractModel}}}"

    # the docstring of the Simulation struct (TYPEDSIGNATURES lists its constructors) can be parsed
    multidoc = Docs.meta(SpeedyWeather)[Docs.Binding(SpeedyWeather, :Simulation)]
    for docstr in values(multidoc.docs)
        @test Docs.parsedoc(docstr) isa Any
    end
end
