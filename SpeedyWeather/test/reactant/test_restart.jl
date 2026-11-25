# Restart files of Reactant simulations are written on the host (arrays, dates and numbers instead of Reactant's
# types), so that they can be read back in another process and on any architecture. Stored as Reactant types,
# `StartFromFile` failed with "AssertionError: buffer.buffer !== C_NULL": the device buffers are gone.

# all arrays of a (nested) NamedTuple of variables on the host, with their path, for comparison
host_arrays(x::AbstractArray, path) = [path => SpeedyWeather.on_architecture(SpeedyWeather.CPU(), x)]
host_arrays(x::NamedTuple, path) = reduce(vcat, [host_arrays(getfield(x, key), "$path.$key") for key in keys(x)]; init = Pair{String, Any}[])
host_arrays(x, path) = Pair{String, Any}[]

# whether two simulations have the same prognostic variables and clock (time)
function same_prognostic_variables(variables, variables_reference)
    arrays = Dict(host_arrays(variables.prognostic, "prognostic"))
    arrays_reference = Dict(host_arrays(variables_reference.prognostic, "prognostic"))
    same_arrays = keys(arrays) == keys(arrays_reference) && all(arrays[key] == arrays_reference[key] for key in keys(arrays))
    same_time = convert(DateTime, variables.prognostic.clock.time) == convert(DateTime, variables_reference.prognostic.clock.time)
    return same_arrays && same_time
end

restart_spectral_grid(architecture) = SpectralGrid(; architecture, truncation = TRUNCATION, nlayers = 8)

@testset "Restart files with Reactant" begin
    path = mktempdir(pwd(), prefix = "tmp_reactant_restart_")

    # a Reactant simulation with a state that differs from the initial conditions, random vorticity
    # without the extra degree, which initialize! removes from the prognostic variables of the restart
    model = PrimitiveWetModel(restart_spectral_grid(ReactantDevice()); feedback = nothing)
    simulation = initialize!(model)
    spectrum = restart_spectral_grid(SpeedyWeather.CPU()).spectrum
    vorticity = rand(ComplexF32, spectrum, size(simulation.variables.prognostic.vorticity)[2:end]...)
    LowerTriangularArrays.truncate_extra_degree!(vorticity)
    copyto!(simulation.variables.prognostic.vorticity, SpeedyWeather.on_architecture(ReactantDevice(), vorticity.data))
    simulation.variables.prognostic.clock.time = DateTime(2001, 3, 4, 5, 6)

    restart = WriteVariablesRestartFile(path = joinpath(path, "run_0001"), write_only_with_output = false)
    SpeedyWeather.finalize!(restart, simulation.variables, model)

    # written on the host: no Reactant types in the file
    SpeedyWeather.jldopen(joinpath(path, "run_0001", "restart.jld2")) do file
        prognostic = file["variables.prognostic"]
        @test prognostic.clock.time isa DateTime
        @test prognostic.clock.Δt isa Millisecond
        @test prognostic.vorticity.data isa Array
        @test prognostic.scale isa Base.RefValue{<:AbstractFloat}
        @test !occursin("Concrete", string(typeof(prognostic)))
    end

    # read back on Reactant and on the CPU
    for architecture in (ReactantDevice(), SpeedyWeather.CPU())
        model_restarted = PrimitiveWetModel(restart_spectral_grid(architecture); initial_conditions = StartFromFile(; path), feedback = nothing)
        simulation_restarted = initialize!(model_restarted)
        @test same_prognostic_variables(simulation_restarted.variables, simulation.variables)
    end

    # and a restart file written on the CPU read on Reactant
    model_cpu = PrimitiveWetModel(restart_spectral_grid(SpeedyWeather.CPU()); initial_conditions = StartFromFile(; path), feedback = nothing)
    simulation_cpu = initialize!(model_cpu)
    simulation_cpu.variables.prognostic.clock.time = DateTime(2002, 7, 8)
    restart_cpu = WriteVariablesRestartFile(path = joinpath(path, "run_0002"), write_only_with_output = false)
    SpeedyWeather.finalize!(restart_cpu, simulation_cpu.variables, model_cpu)

    model_restarted = PrimitiveWetModel(restart_spectral_grid(ReactantDevice()); initial_conditions = StartFromFile(; path, run_number = 2), feedback = nothing)
    simulation_restarted = initialize!(model_restarted)
    @test same_prognostic_variables(simulation_restarted.variables, simulation_cpu.variables)
end
