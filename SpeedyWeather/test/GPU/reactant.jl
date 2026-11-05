# Reactant on the GPU: the BarotropicModel steps without NaNs. Not a correctness test (too expensive
# here, see test/reactant for the correctness tests on Reactant's CPU backend). A standard model run,
# Reactant raises the kernels by default (raise = true) in `initialize!` (`@maybe_jit`) and `run!`

using Reactant, CUDA # we need CUDA even on other devices for the kernel raising

Reactant.set_default_backend("gpu")

@testset "Reactant GPU: BarotropicModel with $Transform steps without NaNs" for Transform in (MatrixSpectralTransform, SpectralTransform)
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice(), nlayers = 1)
    spectral_transform = Transform(spectral_grid)
    model = BarotropicModel(spectral_grid; spectral_transform, feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 10)

    @test !any(isnan, on_architecture(SpeedyWeather.CPU(), simulation.variables.prognostic.vorticity))
end

# `fill_fused_parents!`, `fused_members` and `variables_on_cpu` are defined in copy_variables.jl
@testset "Reactant GPU: deepcopy keeps every fused variable a view of its copied parent" begin
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice())
    vars = Variables(PrimitiveWetModel(spectral_grid; feedback = nothing))
    # mark the original's parents with 1000 + i before the deepcopy, so the copy starts with these values
    members = fused_members(fill_fused_parents!(vars, 1000), 1000)
    @test !isempty(members)

    # mark the copy's parents with 2000 + i: a view of a copied parent then holds 2000 + i, a detached
    # copy of a fused variable still the 1000 + i it was copied with
    copied = deepcopy(vars)
    fill_fused_parents!(copied, 2000)
    @test fused_members(copied, 2000) == members    # every fused variable of the copy is a view of its parent
    @test fused_members(vars, 1000) == members      # the original still holds 1000 + i: no shared memory
end

#= # TODO: Reactative test, the compiled_call just takes too long in the CI with consumer GPUs :(
# with both spectral transforms (the FFT-based one uses cuFFT through XLA on the GPU)
@testset "Reactant GPU: PrimitiveWetModel with $Transform steps without NaNs" for Transform in (SpectralTransform, ) # SpectralTransform
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice())
    spectral_transform = Transform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform, feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 10)

    (; prognostic) = simulation.variables
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
        @test !any(isnan, on_architecture(SpeedyWeather.CPU(), getproperty(prognostic, name)))
    end
end
=#

@testset "Reactant GPU: copy! between CPU and Reactant variables" begin
    model = PrimitiveWetModel(SpectralGrid(); feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 3)                     # non-trivial data in all groups
    vars = simulation.variables
    reference = variables_on_cpu(vars)

    reactant_model = PrimitiveWetModel(SpectralGrid(architecture = SpeedyWeather.ReactantDevice()); feedback = nothing)
    vars_reactant = copy!(Variables(reactant_model), vars)  # CPU → Reactant
    @test variables_on_cpu(vars_reactant) == reference
    vars_back = copy!(Variables(model), vars_reactant)      # Reactant → CPU
    @test variables_on_cpu(vars_back) == reference
end
