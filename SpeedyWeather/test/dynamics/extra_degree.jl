# All spectral variables carry one more degree than the truncation (l = lmax, 1-based) as vector
# quantities need it for the meridional derivative. The prognostic variables are scalars that must
# not use it: once non-zero it would never change again (all tendencies are zero there) but still
# enter the grid-point fields every time step. This callback records the largest absolute value
# in that extra degree of all spectral prognostic variables (incl. tracers) after every time step.
mutable struct ExtraDegreeLog <: SpeedyWeather.AbstractCallback
    max_extra_degree::Vector{Float64}
end

ExtraDegreeLog() = ExtraDegreeLog(Float64[])

# largest absolute value in the extra degree of all spectral variables in a (nested) NamedTuple
max_extra_degree(vars::NamedTuple) = maximum(max_extra_degree, values(vars); init = 0.0)
max_extra_degree(::Any) = 0.0       # grid variables, clock, ...

function max_extra_degree(L::LowerTriangularArray)
    (; lmax, l_indices) = L.spectrum
    extra = selectdim(L.data, 1, l_indices .== lmax)
    return Float64(maximum(abs, extra; init = 0))
end

# initial conditions (incl. any set! after initialize!) before the first time step
SpeedyWeather.initialize!(log::ExtraDegreeLog, vars, model) =
    (log.max_extra_degree = [max_extra_degree(vars.prognostic)]; nothing)
SpeedyWeather.callback!(log::ExtraDegreeLog, vars, model) =
    (push!(log.max_extra_degree, max_extra_degree(vars.prognostic)); nothing)
SpeedyWeather.finalize!(::ExtraDegreeLog, args...) = nothing

@testset "No extra degree in prognostic variables" begin
    @testset for Model in (BarotropicModel, ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel)
        nlayers = Model <: Union{BarotropicModel, ShallowWaterModel} ? 1 : 4
        spectral_grid = SpectralGrid(truncation = 32; nlayers)
        model = Model(spectral_grid)
        add!(model, Tracer(:abc))

        log = ExtraDegreeLog()
        add!(model, log)
        simulation = initialize!(model)

        # a tracer with a sharp edge, which has a large projection onto the extra degree
        set!(simulation, abc = (λ, φ, σ) -> φ > 0 ? 1 : 0, namespace = :tracers)

        run!(simulation, period = Day(1))
        @test !simulation.model.feedback.nans_detected
        @test length(log.max_extra_degree) == simulation.variables.prognostic.clock.n_steps + 1
        @test all(iszero, log.max_extra_degree)
    end
end
