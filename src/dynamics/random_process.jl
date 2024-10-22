abstract type AbstractRandomProcess <: AbstractModelComponent end

"""$(TYPEDSIGNATURES)
General transform for `random processes <: AbstractRandomProcess`.
Takes the spectral `random_pattern` in the prognostic variables
and transforms it to spectral space in `diagn.grid.random_pattern`."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    random_process::AbstractRandomProcess,
    spectral_transform::SpectralTransform,
)
    grid = diagn.grid.random_pattern
    spec = progn.random_pattern
    transform!(grid, spec, spectral_transform)

    if :clamp in fieldnames(typeof(random_process))
        clamp!(grid, random_process.clamp...)
    end
end

export NoRandomProcess
"""Dummy type for no random process."""
struct NoRandomProcess <: AbstractRandomProcess end
NoRandomProcess(::SpectralGrid) = NoRandomProcess()

"""$(TYPEDSIGNATURES)
`NoRandomProcess` does not need to transform any random pattern from
spectral to grid space."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    random_process::NoRandomProcess,
    spectral_transform::SpectralTransform,
)
    return nothing
end

initialize!(process::NoRandomProcess, model::AbstractModel) = nothing
random_process!(progn::PrognosticVariables, process::NoRandomProcess) = nothing

export SpectralAR1Process

"""First-order auto-regressive random process (AR1) in spectral space,
evolving `wavenumbers` with respectice `time_scales` and `standard_deviations`
independently. Transformed after every time step to grid space with a
`clamp` applied to limit extrema. For reproducability `seed` can be
provided and an independent `random_number_generator` is used
that is reseeded on every `initialize!`. Fields are $(TYPEDFIELDS)"""
@kwdef struct SpectralAR1Process{NF} <: AbstractRandomProcess
    "[OPTION] Time scales of every AR1 process respectively"
    time_scales::Vector{Second} = [Hour(6 + 2h) for h in 0:10]

    "[OPTION] Wavenumbers of every AR1 process respectively"
    wavenumbers::Vector{Int} = [26-h for h in 0:10]

    "[OPTION] Standard deviations of every AR1 process respectively"
    standard_deviations::Vector{NF} = [0.4-h/20 for h in 0:10]

    "[OPTION] Range to clamp values into after every transform into grid space"
    clamp::NTuple{2, NF} = (-1, 1)

    "[OPTION] Function to be called for random number generation"
    rand_function::Function = randn

    "[OPTION] Random number generator seed"
    seed::Int = 123

    "Independent random number generator for this random process"
    random_number_generator::Random.Xoshiro = Random.Xoshiro(seed)

    "Precomputed auto-regressive factors [1], function of time scale"
    autoregressive_factors::Vector{NF} = zeros(NF, length(time_scales))

    "Precomputed noise factors [1], function of time scale"
    noise_factors::Vector{NF} = zeros(NF, length(time_scales))
end

# generator function
SpectralAR1Process(SG::SpectralGrid, kwargs...) = SpectralAR1Process{SG.NF}(; kwargs...)

function initialize!(
    process::SpectralAR1Process,
    model::AbstractModel,
)
    (; time_scales, wavenumbers, standard_deviations) = process
    @assert maximum(wavenumbers) <= model.spectral_grid.trunc
    @assert length(time_scales) == length(wavenumbers) == length(standard_deviations) || throw(DimensionMismatch)

    dt = model.time_stepping.Δt_sec         # in seconds

    for (i, (τ, σ)) in enumerate(zip(process.time_scales, process.standard_deviations))
        process.autoregressive_factors[i] = exp(-dt/Second(τ).value)
        process.noise_factors[i] = σ*sqrt(1 - exp(-dt/Second(τ).value))
    end

    # reseed the random number generator
    Random.seed!(process.random_number_generator, process.seed)
    return nothing
end

function random_process!(
    progn::PrognosticVariables,
    process::SpectralAR1Process{NF},
) where NF

    (; random_pattern) = progn
    (; wavenumbers, autoregressive_factors, noise_factors) = process
    (; rand_function) = process

    for (l, a, ξ) in zip(wavenumbers, autoregressive_factors,  noise_factors)
        for m in 0:l-1
            r = 2rand(process.random_number_generator, Complex{NF}) - (1 + im)
            random_pattern[l+1, m+1] = a*random_pattern[l+1, m+1] + ξ*r
        end
    end
end