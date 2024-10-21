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

@kwdef struct SpectralAR1Process{NF} <: AbstractRandomProcess
    time_scales::Vector{Second} = [Hour(6), Day(3), Day(30)]
    wavenumbers::Vector{Int} = [26, 16, 8]
    standard_deviations::Vector{NF} = [0.52, 0.18, 0.06]
    clamp::NTuple{2, NF} = (-1, 1)
    rand_function::Function = randn

    seed::Int = 123
    random_number_generator::Random.Xoshiro = Random.Xoshiro(seed)

    autoregressive_factors::Vector{NF} = zeros(NF, length(time_scales))
    noise_factors::Vector{NF} = zeros(NF, length(time_scales))
end

SpectralAR1Process(SG::SpectralGrid, kwargs...) = SpectralAR1Process{SG.NF}(; kwargs...)

function initialize!(
    process::SpectralAR1Process,
    model::AbstractModel,
)
    @assert maximum(process.wavenumbers) <= model.spectral_grid.trunc

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
            r = rand_function(process.random_number_generator, Complex{NF})
            random_pattern[l+1, m+1] = a*random_pattern[l+1, m+1] + ξ*r
        end
    end
end