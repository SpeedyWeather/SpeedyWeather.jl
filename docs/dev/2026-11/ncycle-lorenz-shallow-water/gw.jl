using SpeedyWeather, Random, Printf
mutable struct EtaVar <: SpeedyWeather.AbstractCallback
    t::Vector{Float64}; v::Vector{Float64}
end
EtaVar() = EtaVar(Float64[], Float64[])
function etavar(vars)
    η = get_step(vars.prognostic.η, 1)
    s = 0.0
    for (lm, (l, m)) in enumerate(zip(η.spectrum.l_indices, η.spectrum.m_indices))
        l == 1 && continue    # skip mean (1-based l)
        s += (m == 1 ? 1 : 2) * abs2(η[lm])
    end
    s
end
SpeedyWeather.initialize!(c::EtaVar, vars, model) = (empty!(c.t); empty!(c.v); push!(c.t, 0); push!(c.v, etavar(vars)))
function SpeedyWeather.callback!(c::EtaVar, vars, model)
    push!(c.t, vars.prognostic.clock.time_step_counter*model.time_stepping.Δt/3600)
    push!(c.v, etavar(vars))
end
SpeedyWeather.finalize!(::EtaVar, args...) = nothing

function rungw(; θ=nothing, μ=nothing, explicit=false, stepper=:leapfrog, Δt=Minute(40), centering=0.5, steps=4, variant=NCycleLorenzA(), days=2, trunc=128, kw...)
    Random.seed!(1234)
    spectral_grid = SpectralGrid(truncation=trunc, nlayers=1)
    implicit = explicit ? nothing : θ !== nothing ? ImplicitSWθ{Float32}(θ) : μ === nothing ? ImplicitShallowWater(spectral_grid; centering) : ImplicitSWX{Float32}(centering, μ)
    orography = EarthOrography(spectral_grid, smoothing=false)
    initial_conditions = RandomWaves(spectral_grid, lmin=10, lmax=30)
    time_stepping = stepper == :leapfrog ? Leapfrog(spectral_grid; Δt_at_T32=Δt, adjust_with_output=false, kw...) :
        NCycleLorenz(spectral_grid; Δt_at_T32=Δt, steps, variant, adjust_with_output=false)
    model = ShallowWaterModel(spectral_grid; orography, initial_conditions, implicit, time_stepping)
    cb = EtaVar(); add!(model, :etavar => cb)
    simulation = initialize!(model)
    run!(simulation, period=Day(days))
    r = cb.v ./ cb.v[1]
    blown = !all(isfinite, r) || maximum(r) > 10
    return (; ratio_end = blown ? NaN : r[end], rmax = maximum(filter(isfinite, r)), t=cb.t, r)
end
include("impl.jl")
