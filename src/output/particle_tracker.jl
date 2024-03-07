export ParticleTracker
Base.@kwdef mutable struct ParticleTracker{NF} <: AbstractCallback
    "[OPTION] Frequency to track particles at"
    Δt::Second = Hour(1)

    "[OPTION] File name"
    file_name::String = "particles.jld2"
    
    "Number of particles to track"
    n_particles::Int = 0
    
    "Number of timesteps to track particles at"
    n_steps::Int = 0

    "Count up the tracked particle positions"
    counter::Int = 0

    "Execute callback every n time steps"
    every_n_timesteps::Int = 0

    # tracking arrays
    times::Vector{DateTime} = zeros(DateTime, n_steps+1)
    lons::Matrix{NF} = zeros(NF,n_particles, n_steps+1)
    lats::Matrix{NF} = zeros(NF,n_particles, n_steps+1)
    σs::Matrix{NF} = zeros(NF,n_particles, n_steps+1)
end

ParticleTracker(SG::SpectralGrid;kwargs...) = ParticleTracker{SG.NF}(;kwargs...)

function initialize!(
    callback::ParticleTracker{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
) where NF

    (;n_particles) = model.particle_advection
    callback.n_particles = n_particles

    # calculate tracking frequency, adjust tracking time step
    (;time_stepping) = model
    callback.every_n_timesteps = max(1,round(Int, 
                    Millisecond(callback.Δt).value/time_stepping.Δt_millisec.value))
    callback.Δt = Second(round(Int,callback.every_n_timesteps*time_stepping.Δt_sec))
    callback.n_steps = progn.clock.n_timesteps ÷ callback.every_n_timesteps + 1
    (;n_steps) = callback

    # allocate tracking arrays of correct size
    callback.lons  = zeros(NF, n_particles, n_steps+1)
    callback.lats  = zeros(NF, n_particles, n_steps+1)
    callback.σs    = zeros(NF, n_particles, n_steps+1)
    callback.times = zeros(DateTime, n_steps+1)

    # set initial conditions
    for (p,particle) in enumerate(progn.particles)
        callback.times[1]  = progn.clock.time
        callback.lons[p,1] = particle.lon
        callback.lats[p,1] = particle.lat
        callback.σs[p,1]   = particle.σ
    end
    
    # (re)set counter to 1
    callback.counter = 1
end

function callback!(
    callback::ParticleTracker{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
) where NF

    # call callback this timestep?
    progn.clock.timestep_counter % callback.every_n_timesteps == 0 || return nothing
    
    # otherwise call now, but escape if n_steps + 1 already reached as arrays are of that size
    callback.counter += 1  
    callback.counter > callback.n_steps + 1 && return nothing
    i = callback.counter

    for (p,particle) in enumerate(progn.particles)
        callback.times[i] = progn.clock.time
        callback.lons[p,i] = particle.lon
        callback.lats[p,i] = particle.lat
        callback.σs[p,i] = particle.σ
    end
end

function finish!(
    callback::ParticleTracker{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
) where NF

    (;lons, lats, σs, times, n_particles, file_name) = callback
    (;run_path) = model.output

    jldopen(joinpath(run_path,file_name),"w"; compress=true) do f
        f["description"] = "Particle advection simulated with SpeedyWeather.jl"
        f["version"] = model.output.pkg_version
        f["n_particles"] = n_particles
        f["lats"] = lats
        f["lons"] = lons
        f["sigma"] = σs
        f["times"] = times
    end

    # with no output run_path is "" by default
    run_path == "" && @info "ParticleTracker stored data in $(joinpath(pwd(),filename))"
end