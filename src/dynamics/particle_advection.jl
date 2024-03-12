abstract type AbstractParticleAdvection <: AbstractModelComponent end

export NoParticleAdvection
struct NoParticleAdvection <: AbstractParticleAdvection end
NoParticleAdvection(SG::SpectralGrid) = NoParticleAdvection()
initialize!(::NoParticleAdvection, ::ModelSetup) = nothing
particle_advection!(progn, diagn, lf, ::NoParticleAdvection) = nothing

export ParticleAdvection2D
Base.@kwdef struct ParticleAdvection2D{NF} <: AbstractParticleAdvection
    "[OPTION] Execute particle advection every n timesteps"
    every_n_timesteps::Int = 8

    "[OPTION] Advect with velocities from this vertical layer index"
    layer::Int = 1

    "Time step used for particle advection (scaled by radius, converted to degrees) [s*˚/m]"
    Δt::Base.RefValue{NF} = Ref(zero(NF))
end

function ParticleAdvection2D(SG::SpectralGrid; kwargs...)
    SG.n_particles == 0 && @warn "ParticleAdvection2D created but n_particles = 0 in spectral grid."
    ParticleAdvection2D{SG.NF}(; kwargs...)
end

function initialize!(
    particle_advection::ParticleAdvection2D,
    model::ModelSetup,
)
    (; every_n_timesteps) = particle_advection
    # Δt [˚*s/m] is scaled by radius to convert more easily from velocity [m/s]
    # to [˚/s] for particle locations in degree
    particle_advection.Δt[] = every_n_timesteps * model.time_stepping.Δt
    particle_advection.Δt[] *= (360/2π)
end

"""
$(TYPEDSIGNATURES)
Initialize particle locations uniformly in latitude, longitude and in the
vertical σ coordinates. This uses a cosin-distribution in latitude for
an equal-area uniformity."""
function initialize!(
    particles::Vector{P},
    model::ModelSetup,
) where {P<:Particle}
    for i in eachindex(particles)
        # uniform random in lon (360*rand), lat (cos-distribution), σ (rand)
        particles[i] = rand(P)
    end
end

# function barrier
function particle_advection!(progn, diagn, lf, adv::ParticleAdvection2D)
    particle_advection!(progn.particles, diagn, progn.clock, lf, adv)
end

function particle_advection!(
    particles::Vector{Particle{NF}},
    diagn::AbstractVariables,
    clock::Clock,
    lf::Int,    # leapfrog index used to distinguish between 1st and other steps
    particle_advection::ParticleAdvection2D,
) where NF

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    # escape immediately if advection not on this timestep
    clock.timestep_counter % particle_advection.every_n_timesteps == 0 || return nothing

    # don't do particle advection on the 2nd time step (which is half a leapfrog step that isn't counted)
    clock.timestep_counter == 0 && lf == 2 && return nothing

    # also escape if no particle is active
    any_active::Bool = false 
    for particle in particles
        any_active |= active(particle)
    end
    any_active || return nothing

    # HEUN: PREDICTOR STEP, use u, v at previous time step and location
    Δt = particle_advection.Δt[]        # time step [s*˚/m]

    # on the first time step the predictor step won't do anything because uv_old = 0
    # so the corrector step becomes Euler forward when we set Δt_half = Δt actually!
    Δt_half = Δt/lf                     # = Δt/2, but Δt on the 1st time step

    u_old = diagn.particles.u           # from previous time step and location
    v_old = diagn.particles.v           # from previous time step and location
    
    # HACK: reuse u, v arrays (old velocity) on the fly for interpolation
    # as they're not needed anymore after new (predicted) location is found
    # same is true for the corrector step, interpolating velocities for the
    # next time step of the particle advection
    lats = diagn.particles.u
    lons = diagn.particles.v

    for i in eachindex(particles, u_old, v_old)
        # sum up Heun's first term in 1/2*Δt*(uv_old + uv_new) on the fly
        # on first time step old u=v=0, so we just modulo all particles
        # so that one could start with a particle at -120˚E => 240˚E here 
        particles[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)
        
        # predictor step, used to evaluate u_new, v_new 
        diagn.particles.locations[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)

        # reuse work arrays on the fly for new (predicted) locations
        lats[i] = diagn.particles.locations[i].lat
        lons[i] = diagn.particles.locations[i].lon
    end

    # CORRECTOR STEP, use u, v at new location and new time step
    k = particle_advection.layer
    (;u_grid, v_grid) = diagn.layers[k].grid_variables
    (;interpolator) = diagn.particles
    RingGrids.update_locator!(interpolator, lats, lons)

    # interpolate new velocity on predicted new locations
    u_new = diagn.particles.u
    v_new = diagn.particles.v
    interpolate!(u_new, u_grid, interpolator)
    interpolate!(v_new, v_grid, interpolator)

    for i in eachindex(particles, u_new, v_new)
        # sum up Heun's 2nd term in 1/2*Δt*(uv_old + uv_new) on the fly
        particles[i] = advect_2D(particles[i], u_new[i], v_new[i], Δt_half)

        # reuse work arrays on the fly for new (correct) locations
        lats[i] = particles[i].lat
        lons[i] = particles[i].lon
    end

    # store new velocities at new (corrected locations) to be used on
    # next particle advection time step
    RingGrids.update_locator!(interpolator, lats, lons)
    interpolate!(u_new, u_grid, interpolator)
    interpolate!(v_new, v_grid, interpolator)
    return nothing
end

function advect_2D(
    particle::Particle{NF,true},    # particle to advect
    u::NF,                          # zonal velocity [m/s]
    v::NF,                          # meridional velocity [m/s]
    dt::NF,                         # scaled time step [s*˚/m]    
) where NF

    dlat = v * dt                               # increment in latitude [˚N]
    coslat = max(cosd(particle.lat), eps(NF))   # prevents division by zero
    dlon = u * dt/coslat                        # increment in longitude [˚E]
    return mod(move(particle,dlon,dlat))        # move, mod back to [0, 360˚E], [-90, 90˚N]
end

@inline advect_2D(p::Particle{NF, false}, args...) where NF = p