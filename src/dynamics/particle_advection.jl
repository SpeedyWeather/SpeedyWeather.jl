abstract type AbstractParticleAdvection <: AbstractModelComponent end

# function barrier for all particle advections, dispatch by model.particle_advection
# 1. initial conditions for particles
initialize!(particles::AbstractVector{P}, progn, diagn, model) where {P<:Particle}=
    initialize!(particles, progn, diagn, model.particle_advection, model)

# 2. initialize the particle advection work arrays
function initialize!(
    diagn::DiagnosticVariables,
    particles::AbstractVector{P},   # for dispatch to distinguish from other initialize! functions
    progn::PrognosticVariables,
    model::AbstractModel,
) where {P <: Particle}
    # dispatch by model.particle_advection
    initialize!(diagn, particles, progn, model.particle_advection, model)
end

# 3. the repeated call to actually advect particles
particle_advection!(progn, diagn, model) = particle_advection!(progn, diagn, model.particle_advection, model)

# no particle advection
particle_advection!(progn, diagn, ::Nothing, ::AbstractModel) = nothing

export ParticleAdvection2D
@kwdef struct ParticleAdvection2D{NF} <: AbstractParticleAdvection
    "[OPTION] Execute particle advection every n timesteps"
    every_n_timesteps::Int = 6

    "[OPTION] Advect with velocities from this vertical layer index"
    layer::Int = 1

    "Time step used for particle advection (scaled by radius, converted to degrees) [s*˚/m]"
    Δt::Base.RefValue{NF} = Ref(zero(NF))
end

function ParticleAdvection2D(SG::SpectralGrid; kwargs...)
    SG.nparticles == 0 && @warn "ParticleAdvection2D created but nparticles = 0 in spectral grid."
    ParticleAdvection2D{SG.NF}(; kwargs...)
end

function initialize!(
    particle_advection::ParticleAdvection2D,
    model::AbstractModel,
)
    (; nlayers) = model.spectral_grid
    (; layer) = particle_advection
    nlayers < layer && @warn "Particle advection on layer $layer on spectral grid with nlayers=$nlayers."

    (; every_n_timesteps) = particle_advection
    # Δt [˚*s/m] is scaled by radius to convert more easily from velocity [m/s]
    # to [˚/s] for particle locations in degree
    particle_advection.Δt[] = every_n_timesteps * model.time_stepping.Δt
    particle_advection.Δt[] *= (360/2π)
end

"""
$(TYPEDSIGNATURES)
Initialize particle locations uniformly in latitude, longitude and in the
vertical σ coordinates. This uses a cosin-distribution in latitude for
an equal-area uniformity."""
function initialize!(
    particles::AbstractVector{P},
    progn::PrognosticVariables,     # used for dispatch as all sub components
    diagn::DiagnosticVariables,     # have this function signature
    particle_advection::ParticleAdvection2D,
    model::AbstractModel,
) where {P <: Particle}
    for i in eachindex(particles)
        # uniform random in lon (360*rand), lat (cos-distribution), σ (rand)
        particles[i] = rand(P)
    end
end

"""$(TYPEDSIGNATURES)
Initialize particle advection time integration: Store u,v interpolated initial conditions
in `diagn.particles.u` and `.v`  to be used when particle advection actually executed for first time."""
function initialize!(
    diagn::DiagnosticVariables,
    particles::AbstractVector{P},
    progn::PrognosticVariables,
    particle_advection::ParticleAdvection2D,
    model::AbstractModel,
) where {P<:Particle}

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    k = particle_advection.layer
    u_grid = field_view(diagn.grid.u_grid, :, k)
    v_grid = field_view(diagn.grid.v_grid, :, k)
    (; interpolator) = diagn.particles
    
    # interpolate initial velocity on initial locations
    lats = diagn.particles.u    # reuse u,v arrays as only used for u, v
    lons = diagn.particles.v    # after update_locator!
    
    for i in eachindex(particles)
        # modulo all particles here
        # i.e. one can start with a particle at -120˚E which moduloed to 240˚E here
        particles[i] = mod(particles[i])
        lons[i] = particles[i].lon
        lats[i] = particles[i].lat
    end

    RingGrids.update_locator!(interpolator, lons, lats)
    u0 = diagn.particles.u      # now reused arrays are actually u, v
    v0 = diagn.particles.v
    interpolate!(u0, u_grid, interpolator)
    interpolate!(v0, v_grid, interpolator)
end

# function barrier, unpack what's needed
function particle_advection!(progn, diagn, adv::ParticleAdvection2D, model::AbstractModel)
    particle_advection!(progn.particles, diagn, progn.clock, adv)
end

function particle_advection!(
    particles::AbstractVector{P},
    diagn::AbstractVariables,
    clock::Clock,
    particle_advection::ParticleAdvection2D,
) where {P<:Particle}

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    # decide whether to execute on this time step:
    # execute always on last time step *before* time step is divisible by
    # `particle_advection.every_n_timesteps`, e.g. 7, 15, 23, ... for n=8 which
    # already contains u, v at i=8, 16, 24, etc as executed after `transform!`
    # even though the clock hasn't be step forward yet, this means time = time + Δt here

    # should not be called on the 1st step in first_timesteps, which is excluded
    # with a lf2 == 2 check before this function is called
        
    # escape immediately if advection not on this timestep
    n = particle_advection.every_n_timesteps
    clock.timestep_counter % n == (n-1) || return nothing   

    # also escape if no particle is active
    any(isactive.(particles)) || return nothing

    # HEUN: PREDICTOR STEP, use u, v at previous time step and location
    Δt = particle_advection.Δt[]        # time step [s*˚/m]
    Δt_half = Δt/2                      # /2 because Heun is average of Euler+corrected step

    u_old = diagn.particles.u           # from previous time step and location
    v_old = diagn.particles.v           # from previous time step and location
    
    # HACK: reuse u, v arrays (old velocity) on the fly for interpolation
    # as they're not needed anymore after new (predicted) location is found
    # same is true for the corrector step, interpolating velocities for the
    # next time step of the particle advection
    lons = diagn.particles.u
    lats = diagn.particles.v

    for i in eachindex(particles, u_old, v_old)
        isactive(particles[i]) || continue
        # sum up Heun's first term in 1/2*Δt*(uv_old + uv_new) on the fly
        # use only Δt/2
        particles[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)
        
        # predictor step, used to evaluate u_new, v_new
        # now again with Δt/2 to have an Euler timestep with Δt together with prev line
        diagn.particles.locations[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)

        # reuse work arrays on the fly for new (predicted) locations
        lons[i] = diagn.particles.locations[i].lon
        lats[i] = diagn.particles.locations[i].lat
    end

    # CORRECTOR STEP, use u, v at new location and new time step
    k = particle_advection.layer
    u_grid = field_view(diagn.grid.u_grid, :, k)
    v_grid = field_view(diagn.grid.v_grid, :, k)
    (; interpolator) = diagn.particles
    RingGrids.update_locator!(interpolator, lons, lats)

    # interpolate new velocity on predicted new locations
    u_new = diagn.particles.u
    v_new = diagn.particles.v
    interpolate!(u_new, u_grid, interpolator)
    interpolate!(v_new, v_grid, interpolator)

    for i in eachindex(particles, u_new, v_new)
        isactive(particles[i]) || continue

        # sum up Heun's 2nd term in 1/2*Δt*(uv_old + uv_new) on the fly
        particles[i] = advect_2D(particles[i], u_new[i], v_new[i], Δt_half)

        # reuse work arrays on the fly for new (correct) locations
        lons[i] = particles[i].lon
        lats[i] = particles[i].lat
    end

    # store new velocities at new (corrected locations) to be used on
    # next particle advection time step
    RingGrids.update_locator!(interpolator, lons, lats)
    interpolate!(u_new, u_grid, interpolator)
    interpolate!(v_new, v_grid, interpolator)
    return nothing
end

function advect_2D(
    particle::Particle{NF},         # particle to advect
    u::NF,                          # zonal velocity [m/s]
    v::NF,                          # meridional velocity [m/s]
    dt::NF,                         # scaled time step [s*˚/m]    
) where NF

    dlat = v * dt                               # increment in latitude [˚N]
    coslat = max(cosd(particle.lat), eps(NF))   # prevents division by zero
    dlon = u * dt/coslat                        # increment in longitude [˚E]
    return mod(move(particle, dlon, dlat))      # move, mod back to [0, 360˚E], [-90, 90˚N]
end