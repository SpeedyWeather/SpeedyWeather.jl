abstract type AbstractParticleAdvection <: AbstractModelComponent end

# function barrier for all particle advections, dispatch by model.particle_advection
# 1. initial conditions for particles
initialize!(particles::AbstractVector{P}, vars, model) where {P <: Particle} =
    initialize!(particles, vars, model.particle_advection, model)

# 2. initialize the particle advection work arrays
function initialize!(
        vars::Variables,
        particles::AbstractVector{P},   # for dispatch to distinguish from other initialize! functions
        model::AbstractModel,
    ) where {P <: Particle}
    # dispatch by model.particle_advection
    return initialize!(vars, particles, model.particle_advection, model)
end

# 3. the repeated call to actually advect particles
particle_advection!(vars, model) = particle_advection!(vars, model.particle_advection, model)

# no particle advection
particle_advection!(vars, ::Nothing, ::AbstractModel) = nothing

export ParticleAdvection2D

export ParticleAdvection2D
@kwdef struct ParticleAdvection2D{
        NF,
        GeometryType, # <: AbstractGridGeometry
    } <: AbstractParticleAdvection
    "[OPTION] Number of particles"
    nparticles::Int = 10

    "[OPTION] Execute particle advection every n timesteps"
    every_n_timesteps::Int = 6

    "[OPTION] Advect with velocities from this vertical layer index"
    layer::Int = 1

    "[DERIVED] Time step used for particle advection (scaled by radius, converted to degrees) [s*˚/m]"
    Δt::Base.RefValue{NF} = Ref(zero(NF))

    "[DERIVED] Interpolation geometry used during advection"
    geometry::GeometryType
end

function ParticleAdvection2D(SG::SpectralGrid; kwargs...)
    geometry = GridGeometry(SG.grid; NF = SG.NF)
    return ParticleAdvection2D{SG.NF, typeof(geometry)}(; geometry, kwargs...)
end

function variables(P::ParticleAdvection2D)
    (; nparticles) = P
    return (
        PrognosticVariable(:particles, ParticleVectorDim(nparticles), desc = "Particle locations", units = "˚/1"),
        ParticleVariable(:locations, ParticleVectorDim(nparticles), desc = "Particle locations", units = "˚/1"),
        ParticleVariable(:u, VectorDim(nparticles), desc = "Zonal velocity at particle location", units = "m/s"),
        ParticleVariable(:v, VectorDim(nparticles), desc = "Meridional velocity at particle location", units = "m/s"),
        ParticleVariable(:w, VectorDim(nparticles), desc = "Vertical velocity in σ coordinate at particle location", units = "1"),
        ParticleVariable(:locator, LocatorDim(nparticles), desc = "Vertical σ coordinate of particle location", units = "1"),
    )
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
    particle_advection.Δt[] *= (360 / 2π)
    return particle_advection
end

"""
$(TYPEDSIGNATURES)
Initialize particle locations uniformly in latitude, longitude and in the
vertical σ coordinates. This uses a cosin-distribution in latitude for
an equal-area uniformity."""
function initialize!(
        particles::AbstractVector{P},
        vars::Variables,
        particle_advection::ParticleAdvection2D,
        model::AbstractModel,
    ) where {P <: Particle}
    particles .= on_architecture(architecture(particles), rand(P, length(particles)))
    return particles
end

"""$(TYPEDSIGNATURES)
Initialize particle advection time integration: Store u,v interpolated initial conditions
in `diagn.particles.u` and `.v`  to be used when particle advection actually executed for first time."""
function initialize!(
        vars::Variables,
        particles::AbstractVector{P},
        particle_advection::ParticleAdvection2D,
        model::AbstractModel,
    ) where {P <: Particle}

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    k = particle_advection.layer
    u_grid = field_view(vars.grid.u, :, k)
    v_grid = field_view(vars.grid.v, :, k)
    (; locator) = vars.particles
    (; geometry) = particle_advection

    # interpolate initial velocity on initial locations
    lats = vars.particles.u    # reuse u,v arrays as only used for u, v
    lons = vars.particles.v    # after update_locator!
    σ = Vector(model.geometry.σ_levels_full)[k] # explicitly on CPU

    # modulo particles and extract their coordinates
    launch!(
        architecture(particles), LinearWorkOrder, (length(particles),), _initialize_particles_kernel!,
        particles, lons, lats, σ
    )

    RingGrids.update_locator!(locator, geometry, lons, lats)
    u0 = vars.particles.u      # now reused arrays are actually u, v
    v0 = vars.particles.v
    interpolate!(u0, u_grid, locator, geometry)
    interpolate!(v0, v_grid, locator, geometry)
    return nothing
end

# Kernel to modulo particles and extract their coordinates
@kernel inbounds = true function _initialize_particles_kernel!(
        particles, lons, lats, σ
    )
    i = @index(Global, Linear)

    # modulo all particles here
    # i.e. one can start with a particle at -120˚E which moduloed to 240˚E here
    particles[i] = mod(set(particles[i]; σ = σ))
    lons[i] = particles[i].lon
    lats[i] = particles[i].lat
end

# function barrier, unpack what's needed
function particle_advection!(vars, adv::ParticleAdvection2D, model::AbstractModel)
    return particle_advection!(vars.prognostic.particles, vars, vars.prognostic.clock, adv)
end

function particle_advection!(
        particles::AbstractVector{P},
        vars::Variables,
        clock::Clock,
        particle_advection::ParticleAdvection2D,
    ) where {P <: Particle}

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    (; locator) = vars.particles
    (; geometry) = particle_advection

    # decide whether to execute on this time step:
    # execute always on last time step *before* time step is divisible by
    # `particle_advection.every_n_timesteps`, e.g. 7, 15, 23, ... for n=8 which
    # already contains u, v at i=8, 16, 24, etc as executed after `transform!`
    # even though the clock hasn't be step forward yet, this means time = time + Δt here

    # should not be called on the 1st step in first_timesteps, which is excluded
    # with a lf2 == 2 check before this function is called

    # escape immediately if advection not on this timestep
    n = particle_advection.every_n_timesteps
    clock.timestep_counter % n == (n - 1) || return nothing

    # HEUN: PREDICTOR STEP, use u, v at previous time step and location
    Δt = particle_advection.Δt[]        # time step [s*˚/m]
    Δt_half = Δt / 2                      # /2 because Heun is average of Euler+corrected step

    u_old = vars.particles.u           # from previous time step and location
    v_old = vars.particles.v           # from previous time step and location

    # HACK: reuse u, v arrays (old velocity) on the fly for interpolation
    # as they're not needed anymore after new (predicted) location is found
    # same is true for the corrector step, interpolating velocities for the
    # next time step of the particle advection
    lons = vars.particles.u
    lats = vars.particles.v

    # Launch predictor step kernel
    launch!(
        architecture(u_old), LinearWorkOrder, (length(particles),),
        predictor_step_kernel!, particles, vars.particles.locations, u_old, v_old, lons, lats, Δt_half
    )

    # CORRECTOR STEP, use u, v at new location and new time step
    k = particle_advection.layer
    u_grid = field_view(vars.grid.u, :, k)
    v_grid = field_view(vars.grid.v, :, k)
    RingGrids.update_locator!(locator, geometry, lons, lats)

    # interpolate new velocity on predicted new locations
    u_new = vars.particles.u
    v_new = vars.particles.v
    interpolate!(u_new, u_grid, locator, geometry)
    interpolate!(v_new, v_grid, locator, geometry)

    # Launch corrector step kernel
    launch!(
        architecture(u_new), LinearWorkOrder, (length(particles),),
        corrector_step_kernel!, particles, u_new, v_new, lons, lats, Δt_half
    )

    # store new velocities at new (corrected locations) to be used on
    # next particle advection time step
    RingGrids.update_locator!(locator, geometry, lons, lats)
    interpolate!(u_new, u_grid, locator, geometry)
    interpolate!(v_new, v_grid, locator, geometry)
    return nothing
end

@inline function advect_2D(
        particle::Particle{NF},         # particle to advect
        u::NF,                          # zonal velocity [m/s]
        v::NF,                          # meridional velocity [m/s]
        dt::NF,                         # scaled time step [s*˚/m]
    ) where {NF}

    dlat = v * dt                               # increment in latitude [˚N]
    coslat = max(cosd(particle.lat), eps(NF))   # prevents division by zero
    dlon = u * dt / coslat                        # increment in longitude [˚E]
    return mod(move(particle, dlon, dlat))      # move, mod back to [0, 360˚E], [-90, 90˚N]
end

# Kernel for predictor step in Heun's method
@kernel inbounds = true function predictor_step_kernel!(
        particles, locations, u_old, v_old, lons, lats, @Const(Δt_half)
    )
    i = @index(Global, Linear)

    if isactive(particles[i])
        # sum up Heun's first term in 1/2*Δt*(uv_old + uv_new) on the fly
        # use only Δt/2
        particles[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)

        # predictor step, used to evaluate u_new, v_new
        # now again with Δt/2 to have an Euler timestep with Δt together with prev line
        locations[i] = advect_2D(particles[i], u_old[i], v_old[i], Δt_half)

        # reuse work arrays on the fly for new (predicted) locations
        lons[i] = locations[i].lon
        lats[i] = locations[i].lat
    end
end

# Kernel for corrector step in Heun's method
@kernel inbounds = true function corrector_step_kernel!(
        particles, u_new, v_new, lons, lats, @Const(Δt_half)
    )
    i = @index(Global, Linear)

    if isactive(particles[i])
        # sum up Heun's 2nd term in 1/2*Δt*(uv_old + uv_new) on the fly
        particles[i] = advect_2D(particles[i], u_new[i], v_new[i], Δt_half)

        # reuse work arrays on the fly for new (correct) locations
        lons[i] = particles[i].lon
        lats[i] = particles[i].lat
    end
end
