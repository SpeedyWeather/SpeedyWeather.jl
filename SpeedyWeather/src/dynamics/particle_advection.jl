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

# TODO: remove fallback for reactant when it's compatible with particle advection
initialize!(::Variables, ::Nothing, ::AbstractModel) = nothing

# 3. the repeated call to actually advect particles
particle_advection!(vars, model) = particle_advection!(vars, model.particle_advection, model)

# no particle advection
particle_advection!(vars, ::Nothing, ::AbstractModel) = nothing

export ParticleAdvection2D
@kwdef struct ParticleAdvection2D{
        NF,
        GeometryType, # <: AbstractGridGeometry
        IntType,
    } <: AbstractParticleAdvection

    "[OPTION] Number of particles"
    nparticles::IntType = 10

    "[OPTION] Execute particle advection every n timesteps"
    every_n_time_steps::IntType = 6

    "[OPTION] Advect with velocities from this vertical layer index"
    layer::IntType = 1

    "[DERIVED] Interpolation geometry used during advection"
    geometry::GeometryType
end

function ParticleAdvection2D(SG::SpectralGrid; kwargs...)
    geometry = GridGeometry(SG.grid; NF = SG.NF)
    return ParticleAdvection2D{SG.NF, typeof(geometry), typeof(SG.truncation)}(; geometry, kwargs...)
end

export ParticleAdvection3D

# Continuous 3D particle advection: particles are tracked at their own σ coordinate.
# σ_levels_full is read from model.geometry at runtime.
@kwdef struct ParticleAdvection3D{
        NF,
        GeometryType, # <: AbstractGridGeometry
        IntType,
    } <: AbstractParticleAdvection

    "[OPTION] Number of particles"
    nparticles::IntType = 10

    "[OPTION] Execute particle advection every n timesteps"
    every_n_time_steps::IntType = 6

    "[DERIVED] Interpolation geometry used during advection"
    geometry::GeometryType
end

function ParticleAdvection3D(SG::SpectralGrid; kwargs...)
    geometry = GridGeometry(SG.grid; NF = SG.NF)
    return ParticleAdvection3D{SG.NF, typeof(geometry), typeof(SG.truncation)}(; geometry, kwargs...)
end

variables(P::AbstractParticleAdvection) = variables(typeof(P), P.nparticles)

# Common particle variables shared between 2D and 3D advection (no locator)
function _common_particle_variables(nparticles)
    return (
        PrognosticVariable(:particles, ParticleVectorDim(nparticles), desc = "Particle locations", units = "˚/1"),
        ParticleVariable(:locations, ParticleVectorDim(nparticles), desc = "Particle locations", units = "˚/1"),
        ParticleVariable(:u, VectorDim(nparticles), desc = "Zonal velocity at particle location", units = "m/s"),
        ParticleVariable(:v, VectorDim(nparticles), desc = "Meridional velocity at particle location", units = "m/s"),
    )
end

function variables(::Type{<:ParticleAdvection2D}, nparticles)
    return (
        _common_particle_variables(nparticles)...,
        ParticleVariable(:locator, LocatorDim(), desc = "Particle locator for horizontal interpolation", units = "1"),
    )
end

# 3D advection: same common variables but uses LocatorDim with nlayers set,
# embedding north/south_pole_average arrays in the locator (type-stable, avoids per-call allocation)
function variables(P::ParticleAdvection3D, model::AbstractModel)
    (; nparticles) = P
    (; nlayers) = model.spectral_grid
    return (
        _common_particle_variables(nparticles)...,
        ParticleVariable(:locator, LocatorDim(; nlayers), desc = "Particle locator with embedded pole averages for 3D interpolation", units = "1"),
        ParticleVariable(:w, VectorDim(nparticles), desc = "Vertical velocity dσ/dt at particle location", units = "1/s"),
    )
end

function initialize!(particle_advection::ParticleAdvection2D, model::AbstractModel)
    (; nlayers) = model.spectral_grid
    (; layer) = particle_advection
    nlayers < layer && @warn "Particle advection on layer $layer on spectral grid with nlayers=$nlayers."
    return nothing
end

initialize!(::ParticleAdvection3D, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Initialize particle locations uniformly in latitude, longitude and in the
vertical σ coordinates. This uses a cosin-distribution in latitude for
an equal-area uniformity."""
function initialize!(
        particles::AbstractVector{P},
        vars::Variables,
        particle_advection::AbstractParticleAdvection,
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
    (; time_stepping) = model

    # index the step dimension according to time stepper
    l = which_prognostic_step(vars.grid.u, time_stepping, particle_advection, model)
    u_grid = field_view(vars.grid.u, :, k, l)
    v_grid = field_view(vars.grid.v, :, k, l)
    (; locator) = vars.particles
    (; geometry) = particle_advection

    # interpolate initial velocity on initial locations
    lons = vars.particles.u                     # reuse u,v arrays as only used for u, v
    lats = vars.particles.v                     # after update_locator!
    σ = Vector(model.geometry.σ_levels_full)[k] # explicitly on CPU

    # modulo particles and extract their coordinates
    launch!(
        architecture(particles), LinearWorkOrder, (length(particles),), _initialize_2D_particles_kernel!,
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
# set every particle to vertical σ provided
@kernel inbounds = true function _initialize_2D_particles_kernel!(
        particles, lons, lats, σ
    )
    i = @index(Global, Linear)

    # modulo all particles here
    # i.e. one can start with a particle at -120˚E which moduloed to 240˚E here
    particles[i] = mod(set(particles[i]; σ = σ))
    lons[i] = particles[i].lon
    lats[i] = particles[i].lat
end

"""$(TYPEDSIGNATURES)
Initialize 3D particle advection work arrays: interpolate u, v, w at each particle's
initial 3D position to seed the Heun predictor for the first advection step."""
function initialize!(
        vars::Variables,
        particles::AbstractVector{P},
        particle_advection::ParticleAdvection3D,
        model::AbstractModel,
    ) where {P <: Particle}

    length(particles) == 0 && return nothing

    # index step dimension according to time stepper
    (; time_stepping) = model
    l = which_prognostic_step(vars.grid.u, time_stepping, particle_advection, model)
    u_3d = field_view(vars.grid.u, :, :, l)     # prognostic variables have a step dimension
    v_3d = field_view(vars.grid.v, :, :, l)
    w_3d = vars.dynamics.w                      # vertical velocity is diagnostic in hydrostatic models

    # interpolate initial velocity on initial locations
    σ_levels_full = model.geometry.σ_levels_full
    σ_levels_half = model.geometry.σ_levels_half
    (; locator) = vars.particles
    (; geometry) = particle_advection
    lons = vars.particles.u                     # reuse u,v arrays as only used for u, v
    lats = vars.particles.v                     # after update_locator!

    launch!(
        architecture(particles), LinearWorkOrder, (length(particles),),
        _initialize_3D_particles_kernel!, particles, lons, lats
    )

    RingGrids.update_locator!(locator, geometry, lons, lats)
    u0 = vars.particles.u
    v0 = vars.particles.v
    w0 = vars.particles.w

    interpolate_3D!(u0, u_3d, locator, geometry, particles, σ_levels_full, Center())
    interpolate_3D!(v0, v_3d, locator, geometry, particles, σ_levels_full, Center())
    interpolate_3D!(w0, w_3d, locator, geometry, particles, σ_levels_half, Face())
    return nothing
end

# Modulo particles and extract horizontal coordinates; σ is each particle's own
@kernel inbounds = true function _initialize_3D_particles_kernel!(
        particles, lons, lats
    )
    i = @index(Global, Linear)
    particles[i] = mod(particles[i])
    lons[i] = particles[i].lon
    lats[i] = particles[i].lat
end

function particle_advection!(
        vars::Variables,
        particle_advection::ParticleAdvection2D,
        model::AbstractModel,
    )

    (; particles, clock) = vars.prognostic

    # escape immediately for no particles
    length(particles) == 0 && return nothing

    (; locator) = vars.particles
    (; geometry) = particle_advection
    (; time_stepping) = model

    # decide whether to execute on this time step:
    # execute always on last time step *before* time step is divisible by
    # `particle_advection.every_n_time_steps`, e.g. 7, 15, 23, ... for n=8 which
    # already contains u, v at i=8, 16, 24, etc as executed after `transform!`
    # even though the clock hasn't be step forward yet, this means time = time + Δt here

    # should not be called on the 1st Euler step of Leapfrog, which is excluded
    # by using `clock.time_step_counter` (which doesn't count that spin up step)
    # as long as `every_n_time_steps` > 1

    # escape immediately if advection not on this timestep
    n = particle_advection.every_n_time_steps
    clock.time_step_counter % n == (n - 1) || return nothing
    NF = eltype(eltype(particles))      # number format used for particle locations
    (; radius) = model.planet

    # HEUN: PREDICTOR STEP, use u, v at previous time step and location
    # calculate time step on the fly
    Δt = model.time_stepping.Δt                 # time step [s]
    Δt *= n * convert(NF, 180 / (π * radius))   # scale to [s*°/m] to obtain [˚] when multiplied with velocity in [m/s]
    Δt_half = Δt / 2                    # /2 because Heun is average of Euler+corrected step

    u_old = vars.particles.u            # from previous time step and locationq
    v_old = vars.particles.v            # from previous time step and location

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
    l = which_prognostic_step(vars.grid.u, time_stepping, particle_advection, model)
    u_grid = field_view(vars.grid.u, :, k, l)
    v_grid = field_view(vars.grid.v, :, k, l)
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

function particle_advection!(
        vars::Variables,
        particle_advection::ParticleAdvection3D,
        model::AbstractModel,
    )
    (; particles, clock) = vars.prognostic
    length(particles) == 0 && return nothing

    (; locator) = vars.particles
    (; geometry) = particle_advection
    (; time_stepping) = model

    n = particle_advection.every_n_time_steps
    clock.time_step_counter % n == (n - 1) || return nothing
    NF = eltype(eltype(particles))
    (; radius) = model.planet

    # Calculate time step on the fly
    # horizontal time step, scale to [s*°/m] to obtain [˚] when multiplied with velocity in [m/s]
    Δt = model.time_stepping.Δt * n * convert(NF, 180 / (π * radius))
    Δt_half = Δt / 2            # /2 because Heun is average of Euler+corrected step

    # vertical time step without the degree/radius conversion for dσ = w[s⁻¹]*Δt_vert[s]
    Δt_vert = model.time_stepping.Δt * n
    Δt_vert_half = Δt_vert / 2

    u_old = vars.particles.u
    v_old = vars.particles.v
    w_old = vars.particles.w

    lons = vars.particles.u
    lats = vars.particles.v

    # HEUN: PREDICTOR STEP with old velocities and current particle position
    launch!(
        architecture(u_old), LinearWorkOrder, (length(particles),),
        predictor_step_3D_kernel!,
        particles, vars.particles.locations, u_old, v_old, w_old, lons, lats,
        Δt_half, Δt_vert_half,
    )

    # CORRECTOR STEP: interpolate new u, v, w at predicted locations
    l = which_prognostic_step(vars.grid.u, time_stepping, particle_advection, model)
    u_3d = field_view(vars.grid.u, :, :, l)
    v_3d = field_view(vars.grid.v, :, :, l)
    w_3d = vars.dynamics.w
    σ_levels_full = model.geometry.σ_levels_full
    σ_levels_half = model.geometry.σ_levels_half

    RingGrids.update_locator!(locator, geometry, lons, lats)
    u_new = vars.particles.u
    v_new = vars.particles.v
    w_new = vars.particles.w
    interpolate_3D!(u_new, u_3d, locator, geometry, vars.particles.locations, σ_levels_full, Center())
    interpolate_3D!(v_new, v_3d, locator, geometry, vars.particles.locations, σ_levels_full, Center())
    interpolate_3D!(w_new, w_3d, locator, geometry, vars.particles.locations, σ_levels_half, Face())

    launch!(
        architecture(u_new), LinearWorkOrder, (length(particles),),
        corrector_step_3D_kernel!,
        particles, u_new, v_new, w_new, lons, lats, Δt_half, Δt_vert_half,
    )

    # store new velocities at corrected position for next advection step
    RingGrids.update_locator!(locator, geometry, lons, lats)
    interpolate_3D!(u_new, u_3d, locator, geometry, particles, σ_levels_full, Center())
    interpolate_3D!(v_new, v_3d, locator, geometry, particles, σ_levels_full, Center())
    interpolate_3D!(w_new, w_3d, locator, geometry, particles, σ_levels_half, Face())
    return nothing
end

@inline function advect_2D(
        particle::Particle{NF},                 # particle to advect
        u::NF,                                  # zonal velocity [m/s]
        v::NF,                                  # meridional velocity [m/s]
        Δt::NF,                                 # scaled time step [s*˚/m]
    ) where {NF}

    dlat = v * Δt                                           # increment in latitude [˚N]
    coslat = max(cos(deg2rad(particle.lat)), eps(NF))       # prevents division by zero
    dlon = u * Δt / coslat                                  # increment in longitude [˚E]
    return mod(move(particle, dlon, dlat))      # move, mod back to [0, 360˚E], [-90, 90˚N]
end

@inline function advect_3D(
        particle::Particle{NF},
        u::NF,
        v::NF,
        w::NF,
        Δt::NF,         # horizontal: [s*˚/m]
        Δt_vert::NF,    # vertical: [s]; dσ = w[s⁻¹] * dt_vert[s]
    ) where {NF}

    dlat = v * Δt
    coslat = max(cos(deg2rad(particle.lat)), eps(NF))
    dlon = u * Δt / coslat
    dσ = w * Δt_vert
    return mod(move(particle, dlon, dlat, dσ))
end

# Kernel for predictor step in Heun's method
@kernel inbounds = true function predictor_step_kernel!(
        particles, locations, u_old, v_old, lons, lats, Δt_half
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
        particles, u_new, v_new, lons, lats, Δt_half
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

# 3D predictor: first Heun half-step with old velocity + predict full-step location
@kernel inbounds = true function predictor_step_3D_kernel!(
        particles, locations, u_old, v_old, w_old, lons, lats, Δt_half, Δt_vert_half
    )
    i = @index(Global, Linear)

    if isactive(particles[i])
        particles[i] = advect_3D(particles[i], u_old[i], v_old[i], w_old[i], Δt_half, Δt_vert_half)
        locations[i] = advect_3D(particles[i], u_old[i], v_old[i], w_old[i], Δt_half, Δt_vert_half)
        lons[i] = locations[i].lon
        lats[i] = locations[i].lat
    end
end

# 3D corrector: second Heun half-step with new velocity at predicted location
@kernel inbounds = true function corrector_step_3D_kernel!(
        particles, u_new, v_new, w_new, lons, lats, Δt_half, Δt_vert_half
    )
    i = @index(Global, Linear)

    if isactive(particles[i])
        particles[i] = advect_3D(particles[i], u_new[i], v_new[i], w_new[i], Δt_half, Δt_vert_half)
        lons[i] = particles[i].lon
        lats[i] = particles[i].lat
    end
end
