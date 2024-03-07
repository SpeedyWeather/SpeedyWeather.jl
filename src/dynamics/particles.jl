abstract type AbstractParticle{NF} end

export Particle
struct Particle{
    NF<:AbstractFloat,  # number format of coordinates
    isactive,           # bool indicating whether particle is active
} <: AbstractParticle{NF}
    lon::NF             # longitude [0,360˚]
    lat::NF             # latitude [-90˚,90˚]
    σ::NF               # vertical sigma coordinate [0,1] (top to surface)
end

# particle is by default active, of number format DEFAULT_NF and at 0˚N, 0˚E, σ=0

# keyword constructors
Particle{NF}(;lon,lat,σ=0) where NF = Particle{NF}(lon,lat,σ)
Particle{NF,isactive}(;lon,lat,σ=0) where {NF,isactive} = Particle{NF,isactive}(lon,lat,σ)
Particle(;lon,lat,σ=0) = Particle(lon,lat,σ)
Particle(;lon,lat,σ=0) = Particle(lon,lat,σ)

# empty constructor map to zero
Particle() = zero(Particle)
Particle{NF}() where NF = zero(Particle{NF})
Particle{NF,isactive}() where {NF,isactive} = zero(Particle{NF,isactive})
Particle{NF}(args...) where NF = Particle{NF,true}(args...)

# promotion of arguments
Particle(lon,lat) = Particle(lon,lat,0)
Particle(lon::Integer,lat::Integer) = Particle(lon,lat,0)
Particle(args::Integer...) = Particle{DEFAULT_NF,true}(args...)
Particle(args...) = Particle{promote_type(typeof.(args)...),true}(args...)

# zero generators
Base.zero(::Type{Particle}) = Particle{DEFAULT_NF,true}(0,0,0)
Base.zero(::Type{Particle{NF}}) where NF = Particle{NF,true}(0,0,0)
Base.zero(::Type{Particle{NF,isactive}}) where {NF,isactive} = Particle{NF,isactive}(0,0,0)
Base.zero(::P) where {P<:Particle} = zero(P)

Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle}) = rand(rng,Particle{DEFAULT_NF,true})
Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle{NF}}) where NF = rand(rng,Particle{NF,true})
function Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle{NF,isactive}}) where {NF,isactive}
    lon = 360*rand(rng,NF)
    yπ = convert(NF,π)*(2rand(rng,NF)-1)    # yπ ∈ [-1,1]*π
    lat = sign(yπ)*acosd((1 + cos(yπ))/2)   # cos-distributed latitude
    σ = rand(rng,NF)
    return Particle{NF,isactive}(lon,lat,σ)
end

# equality with same location and same activity
function Base.:(==)(p1::Particle{NF1,active1},p2::Particle{NF2,active2}) where {NF1,active1,NF2,active2}
    return (active1==active2) && (p1.lon == p2.lon) && (p1.lat == p2.lat) && (p1.σ == p2.σ)
end

function Base.show(io::IO,p::Particle{NF,isactive}) where {NF,isactive}
    lat = @sprintf("%6.2f",p.lat)
    lon = @sprintf("%6.2f",p.lon)
    σ = @sprintf("%.2f",p.σ)
    activity = isactive ? "  active" : "inactive"
    print(io,"Particle{$NF, $activity}($(lon)˚E, $(lat)˚N, σ=$σ)")
end

export move
@inline function move(p::Particle{NF,true},dlon,dlat,dσ) where NF
    (;lon, lat, σ) = p
    Particle{NF,true}(lon+dlon,lat+dlat,σ+dσ)
end

@inline function move(p::Particle{NF,true},dlon,dlat) where NF
    (;lon, lat, σ) = p
    Particle{NF,true}(lon+dlon,lat+dlat,σ)
end

# don't move inactive particles
@inline move(p::Particle{NF,false},args...) where NF = p

export activate, deactivate, active
activate(  p::Particle{NF}) where NF = Particle{NF, true}(p.lon,p.lat,p.σ)
deactivate(p::Particle{NF}) where NF = Particle{NF,false}(p.lon,p.lat,p.σ)
active(::Particle{NF,isactive}) where {NF,isactive} = isactive

function Base.mod(p::P) where {P<:Particle}
    (;lon, lat, σ) = p
    crossed_pole = lat > 90 || lat < -90
    lat = lat - 2crossed_pole*((lat + 90) % 180)
    lon = mod(lon + 180*crossed_pole,360)
    σ = clamp(σ,0,1)
    return P(lon,lat,σ)
end

abstract type AbstractParticleAdvection <: AbstractModelComponent end

export NoParticleAdvection
struct NoParticleAdvection <: AbstractParticleAdvection end
n_particles(::NoParticleAdvection) = 0
NoParticleAdvection(SG::SpectralGrid) = NoParticleAdvection()
particle_advection!(particles,diagn,dt,::NoParticleAdvection) = nothing

export ParticleAdvection
Base.@kwdef struct ParticleAdvection{
    NF<:AbstractFloat,
    Grid<:AbstractGrid,
} <: AbstractParticleAdvection

    "Number of latitude rings on one hemisphere, Equator incl. Resolution parameter."
    nlat_half::Int

    "[OPTION] Number of particles to advect"
    n_particles::Int = 10

    "Work array"
    a::Vector{NF} = zeros(n_particles)
    b::Vector{NF} = zeros(n_particles)

    "Interpolator to interpolate velocity fields onto particle positions"
    interpolator::AnvilInterpolator{NF,Grid} = AnvilInterpolator(NF, Grid, nlat_half, n_particles)
end

ParticleAdvection(SG::SpectralGrid;kwargs...) = ParticleAdvection{SG.NF,SG.Grid}(;nlat_half=SG.nlat_half,kwargs...)
n_particles(p::ParticleAdvection) = p.n_particles

function initialize!(
    particles::Vector{P},
    model::ModelSetup,
) where {P<:Particle}

    for i in eachindex(particles)
        particles[i] = rand(P)
    end
end

function particle_advection!(
    particles::Vector{Particle{NF}},
    diagn,
    dt::Real,
    particle_advection::ParticleAdvection,
) where NF

    # escape immediately for no particles
    particle_advection.n_particles == 0 && return nothing

    # also escape if no particle is active
    any_active::Bool = false 
    for particle in particles
        any_active |= active(particle)
    end
    any_active || return nothing

    (;u_grid, v_grid) = diagn.layers[1].grid_variables

    lats = particle_advection.a
    lons = particle_advection.b

    for i in eachindex(particles,lats,lons)
        lats[i] = particles[i].lat
        lons[i] = particles[i].lon
    end

    (;interpolator) = particle_advection
    RingGrids.update_locator!(interpolator,lats,lons)

    # effective time step in seconds * degree to move particles with coordinates in degrees
    # technically 360/2πr with radius r, but timestep dt is already scaled as dt/r
    degrees_per_meter = 360/2π
    dt = convert(NF,dt*degrees_per_meter)
    u_interp = particle_advection.a
    v_interp = particle_advection.b

    interpolate!(u_interp,u_grid,interpolator)
    interpolate!(v_interp,v_grid,interpolator)

    for i in eachindex(particles,u_interp,v_interp)
        particle = particles[i]
        u, v = u_interp[i], v_interp[i]
        particles[i] = advect_2D(particle,u,v,dt)
    end
end

function advect_2D(particle::Particle{NF,true},u::NF,v::NF,dt::NF) where NF
    dlat = v * dt
    dlon = u * dt/cosd(particle.lat)
    return mod(move(particle,dlon,dlat))
end

@inline advect_2D(p::Particle{NF,false},args...) where NF = p