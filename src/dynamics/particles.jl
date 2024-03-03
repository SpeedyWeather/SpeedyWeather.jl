abstract type AbstractParticle{NF} end

export Particle
Base.@kwdef struct Particle{NF<:AbstractFloat} <: AbstractParticle{NF}
    lat::NF = 0         # [-90˚,90˚]
    lon::NF = 0         # [0,360˚]
    σ::NF   = 0         # [0,1]
    active::Bool = true # deactivated particles don't move
end

function Base.show(io::IO,p::Particle)
    lat = @sprintf("%6.2f",p.lat)
    lon = @sprintf("%7.2f",p.lon)
    σ = @sprintf("%.2f",p.σ)
    print(io,"$(typeof(p))($(lat)˚N, $(lon)˚E, σ=$σ)")
end

Particle() = Particle{DEFAULT_NF}()
Base.zero(::Type{Particle}) = Particle{DEFAULT_NF}()
Base.zero(::Type{Particle{NF}}) where NF = Particle{NF}()
Base.zero(::T) where {T<:Particle} = zero(T)

export activate, deactivate
activate(p::P) where {P<:Particle} = P(p.lon,p.lat,p.σ,true)
deactivate(p::P) where {P<:Particle} = P(p.lon,p.lat,p.σ,false)
function Base.mod(p::P) where {P<:Particle}
    (;lon, lat, σ, active) = p
    crossed_pole = lat > 90 || lat < -90
    lat = lat - 2crossed_pole*((lat + 90) % 180)
    lon = (lon + 180*crossed_pole) % 360
    σ = clamp(σ,0,1)
    return P(;lat,lon,σ,active)
end

abstract type AbstractParticleAdvection <: AbstractModelComponent end

export ParticleAdvection
Base.@kwdef struct ParticleAdvection{
    NF<:AbstractFloat,
    Grid<:AbstractGrid,
} <: AbstractParticleAdvection

    "Number of latitude rings on one hemisphere, Equator incl. Resolution parameter."
    nlat_half::Int

    "Number of particles to advect"
    n_particles::Int = 10

    "Work array"
    a::Vector{NF} = zeros(n_particles)
    b::Vector{NF} = zeros(n_particles)

    "Interpolator to interpolate velocity fields onto particle positions"
    interpolator::AnvilInterpolator{NF,Grid} = AnvilInterpolator(NF, Grid, nlat_half, n_particles)
end

ParticleAdvection(SG::SpectralGrid;kwargs...) = ParticleAdvection{SG.NF,SG.Grid}(;nlat_half=SG.nlat_half,kwargs...)

function initialize!(
    particles::Vector{P},
    model::ModelSetup,
) where {P<:Particle}

    for i in eachindex(particles)
        particles[i] = P(lat=5*randn(),lon=5*randn())
    end
end

function particle_advection!(
    particles::Vector{Particle{NF}},
    diagn,
    dt::Real,
    particle_advection::AbstractParticleAdvection,
) where NF

    (;u_grid, v_grid) = diagn.layers[1].grid_variables

    lats = particle_advection.a
    lons = particle_advection.b

    for i in eachindex(particles,lats,lons)
        lats[i] = particles[i].lat
        lons[i] = particles[i].lon
    end

    (;interpolator) = particle_advection
    RingGrids.update_locator!(interpolator,lats,lons)

    degrees_per_meter = convert(NF,360/2π)
    dt_NF = convert(NF,dt)
    u_interp = particle_advection.a
    v_interp = particle_advection.b

    interpolate!(u_interp,u_grid,interpolator)
    interpolate!(v_interp,v_grid,interpolator)

    for i in eachindex(particles,u_interp,v_interp)
        (;lon,lat) = particles[i]
        lon += u_interp[i] * dt_NF * degrees_per_meter/cosd(lat)
        lat += v_interp[i] * dt_NF * degrees_per_meter
        particles[i] = mod(Particle(lat,lon,zero(lat),true))
    end

    return nothing
end