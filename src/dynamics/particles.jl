abstract type AbstractParticle{NF} end

export Particle

"""
Particle with location lon (longitude), lat (latitude) and σ (vertical coordinate).
Longitude is assumed to be in [0,360˚E), latitude in [-90˚,90˚N] and σ in [0,1] but
not strictly enforced at creation, see `mod(::Particle)` and `ismod(::Particle)`.
A particle is either active or inactive, determined by the Boolean in it's 2nd
type parameter. By default, a particle is active, of number format DEFAULT_NF
and at 0˚N, 0˚E, σ=0.
$(TYPEDFIELDS)"""
struct Particle{
    NF<:AbstractFloat,  # number format of coordinates
    isactive,           # bool indicating whether particle is active
} <: AbstractParticle{NF}
    "longitude in [0,360˚]"
    lon::NF
    "latitude [-90˚,90˚]"
    lat::NF
    "vertical sigma coordinate [0 (top) to 1 (surface)]"
    σ::NF
end

# keyword constructors
Particle{NF}(;lon, lat, σ=0) where NF = Particle{NF}(lon, lat, σ)
Particle{NF, isactive}(; lon, lat, σ=0) where {NF,isactive} = Particle{NF,isactive}(lon, lat, σ)
Particle(;lon, lat, σ=0) = Particle(lon, lat, σ)

# parametric constructors
Particle{NF}(lon, lat) where NF = Particle{NF,true}(lon, lat, 0)
Particle{NF, isactive}(lon, lat) where {NF, isactive} = Particle{NF, isactive}(lon, lat, 0)
Particle{NF}(lon, lat, σ) where NF = Particle{NF,true}(lon, lat, σ)

# promotion of arguments if NF not provided
Particle(lon, lat) = Particle(lon, lat, 0)
Particle(lon, lat, σ) = Particle{promote_type(typeof.((lon, lat, σ))...), true}(lon, lat, σ)
Particle(lon::Integer, lat::Integer) = Particle(lon,lat,0)
Particle(lon::Integer, lat::Integer, σ::Integer) = Particle{DEFAULT_NF, true}(lon, lat, σ)

# zero generators
Base.zero(::Type{Particle}) = Particle{DEFAULT_NF,true}(0,0,0)
Base.zero(::Type{Particle{NF}}) where NF = Particle{NF,true}(0,0,0)
Base.zero(::Type{Particle{NF,isactive}}) where {NF,isactive} = Particle{NF,isactive}(0,0,0)
Base.zero(::P) where {P<:Particle} = zero(P)

Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle}) = rand(rng, Particle{DEFAULT_NF,true})
Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle{NF}}) where NF = rand(rng, Particle{NF,true})

# rand uniformly distributed over the globe with cos-distribution for poles
function Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle{NF,isactive}}) where {NF,isactive}
    lon = 360*rand(rng,NF)                  # ∈ [0,360˚E]
    lat = asind(2rand(rng,NF)-1)            # cos-distributed latitude
    σ = rand(rng,NF)
    return Particle{NF,isactive}(lon,lat,σ)
end

# equality with same location and same activity
function Base.:(==)(p1::Particle{NF1,active1},p2::Particle{NF2,active2}) where {NF1,active1,NF2,active2}
    return  (active1 == active2) &&     # both active or both inactive
            (p1.lat == p2.lat) &&       # same latitude
            (p1.σ == p2.σ) &&           # same elevation
            ((p1.lon == p2.lon) ||      # same longitude OR
            (p1.lat*p2.lat == 8100))    # both at the north/south pole, because 90˚N, 0˚E == 90˚N, 10˚E
end

function Base.isapprox(
    p1::Particle{NF1,active1},
    p2::Particle{NF2,active2};
    kwargs...,
) where {NF1,active1,NF2,active2}
    b = (active1 == active2)                    # both active or both inactive
    b &= isapprox(p1.lat, p2.lat; kwargs...)    # same latitude
    b &= isapprox(p1.σ, p2.σ; kwargs...)        # same elevation
    c =  isapprox(p1.lon, p2.lon; kwargs...)    # same longitude OR
    b &= c | isapprox(p1.lat * p2.lat, 8100; kwargs...) # both at the north/south pole
                                                        # because 90˚N, 0˚E == 90˚N, 10˚E
    return b
end

Base.convert(::Type{Particle{NF}}, p::Particle) where NF = Particle{NF, active(p)}(p.lon, p.lat, p.σ)
Base.convert(::Type{Particle{NF, isactive}}, p::Particle) where {NF, isactive} =
                    Particle{NF, isactive}(p.lon, p.lat, p.σ)

function Base.show(io::IO,p::Particle{NF,isactive}) where {NF,isactive}
    lat = @sprintf("%6.2f",p.lat)
    lon = @sprintf("%6.2f",p.lon)
    σ = @sprintf("%.2f",p.σ)
    activity = isactive ? "  active" : "inactive"
    print(io,"Particle{$NF, $activity}($(lon)˚E, $(lat)˚N, σ = $σ)")
end

export move

"""$(TYPEDSIGNATURES)
Move a particle with increments (dlon, dlat, dσ) in those respective coordinates.
Only active particles are moved."""
@inline function move(p::Particle{NF,true}, dlon, dlat, dσ) where NF
    (;lon, lat, σ) = p
    Particle{NF,true}(lon+dlon,lat+dlat,σ+dσ)
end

"""$(TYPEDSIGNATURES)
Move a particle with increments (dlon, dlat) in 2D. No movement in vertical σ.
Only active particles are moved."""
@inline function move(p::Particle{NF,true}, dlon, dlat) where NF
    (;lon, lat, σ) = p
    Particle{NF,true}(lon+dlon,lat+dlat,σ)
end

"""$(TYPEDSIGNATURES)
Inactive particles are not moved."""
@inline move(p::Particle{NF,false}, args...) where NF = p

export activate, deactivate, active

"""$(TYPEDSIGNATURES)
Activate particle. Active particles can move."""
activate(  p::Particle{NF}) where NF = Particle{NF, true}(p.lon, p.lat, p.σ)

"""$(TYPEDSIGNATURES)
Deactivate particle. Inactive particles cannot move."""
deactivate(p::Particle{NF}) where NF = Particle{NF, false}(p.lon, p.lat, p.σ)

"""$(TYPEDSIGNATURES)
Check whether particle is active."""
active(::Particle{NF, isactive}) where {NF, isactive} = isactive

"""
$(TYPEDSIGNATURES)
Modulo operator for particle locations to map them back into [0,360˚E) and [-90˚,90˚N],
in the horizontal and to clamp vertical σ coordinates into [0,1]."""
@inline function Base.mod(p::P) where {P<:Particle}
    (;lon, lat, σ) = p
    pole_crossed = isodd((abs(lat) - 2eps(lat) + 90) ÷ 180)
    lat = 90 - abs(mod(lat+90,360) - 180)   # new latitude is wrapped around poles
    lon = mod(lon + 180*pole_crossed, 360)  # mod lon into [0,360˚E] but +180 for pole crossings
    σ = clamp(σ, 0, 1)      # particle above top (σ<0) stays at top, ground (σ>1) stays on ground
    return P(lon, lat, σ)
end

"""
$(TYPEDSIGNATURES)
Check that a particle is in longitude [0,360˚E), latitude [-90˚,90˚N], and σ in [0,1]."""
function ismod(p::Particle)
    valid::Bool = true
    valid &= -90 <= p.lat <= 90     # poles included
    valid &= 0 <= p.lon < 360       # 360˚E excluded (=0˚)
    valid &= 0 <= p.σ <= 1          # top and ground included
    return valid
end