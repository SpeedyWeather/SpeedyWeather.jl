abstract type AbstractParticle{NF} end

export Particle

"""
Particle with location lon (longitude), lat (latitude) and σ (vertical coordinate).
Longitude is assumed to be in [0,360˚E), latitude in [-90˚,90˚N] and σ in [0,1] but
not strictly enforced at creation, see `mod(::Particle)` and `ismod(::Particle)`.
A particle is either active or inactive, determined by the `active` field.
By default, a particle is active, of number format DEFAULT_NF and at 0˚N, 0˚E, σ=0.
$(TYPEDFIELDS)"""
struct Particle{
    NF<:AbstractFloat,  # number format of coordinates
} <: AbstractParticle{NF}
    "activity"
    active::Bool
    "longitude in [0,360˚]"
    lon::NF
    "latitude [-90˚,90˚]"
    lat::NF
    "vertical sigma coordinate [0 (top) to 1 (surface)]"
    σ::NF
end

# keyword constructors
Particle{NF}(;lon, lat, σ=0) where NF = Particle{NF}(true, lon, lat, σ)
Particle(;lon, lat, σ=0) = Particle(lon, lat, σ)

# parametric constructors
Particle{NF}(lon, lat) where NF = Particle{NF}(true, lon, lat, 0)
Particle{NF}(lon, lat, σ) where NF = Particle{NF}(true, lon, lat, σ)

# promotion of arguments if NF not provided
Particle(lon, lat) = Particle(lon, lat, 0)
Particle(lon, lat, σ) = Particle{promote_type(typeof.((lon, lat, σ))...)}(true, lon, lat, σ)
Particle(lon::Integer, lat::Integer) = Particle(lon,lat,0)
Particle(lon::Integer, lat::Integer, σ::Integer) = Particle{DEFAULT_NF}(true, lon, lat, σ)

# zero generators
Base.zero(::Type{Particle}) = Particle{DEFAULT_NF}(true,0,0,0)
Base.zero(::Type{Particle{NF}}) where NF = Particle{NF}(true,0,0,0)
Base.zero(::P) where {P<:Particle} = zero(P)
function Base.zeros(ArrayType::Type{<:AbstractArray{P}}, n::Int...) where {P<:Particle}
    z = ArrayType(undef, n...)
    fill!(z, zero(P))
end

Base.eltype(::Type{Particle{NF}}) where NF = NF
Base.eltype(::Particle{NF}) where NF = NF

Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle}) = rand(rng, Particle{DEFAULT_NF})

# rand uniformly distributed over the globe with cos-distribution for poles
function Base.rand(rng::Random.AbstractRNG, ::Random.Sampler{Particle{NF}}) where {NF}
    lon = 360*rand(rng,NF)                  # ∈ [0,360˚E]
    lat = asind(2rand(rng,NF)-1)            # cos-distributed latitude
    σ = rand(rng,NF)
    return Particle{NF}(lon,lat,σ)
end

# equality with same location and same activity
function Base.:(==)(p1::Particle{NF1},p2::Particle{NF2}) where {NF1,NF2}
    return  (p1.active == p2.active) && # same activity
            (p1.lat == p2.lat) &&       # same latitude
            (p1.σ == p2.σ) &&           # same elevation
            ((p1.lon == p2.lon) ||      # same longitude OR
            (p1.lat*p2.lat == 8100))    # both at the north/south pole, because 90˚N, 0˚E == 90˚N, 10˚E
end

function Base.isapprox(
    p1::Particle{NF1},
    p2::Particle{NF2};
    kwargs...,
) where {NF1,NF2}
    b = p1.active == p2.active                  # same activity
    b &= isapprox(p1.lat, p2.lat; kwargs...)    # same latitude
    b &= isapprox(p1.σ, p2.σ; kwargs...)        # same elevation
    c =  isapprox(mod(p1.lon, 360), mod(p2.lon, 360); kwargs...)    # same longitude OR
    b &= c | isapprox(p1.lat * p2.lat, 8100; kwargs...) # both at the north/south pole
                                                        # because 90˚N, 0˚E == 90˚N, 10˚E
    return b
end

Base.convert(::Type{Particle{NF}}, p::Particle) where NF = Particle{NF}(p.lon, p.lat, p.σ)

function Base.show(io::IO,p::Particle{NF}) where {NF}
    lat = @sprintf("%6.2f",p.lat)
    lon = @sprintf("%6.2f",p.lon)
    σ = @sprintf("%.2f",p.σ)
    print(io,"Particle{$NF}($(lon)˚E, $(lat)˚N, σ = $σ, $(p.active) ? active : inactive)")
end

export move

"""$(TYPEDSIGNATURES)
Move a particle with increments (dlon, dlat, dσ) in those respective coordinates."""
@inline function move(p::Particle{NF}, dlon, dlat, dσ) where NF
    if isactive(p)
        (;lon, lat, σ) = p
        return Particle{NF}(p.active, lon+dlon, lat+dlat, σ+dσ)
    else 
        return p
    end
end

"""$(TYPEDSIGNATURES)
Move a particle with increments (dlon, dlat) in 2D. No movement in vertical σ."""
@inline function move(p::Particle{NF}, dlon, dlat) where NF
    if isactive(p)
        (;lon, lat, σ) = p
        return Particle{NF}(p.active, lon+dlon, lat+dlat, σ)
    else 
        return p
    end
end

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
    return P(p.active, lon, lat, σ)
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

export activate, deactivate, isactive

"""$(TYPEDSIGNATURES)
Activate particle. Active particles can move."""
activate(  p::Particle{NF}) where NF = Particle{NF}(true, p.lon, p.lat, p.σ)

"""$(TYPEDSIGNATURES)
Deactivate particle. Inactive particles cannot move."""
deactivate(p::Particle{NF}) where NF = Particle{NF}(false, p.lon, p.lat, p.σ)

"""$(TYPEDSIGNATURES)
Check whether particle is active."""
isactive(p::Particle{NF}) where NF = p.active
