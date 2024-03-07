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

# empty constructor map to zero
Particle() = zero(Particle)
Particle{NF}() where NF = zero(Particle{NF})
Particle{NF,isactive}() where {NF,isactive} = zero(Particle{NF,isactive})
Particle{NF}(args...) where NF = Particle{NF,true}(args...)

# promotion of arguments
Particle(lon,lat) = Particle(lon,lat,0)
Particle(lon::Integer,lat::Integer) = Particle(lon,lat,0)
Particle(lon,lat,σ) = Particle{promote_type(typeof.((lon,lat,σ))...),true}(lon,lat,σ)
Particle(lon::Integer,lat::Integer,σ::Integer) = Particle{DEFAULT_NF,true}(lon,lat,σ)

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
activate(  p::Particle{NF}) where NF = Particle{NF, true}(p.lon, p.lat, p.σ)
deactivate(p::Particle{NF}) where NF = Particle{NF,false}(p.lon, p.lat, p.σ)
active(::Particle{NF,isactive}) where {NF,isactive} = isactive

@inline function Base.mod(p::P) where {P<:Particle}
    (;lon, lat, σ) = p
    crossed_pole = lat > 90 || lat < -90
    lat = lat - 2crossed_pole*((lat + 90) % 180)
    lon = mod(lon + 180*crossed_pole, 360)
    σ = clamp(σ, 0, 1)
    return P(lon, lat, σ)
end