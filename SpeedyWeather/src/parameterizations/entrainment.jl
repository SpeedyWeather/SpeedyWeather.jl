abstract type AbstractEntrainment <: AbstractParameterization end

export NoEntrainment, LinearEntrainment, ConstantEntrainment

"""No entrainment: the rising parcel in the Betts-Miller convection scheme does not mix with
its environment, reproducing the scheme without entrainment. The default for
[`BettsMillerDryConvection`](@ref); [`BettsMillerConvection`](@ref) defaults to
[`LinearEntrainment`](@ref) instead."""
struct NoEntrainment <: AbstractEntrainment end
NoEntrainment(::SpectralGrid; kwargs...) = NoEntrainment()
Adapt.@adapt_structure NoEntrainment

"""$(TYPEDSIGNATURES) No entrainment at any σ level."""
@inline (::NoEntrainment)(σ) = zero(σ)
@inline (::NoEntrainment)(σ, Δσ) = zero(σ)

"""$(TYPEDSIGNATURES)
Fraction of environmental air mixed into the rising parcel over a σ distance `Δσ` at
level `σ`. The entrainment profile `E(σ)` is the fraction mixed over the reference distance
`E.reference_layer_thickness` (Δσ_ref), applied as a rate per unit σ so that the dilution through a
column is independent of the number of layers: the fraction of the parcel that remains is
`(1 - E(σ))^(Δσ/Δσ_ref)` which multiplies to the same over many thin or few thick layers, see #1426."""
@inline function (E::AbstractEntrainment)(σ, Δσ)
    remaining = (1 - min(E(σ), one(σ)))^(Δσ / E.reference_layer_thickness)
    return 1 - remaining
end

"""Whether an entrainment profile actually mixes environmental air into the rising parcel.
`false` for [`NoEntrainment`](@ref) so that the mixing step in the ascent can be skipped
(compiled away) entirely for the default case; `true` for all other entrainment profiles."""
entraining(::AbstractEntrainment) = true
entraining(::NoEntrainment) = false

"""Linear entrainment profile: entrainment rate decreases linearly from `surface_entrainment`
at the surface (σ=1) to zero at `σ_entrainment`. Above (σ < σ_entrainment), entrainment is
zero. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct LinearEntrainment{NF} <: AbstractEntrainment
    "[OPTION] Sigma level at which entrainment becomes zero [1]"
    @param σ_entrainment::NF = 0.5 (bounds = 0 .. 0.99,)

    "[OPTION] Entrainment rate at the surface [1] per reference_layer_thickness"
    @param surface_entrainment::NF = 0.5 (bounds = 0 .. 1,)

    "[OPTION] Reference σ distance [1] over which the entrainment rate is the mixed fraction"
    @param reference_layer_thickness::NF = 1 / 8 (bounds = Positive,)
end

Adapt.@adapt_structure LinearEntrainment
LinearEntrainment(SG::SpectralGrid; kwargs...) = LinearEntrainment{SG.NF}(; kwargs...)
LinearEntrainment(::Type{NF}; kwargs...) where {NF} = LinearEntrainment{NF}(; kwargs...)

"""$(TYPEDSIGNATURES) Entrainment rate at sigma level σ for the linear profile."""
@inline (E::LinearEntrainment)(σ) = E.surface_entrainment * max(zero(σ), σ - E.σ_entrainment) / (1 - E.σ_entrainment)

"""Constant entrainment profile: entrainment rate is the same at every model level.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantEntrainment{NF} <: AbstractEntrainment
    "[OPTION] Constant entrainment rate [1] per reference_layer_thickness"
    @param entrainment_rate::NF = 0.2 (bounds = 0 .. 1,)

    "[OPTION] Reference σ distance [1] over which the entrainment rate is the mixed fraction"
    @param reference_layer_thickness::NF = 1 / 8 (bounds = Positive,)
end

Adapt.@adapt_structure ConstantEntrainment
ConstantEntrainment(SG::SpectralGrid; kwargs...) = ConstantEntrainment{SG.NF}(; kwargs...)
ConstantEntrainment(::Type{NF}; kwargs...) where {NF} = ConstantEntrainment{NF}(; kwargs...)

"""$(TYPEDSIGNATURES) Entrainment rate at sigma level σ for the constant profile."""
@inline (E::ConstantEntrainment)(σ) = E.entrainment_rate
