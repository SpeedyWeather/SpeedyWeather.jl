export SemiImplicit, ImplicitShallowWater, ImplicitPrimitiveEquation
export IsothermalReferenceTemperature, GlobalMeanReferenceTemperature, MaximumReferenceTemperature

abstract type AbstractReferenceTemperature end

"""$(TYPEDSIGNATURES)
Isothermal reference temperature profile for the semi-implicit time stepping, the
implicit operators are therefore only reinitialized when the time step changes.
$(TYPEDFIELDS)"""
@kwdef struct IsothermalReferenceTemperature{NF} <: AbstractReferenceTemperature
    "[OPTION] Temperature of the isothermal reference profile [K]"
    temperature::NF = 300
end

"""$(TYPEDSIGNATURES)
Global-mean temperature profile plus `offset` as reference temperature profile for the
semi-implicit time stepping, recomputed on the `schedule` of the implicit solver.
$(TYPEDFIELDS)"""
@kwdef struct GlobalMeanReferenceTemperature{NF} <: AbstractReferenceTemperature
    "[OPTION] Offset [K] added to the global-mean temperature profile"
    offset::NF = 0
end

"""$(TYPEDSIGNATURES)
Isothermal reference temperature for the semi-implicit time stepping from the current
maximum temperature of the model plus `margin`, recomputed on the `schedule` of the implicit
solver. The semi-implicit scheme is stable for a warm isothermal reference temperature, warmer
than the atmosphere (Simmons, Hoskins and Burridge, 1978), but the warmer the reference the more
the gravity waves treated implicitly are slowed down. In contrast to `IsothermalReferenceTemperature`
this follows the model's climate, e.g. for a snowball Earth or a hot climate.
$(TYPEDFIELDS)"""
@kwdef struct MaximumReferenceTemperature{NF} <: AbstractReferenceTemperature
    "[OPTION] Margin [K] added to the maximum temperature"
    margin::NF = 30
end

"""$(TYPEDSIGNATURES)
Global-mean temperature profile (broadcastable, one value per layer) from the l=m=0 harmonic
of the spectral temperature `temp` (full temperature, not the anomaly), using `norm_sphere` of
the spectral transform."""
global_mean_temperature(temp::LowerTriangularArray, norm_sphere) =
    Base.Broadcast.broadcasted(t -> real(t) / norm_sphere, view(temp.data, 1, :))

# reference temperature profile (broadcastable) from the spectral temperature `temp` (for its
# global-mean profile) or the grid temperature `temp_grid` (full temperature, not the anomaly)
reference_temperature(ref::IsothermalReferenceTemperature, temp, temp_grid, norm_sphere) = ref.temperature
reference_temperature(ref::GlobalMeanReferenceTemperature, temp, temp_grid, norm_sphere) =
    Base.Broadcast.broadcasted(+, global_mean_temperature(temp, norm_sphere), ref.offset)
reference_temperature(ref::MaximumReferenceTemperature, temp, temp_grid, norm_sphere) =
    maximum(temp_grid) + ref.margin

# whether the reference temperature profile has to be recomputed on the schedule
rescheduled(::IsothermalReferenceTemperature) = false
rescheduled(::GlobalMeanReferenceTemperature) = true
rescheduled(::MaximumReferenceTemperature) = true

"""
Semi-implicit solver for the gravity waves, holding only its options. The same solver is used
for all models, its behaviour is defined by dispatch on the model type: no correction for the
`BarotropicModel` (no gravity waves), the shallow water gravity wave terms for the
`ShallowWaterModel` and the linearised gravity wave terms about a reference temperature profile
for the primitive equation models, whose precomputed operators are variables of the model in
`vars.dynamics.implicit`. `reference_temperature` and `schedule` are only used by the primitive
equation models. `centering` is a number (one centering for every (sub)step) or, for the Lorenz
N-cycle, a tuple with a centering per substep of the cycle, see `NCycleLorenz`.
$(TYPEDFIELDS)"""
struct SemiImplicit{NF, CenteringType, ScheduleType, ReferenceTemperatureType} <: AbstractImplicit
    "[OPTION] Centering coefficient: 0.5 = Crank-Nicolson, 1 = backward Euler, or a tuple with a centering per substep of the Lorenz N-cycle"
    centering::CenteringType

    "[OPTION] Reference temperature profile to linearise about, isothermal or the global mean (plus offset), primitive equations only"
    reference_temperature::ReferenceTemperatureType

    "[OPTION] Schedule to recompute the reference temperature profile (and the implicit operators) from the current global-mean temperature, primitive equations only"
    schedule::ScheduleType
end

"""$(TYPEDSIGNATURES)
Generator with number format `NF`."""
function SemiImplicit{NF}(;
        centering = 0.5,
        reference_temperature::AbstractReferenceTemperature = IsothermalReferenceTemperature{NF}(),
        schedule::AbstractSchedule = PeriodicSchedule(every = Day(1)),
    ) where {NF}
    centering = centering isa Tuple ? NF.(centering) : NF(centering)
    return SemiImplicit{NF, typeof(centering), typeof(schedule), typeof(reference_temperature)}(centering, reference_temperature, schedule)
end

SemiImplicit(; kwargs...) = SemiImplicit{DEFAULT_NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Generator using the number format from SpectralGrid."""
SemiImplicit(spectral_grid::SpectralGrid; kwargs...) = SemiImplicit{spectral_grid.NF}(; kwargs...)

"""Previous name of `SemiImplicit` for the shallow water model, the solver is shared by all models."""
const ImplicitShallowWater = SemiImplicit

"""Previous name of `SemiImplicit` for the primitive equation models, the solver is shared by all models."""
const ImplicitPrimitiveEquation = SemiImplicit

# defaults for all models, extended for the primitive equations which need variables and (re)initialization
variables(::SemiImplicit, ::AbstractModel) = ()
initialize!(::SemiImplicit, ::AbstractModel) = nothing
reinitialize!(::SemiImplicit, ::AbstractModel, ::Variables) = nothing

# no gravity waves in the barotropic model, nothing to correct
implicit_correction!(::Variables, ::SemiImplicit, ::AbstractTimeStepper, ::Barotropic) = nothing

"""$(TYPEDSIGNATURES)
Centering of the semi-implicit solver for (sub)step `k` (0-based) of a cycle, the same on every
step if `centering` is a number, otherwise `centering[k+1]` (a centering per substep)."""
substep_centering(centering::Number, k) = centering
substep_centering(centering::Tuple, k) = select_substep(centering, k)

# branch-free indexing `t[k + 1]` of a tuple with a (possibly traced, e.g. Reactant) 0-based index k
@inline select_substep(t::Tuple, k) = sum(ntuple(i -> ifelse(k == i - 1, t[i], zero(t[i])), Val(length(t))))


"""$(TYPEDSIGNATURES)
Check the semi-implicit solver for `model`: a tuple of centerings (one per substep) requires the
Lorenz N-cycle with as many steps per cycle, and `reference_temperature` and `schedule` are only
used by the primitive equation models."""
function check_implicit(implicit::SemiImplicit, model::AbstractModel)
    check_centering(implicit.centering, model.time_stepping)
    if !(model isa PrimitiveEquation) && !(implicit.reference_temperature isa IsothermalReferenceTemperature)
        @warn "SemiImplicit: reference_temperature = $(nameof(typeof(implicit.reference_temperature))) is only used by " *
            "the primitive equation models and ignored for $(nameof(typeof(model)))."
    end
    return nothing
end

check_centering(::Number, ::AbstractTimeStepper) = nothing
check_centering(centering::Tuple, time_stepping::AbstractTimeStepper) = throw(
    ArgumentError(
        "SemiImplicit: a centering per substep (tuple of length $(length(centering))) " *
            "requires NCycleLorenz, use a number for $(nameof(typeof(time_stepping)))."
    )
)
