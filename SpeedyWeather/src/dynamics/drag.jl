abstract type AbstractDrag <: AbstractModelComponent end

# function barrier for all drags to unpack model.drag
function drag!(
        vars::Variables,
        lf::Integer,
        model::AbstractModel,
    )
    return drag!(vars, model.drag, lf, model)
end

# NO DRAG
drag!(vars, drag::Nothing, args...) = nothing

export LinearDrag

"""Linear boundary layer drag following Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
@parameterized @kwdef struct LinearDrag{NF, VectorType} <: AbstractDrag
    "[OPTION] Sigma coordinate below which linear drag is applied"
    @param σb::NF = 0.7 (bounds = 0 .. 1,)

    "[OPTION] Time scale for linear drag coefficient at σ=1 (=1/kf in HS96)"
    time_scale::Second = Day(1)

    "[DERIVED] Precomputed drag coefficients for each layer"
    drag_coefs::VectorType
end

"""
$(TYPEDSIGNATURES)
Generator function using `nlayers` from `SG::SpectralGrid`"""
LinearDrag(SG::SpectralGrid; kwargs...) = LinearDrag{SG.NF, SG.VectorType}(drag_coefs = zeros(SG.NF, SG.nlayers); kwargs...)

"""
$(TYPEDSIGNATURES)
Precomputes the drag coefficients for the `LinearDrag` scheme."""
function initialize!(drag::LinearDrag, model::PrimitiveEquation)

    (; σ_levels_full) = model.geometry
    (; σb, time_scale, drag_coefs) = drag
    kf = 1 / time_scale.value

    # drag only below σb, lin increasing to kf at σ=1
    @. drag_coefs = kf * max(0, (σ_levels_full - σb) / (1 - σb))
    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute tendency for boundary layer drag of a `column` and add to its tendencies fields"""
function drag!(
        vars::Variables,
        drag::LinearDrag,
        lf::Integer,
        model::AbstractModel,
    )
    (; u, v) = vars.grid
    Fu = vars.tendencies.grid.u
    Fv = vars.tendencies.grid.v

    # include radius scaling
    c = vars.prognostic.scale[]
    @. Fu -= c * drag.drag_coefs' .* u
    @. Fv -= c * drag.drag_coefs' .* v

    return nothing
end

# Quadratic drag
export QuadraticDrag
@parameterized @kwdef mutable struct QuadraticDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1]"
    @param c_D::NF = 1.0e-12 (bounds = Nonnegative,)    # TODO is this a good default?
end

QuadraticDrag(SG::SpectralGrid; kwargs...) = QuadraticDrag{SG.NF}(; kwargs...)

initialize!(::QuadraticDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u, v)|*(u, v)

with `c_D` the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`.
`c_D` and layer thickness `H` are precomputed in `initialize!(::QuadraticDrag, ::AbstractModel)`
and scaled by the radius as are the momentum equations."""
function drag!(
        vars::Variables,
        drag::QuadraticDrag,
        lf::Integer,
        model::AbstractModel,
    )
    k = size(vars.grid.u, 2)            # drag only on surface layer
    u = field_view(vars.grid.u, :, k)
    v = field_view(vars.grid.v, :, k)

    Fu = field_view(vars.tendencies.grid.u, :, k)
    Fv = field_view(vars.tendencies.grid.v, :, k)

    # total drag coefficient with radius scaling
    c = drag.c_D / model.atmosphere.layer_thickness
    c *= vars.prognostic.scale[]^2

    return launch!(
        architecture(Fu), LinearWorkOrder, size(Fu), quadratic_drag_kernel!,
        Fu, Fv, u, v, c
    )
end

@kernel inbounds = true function quadratic_drag_kernel!(
        Fu, Fv, u, v, @Const(c)
    )
    ij = @index(Global, Linear)

    # Calculate speed at surface layer k
    speed = sqrt(u[ij]^2 + v[ij]^2)

    # Apply quadratic drag, -= as the tendencies already contain forcing
    Fu[ij] -= c * speed * u[ij]
    Fv[ij] -= c * speed * v[ij]
end

export LinearVorticityDrag
@parameterized @kwdef mutable struct LinearVorticityDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1/s]"
    @param c::NF = 1.0e-7 (bounds = Nonnegative,)
end

LinearVorticityDrag(SG::SpectralGrid; kwargs...) = LinearVorticityDrag{SG.NF}(; kwargs...)

initialize!(::LinearVorticityDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Linear drag for the vorticity equations of the form F = -cξ
with c drag coefficient [1/s]."""
function drag!(
        vars::Variables,
        drag::LinearVorticityDrag,
        lf::Integer,
        model::AbstractModel,
    )
    vor_tend = vars.tendencies.vor
    vor = get_step(vars.prognostic.vor, lf)

    # scale by radius (but only once, the second radius is in vor)
    c = drag.c * vars.prognostic.scale[]
    vor_tend.data .-= c * vor.data      # use .data to bypass conflicting broadcasting

    return nothing
end

export JetDrag
@parameterized @kwdef struct JetDrag{NF, SpectralVariable2D} <: SpeedyWeather.AbstractDrag
    "[OPTION] Relaxation time scale τ"
    time_scale::Second = Day(6)

    "[OPTION] Jet strength [m/s]"
    @param u₀::NF = 20

    "[OPTION] latitude of Gaussian jet [˚N]"
    @param latitude::NF = 30 (bounds = -90 .. 90,)

    "[OPTION] Width of Gaussian jet [˚]"
    @param width::NF = 6 (bounds = Positive,)

    # TO BE INITIALISED
    "Relaxation back to reference vorticity"
    ζ₀::SpectralVariable2D
end

function JetDrag(SG::SpectralGrid; kwargs...)
    ζ₀ = zeros(Complex{SG.NF}, SG.spectrum)
    return JetDrag{SG.NF, SG.SpectralVariable2D}(; ζ₀, kwargs...)
end

function initialize!(drag::JetDrag, model::AbstractModel)
    (; spectral_grid, geometry) = model
    (; grid, NF) = spectral_grid
    u = zeros(NF, grid)

    lat = geometry.latds

    for ij in eachindex(u)
        u[ij] = drag.u₀ * exp(-(lat[ij] - drag.latitude)^2 / (2 * drag.width^2))
    end

    û = transform(u, model.spectral_transform)
    v̂ = zero(û)
    curl!(drag.ζ₀, û, v̂, model.spectral_transform)
    return nothing
end

function drag!(
        vars::Variables,
        drag::JetDrag,
        lf::Integer,
        model::AbstractModel,
    )
    vor = get_step(vars.prognostic.vor, lf)
    vor_tend = vars.tendencies.vor
    (; ζ₀) = drag

    # scale by radius as is vorticity
    s = vars.prognostic.scale[]
    r = s / drag.time_scale.value

    k = size(vor, 2)   # drag only on surface layer

    # GPU kernel launch
    arch = architecture(vars.grid.u)
    return launch!(
        arch, LinearWorkOrder, (size(vor_tend, 1),), jet_drag_kernel!,
        vor_tend, vor, ζ₀, r, k
    )
end

@kernel inbounds = true function jet_drag_kernel!(
        vor_tend, vor, ζ₀, @Const(r), @Const(k)
    )
    lm = @index(Global, Linear)
    vor_tend[lm, k] -= r * (vor[lm, k] - ζ₀[lm])
end
