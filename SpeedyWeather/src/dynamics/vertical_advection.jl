# LAYER THICKNESS WEIGHTS FOR THE VERTICAL ADVECTION STENCIL ----------------------------
#
# The vertical advection kernel needs 1/δ_k (Δp_k/pₛ, pressure_thickness_ratio) instead of
# the constant 1/Δσ_k for hybrid sigma-pressure coordinates, but δ_k is spatially varying
# (needs `ij` and `surface_pressure`) while Δσ_k is a plain per-layer constant. These two
# tiny wrapper types let `vertical_advection_kernel!` stay a single method for both
# coordinates: `layer_thickness(thickness, ij, k)` dispatches to a constant-array lookup for
# `SigmaCoordinates` (compiling to exactly what it did before) or to
# `pressure_thickness_ratio` for `SigmaPressureCoordinates`. Both wrappers are `isbits` and
# `Adapt`-able so they pass through `launch!`/kernel argument adaption unchanged.
struct SigmaLayerThickness{VectorType}
    "Nominal σ layer thickness Δσ_k, one value per layer"
    Δσ::VectorType
end

Adapt.@adapt_structure SigmaLayerThickness

struct HybridLayerThickness{C, F}
    "Hybrid sigma-pressure vertical coordinate"
    coordinate::C
    "Surface pressure [Pa] at the dynamical core time step, one value per grid point"
    surface_pressure::F
end

Adapt.@adapt_structure HybridLayerThickness

"""$(TYPEDSIGNATURES)
Layer thickness weight δ_k = Δp_k/pₛ at grid point `ij`, layer `k`. Equals the constant
Δσ_k for `SigmaLayerThickness` (sigma coordinates); is pₛ-dependent for
`HybridLayerThickness` (hybrid sigma-pressure coordinates)."""
@inline layer_thickness(t::SigmaLayerThickness, ij, k) = t.Δσ[k]
@inline layer_thickness(t::HybridLayerThickness, ij, k) =
    pressure_thickness_ratio(k, t.surface_pressure[ij], t.coordinate)

# construct the appropriate wrapper from the model's vertical coordinate
layer_thickness_weight(::SigmaCoordinates, geometry::Geometry, surface_pressure) =
    SigmaLayerThickness(geometry.σ_levels_thick)
layer_thickness_weight(coordinate::SigmaPressureCoordinates, geometry::Geometry, surface_pressure) =
    HybridLayerThickness(coordinate, surface_pressure)

abstract type AbstractVerticalAdvection <: AbstractModelComponent end
abstract type VerticalAdvection{NF, B} <: AbstractVerticalAdvection end

# Dispersive and diffusive advection schemes `NF` is the type, `B` the half-stencil size
abstract type DiffusiveVerticalAdvection{NF, B} <: VerticalAdvection{NF, B} end
abstract type DispersiveVerticalAdvection{NF, B} <: VerticalAdvection{NF, B} end

export UpwindVerticalAdvection, WENOVerticalAdvection, CenteredVerticalAdvection
struct UpwindVerticalAdvection{NF, B} <: DiffusiveVerticalAdvection{NF, B} end
struct WENOVerticalAdvection{NF} <: DiffusiveVerticalAdvection{NF, 3} end
struct CenteredVerticalAdvection{NF, B} <: DispersiveVerticalAdvection{NF, B} end

CenteredVerticalAdvection(spectral_grid; order = 2) = CenteredVerticalAdvection{spectral_grid.NF, order ÷ 2}()
UpwindVerticalAdvection(spectral_grid; order = 5) = UpwindVerticalAdvection{spectral_grid.NF, (order + 1) ÷ 2}()
WENOVerticalAdvection(spectral_grid) = WENOVerticalAdvection{spectral_grid.NF}()

@inline function retrieve_stencil(k, nlayers, ::VerticalAdvection{NF, B}) where {NF, B}
    # creates allocation-free tuples for k-B:k+B but clamped into (1, nlayers)
    # e.g. (1, 1, 2), (1, 2, 3), (2, 3, 4) ... (for k=1, 2, 3; B=1)
    return ntuple(i -> clamp(i + k - B - 1, 1, nlayers), 2B + 1)
end

function vertical_advection!(vars::Variables, model)

    # δ_k(pₛ) = Δp_k/pₛ for hybrid sigma-pressure coordinates, the constant Δσ_k for sigma
    # coordinates (see the SigmaLayerThickness/HybridLayerThickness wrapper types above).
    thickness = layer_thickness_weight(model.geometry.vertical_coordinates, model.geometry, vars.dynamics.surface_pressure)
    advection_scheme = model.vertical_advection
    (; w) = vars.dynamics

    # unrolled over compile-time variable names (instead of a loop over runtime symbols)
    # to avoid Union-typed variables which Enzyme cannot differentiate
    vertical_advection!(Val(:u), vars, w, thickness, advection_scheme, model)
    vertical_advection!(Val(:v), vars, w, thickness, advection_scheme, model)
    vertical_advection!(Val(:temperature), vars, w, thickness, advection_scheme, model)
    vertical_advection!(Val(:humidity), vars, w, thickness, advection_scheme, model)

    for (name, tracer) in model.tracers
        if tracer.active
            ξ_tend = vars.tendencies.grid_tracers[name]
            ξ = vars.grid.tracers[name]
            s_tend = which_tendency_step(ξ_tend, model.time_stepping, advection_scheme)
            s_prog = which_prognostic_step(ξ, model.time_stepping, advection_scheme, model)
            _vertical_advection!(ξ_tend, s_tend, w, ξ, s_prog, thickness, advection_scheme)
        end
    end
    return nothing
end

# var is a compile-time constant so that haskey and getproperty constant-fold to
# concrete variables (type-stable, required for Enzyme differentiability)
@inline function vertical_advection!(::Val{var}, vars::Variables, w, thickness, advection_scheme, model) where {var}
    haskey(vars.tendencies.grid, var) || return nothing
    # Pass the full step-dimensioned fields plus the step index, rather than a `get_*_step`
    # view. The stencil kernel reads ξ at many vertical offsets per point, and indexing a
    # (`SubArray`-backed) view there is ~2x slower than indexing the contiguous parent array
    # (the step folds into the index).
    ξ_tend = vars.tendencies.grid[var]
    ξ = vars.grid[var]
    s_tend = which_tendency_step(ξ_tend, model.time_stepping, advection_scheme)
    s_prog = which_prognostic_step(ξ, model.time_stepping, advection_scheme)
    return _vertical_advection!(ξ_tend, s_tend, w, ξ, s_prog, thickness, advection_scheme)
end

function _vertical_advection!(
        ξ_tend::AbstractField,      # tendency of quantity ξ (full, with step dimension)
        s_tend::Integer,            # step of ξ_tend to write into
        w::AbstractField,           # vertical velocity at k+1/2
        ξ::AbstractField,           # ξ (full, with step dimension)
        s_prog::Integer,            # step of ξ to advect
        thickness,                  # SigmaLayerThickness or HybridLayerThickness, layer thickness weight δ_k(pₛ)
        adv::VerticalAdvection      # vertical advection scheme of order B
    )
    grids_match(ξ_tend, w, ξ) || throw(DimensionMismatch(ξ_tend, w, ξ))

    nlayers = size(ξ, 2)
    arch = architecture(ξ_tend)

    # worksize is the horizontal × vertical iteration space (skip the step dimension)
    launch!(
        arch, RingGridWorkOrder, (size(ξ_tend, 1), size(ξ_tend, 2)),
        vertical_advection_kernel!,
        ξ_tend, s_tend, w, ξ, s_prog, thickness, nlayers, adv
    )
    return nothing
end

@kernel inbounds = true function vertical_advection_kernel!(
        ξ_tend, s_tend, w, ξ, s_prog, thickness, nlayers, adv
    )
    ij, k = @index(Global, NTuple)

    Δσₖ⁻¹ = inv(layer_thickness(thickness, ij, k))

    # for k=1 "above" term (at k-1/2) is 0, for k==nlayers "below" term (at k+1/2) is zero
    k⁻ = max(1, k - 1)
    k⁺ = k

    k_stencil = retrieve_stencil(k, nlayers, adv)

    w⁻ = w[ij, k⁻]
    w⁺ = w[ij, k⁺]

    # `s_prog` selects which time step of the step-dimensioned ξ to advect; indexing the
    # full array as ξ[ij, k, s_prog] keeps the contiguous parent (the constant step folds
    # into the index), which is ~2x faster here than a `get_*_step` SubArray view would be.
    # tail/front instead of [2:end]/[1:end-1] as tuple-range indexing is not type-stable
    ξᶠ⁺ = reconstructed_at_face(ξ, ij, s_prog, Base.tail(k_stencil), w⁺, adv)
    ξᶠ⁻ = reconstructed_at_face(ξ, ij, s_prog, Base.front(k_stencil), w⁻, adv)

    # -= as the tendencies already contain the parameterizations
    ξ_tend[ij, k, s_tend] -= Δσₖ⁻¹ * (w⁺ * ξᶠ⁺ - w⁻ * ξᶠ⁻ - ξ[ij, k, s_prog] * (w⁺ - w⁻))
end

# reconstructed_at_face indexes ξ[ij, k[i], s]: `k` is the vertical stencil (tuple of
# layer indices) and `s` selects the time step of the s-dimensioned field.

# 1st order upwind
@inline reconstructed_at_face(ξ, ij, s, k, u, ::UpwindVerticalAdvection{NF, 1}) where {NF} =
    ifelse(
    u > 0, ξ[ij, k[1], s],
    ξ[ij, k[2], s]
)

# 3rd order upwind
@inline reconstructed_at_face(ξ, ij, s, k, u, ::UpwindVerticalAdvection{NF, 2}) where {NF} =
    ifelse(
    u > 0, (2ξ[ij, k[1], s] + 5ξ[ij, k[2], s] - ξ[ij, k[3], s]) * 1 // 6,
    (2ξ[ij, k[4], s] + 5ξ[ij, k[3], s] - ξ[ij, k[2], s]) * 1 // 6
)

# 5th order upwind
@inline reconstructed_at_face(ξ, ij, s, k, u, ::UpwindVerticalAdvection{NF, 3}) where {NF} =
    ifelse(
    u > 0, (2ξ[ij, k[1], s] - 13ξ[ij, k[2], s] + 47ξ[ij, k[3], s] + 27ξ[ij, k[4], s] - 3ξ[ij, k[5], s]) * 1 // 60,
    (2ξ[ij, k[6], s] - 13ξ[ij, k[5], s] + 47ξ[ij, k[4], s] + 27ξ[ij, k[3], s] - 3ξ[ij, k[2], s]) * 1 // 60
)

# 2nd order centered
@inline reconstructed_at_face(ξ, ij, s, k, u, ::CenteredVerticalAdvection{NF, 1}) where {NF} =
    (ξ[ij, k[1], s] + ξ[ij, k[2], s]) * 1 // 2

# 4th order centered
@inline reconstructed_at_face(ξ, ij, s, k, u, ::CenteredVerticalAdvection{NF, 2}) where {NF} =
    (-ξ[ij, k[1], s] + 7ξ[ij, k[2], s] + 7ξ[ij, k[3], s] - ξ[ij, k[4], s]) * 1 // 12

const ε = 1 // 1_000_000    # = 1e-6 but number format flexible
const d₀ = 3 // 10
const d₁ = 3 // 5
const d₂ = 1 // 10

@inline weight_β₀(S) = 13 // 12 * (S[1] - 2S[2] + S[3])^2 + 1 // 4 * (3S[1] - 4S[2] + S[3])^2
@inline weight_β₁(S) = 13 // 12 * (S[1] - 2S[2] + S[3])^2 + 1 // 4 * (S[1] - S[3])^2
@inline weight_β₂(S) = 13 // 12 * (S[1] - 2S[2] + S[3])^2 + 1 // 4 * (S[1] - 4S[2] + 3S[3])^2

@inline p₀(S) = (2S[1] + 5S[2] - S[3]) * 1 // 6     # downind stencil
@inline p₁(S) = (-S[1] + 5S[2] + 2S[3]) * 1 // 6    # upwind stencil
@inline p₂(S) = (2S[1] - 7S[2] + 11S[3]) * 1 // 6   # extrapolating stencil

@inline τ₅(β₀, β₁, β₂) = abs(β₂ - β₀)

@inline function weno_reconstruction(S₀, S₁, S₂)
    β₀ = weight_β₀(S₀)
    β₁ = weight_β₁(S₁)
    β₂ = weight_β₂(S₂)

    w₀ = d₀ * (1 + (τ₅(β₀, β₁, β₂) / (β₀ + ε))^2)
    w₁ = d₁ * (1 + (τ₅(β₀, β₁, β₂) / (β₁ + ε))^2)
    w₂ = d₂ * (1 + (τ₅(β₀, β₁, β₂) / (β₂ + ε))^2)

    w₀, w₁, w₂ = (w₀, w₁, w₂) ./ (w₀ + w₁ + w₂)

    return p₀(S₀) * w₀ + p₁(S₁) * w₁ + p₂(S₂) * w₂
end

@inline function reconstructed_at_face(ξ, ij, s, k, u, ::WENOVerticalAdvection)
    if u > 0
        S₀ = (ξ[ij, k[3], s], ξ[ij, k[4], s], ξ[ij, k[5], s])
        S₁ = (ξ[ij, k[2], s], ξ[ij, k[3], s], ξ[ij, k[4], s])
        S₂ = (ξ[ij, k[1], s], ξ[ij, k[2], s], ξ[ij, k[3], s])
    else
        S₀ = (ξ[ij, k[4], s], ξ[ij, k[3], s], ξ[ij, k[2], s])
        S₁ = (ξ[ij, k[5], s], ξ[ij, k[4], s], ξ[ij, k[3], s])
        S₂ = (ξ[ij, k[6], s], ξ[ij, k[5], s], ξ[ij, k[4], s])
    end
    return weno_reconstruction(S₀, S₁, S₂)
end
