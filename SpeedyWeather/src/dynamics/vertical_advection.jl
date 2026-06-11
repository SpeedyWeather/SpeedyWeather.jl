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

    Δσ = model.geometry.σ_levels_thick
    advection_scheme = model.vertical_advection
    (; w) = vars.dynamics

    # unrolled over compile-time variable names (instead of a loop over runtime symbols)
    # to avoid Union-typed variables which Enzyme cannot differentiate
    vertical_advection!(Val(:u), vars, w, Δσ, advection_scheme, model)
    vertical_advection!(Val(:v), vars, w, Δσ, advection_scheme, model)
    vertical_advection!(Val(:temperature), vars, w, Δσ, advection_scheme, model)
    vertical_advection!(Val(:humidity), vars, w, Δσ, advection_scheme, model)

    for (name, tracer) in model.tracers
        if tracer.active
            ξ_tend = get_tendency_step(vars.tendencies.grid_tracers[name], model.time_stepping, advection_scheme)
            ξ = get_prognostic_step(vars.grid.tracers[name], model.time_stepping, advection_scheme)
            _vertical_advection!(ξ_tend, w, ξ, Δσ, advection_scheme)
        end
    end
    return nothing
end

# var is a compile-time constant so that haskey and getproperty constant-fold to
# concrete variables (type-stable, required for Enzyme differentiability)
@inline function vertical_advection!(::Val{var}, vars::Variables, w, Δσ, advection_scheme, model) where {var}
    haskey(vars.tendencies.grid, var) || return nothing
    ξ_tend = get_tendency_step(vars.tendencies.grid[var], model.time_stepping, advection_scheme)
    ξ = get_prognostic_step(vars.grid[var], model.time_stepping, advection_scheme)
    return _vertical_advection!(ξ_tend, w, ξ, Δσ, advection_scheme)
end

function _vertical_advection!(
        ξ_tend::AbstractField,      # tendency of quantity ξ
        w::AbstractField,           # vertical velocity at k+1/2
        ξ::AbstractField,           # ξ
        Δσ,                         # layer thickness on σ levels
        adv::VerticalAdvection      # vertical advection scheme of order B
    )
    grids_match(ξ_tend, w, ξ) || throw(DimensionMismatch(ξ_tend, w, ξ))

    nlayers = size(ξ, 2)
    arch = architecture(ξ_tend)

    launch!(
        arch, RingGridWorkOrder, size(ξ_tend),
        vertical_advection_kernel!,
        ξ_tend, w, ξ, Δσ, nlayers, adv
    )
    return nothing
end

@kernel inbounds = true function vertical_advection_kernel!(
        ξ_tend, w, ξ, Δσ, nlayers, adv
    )
    ij, k = @index(Global, NTuple)

    Δσₖ⁻¹ = inv(Δσ[k])

    # for k=1 "above" term (at k-1/2) is 0, for k==nlayers "below" term (at k+1/2) is zero
    k⁻ = max(1, k - 1)
    k⁺ = k

    k_stencil = retrieve_stencil(k, nlayers, adv)

    w⁻ = w[ij, k⁻]
    w⁺ = w[ij, k⁺]

    # tail/front instead of [2:end]/[1:end-1] as tuple-range indexing is not type-stable
    ξᶠ⁺ = reconstructed_at_face(ξ, ij, Base.tail(k_stencil), w⁺, adv)
    ξᶠ⁻ = reconstructed_at_face(ξ, ij, Base.front(k_stencil), w⁻, adv)

    # -= as the tendencies already contain the parameterizations
    ξ_tend[ij, k] -= Δσₖ⁻¹ * (w⁺ * ξᶠ⁺ - w⁻ * ξᶠ⁻ - ξ[ij, k] * (w⁺ - w⁻))
end

# 1st order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 1}) where {NF} =
    ifelse(
    u > 0, ξ[ij, k[1]],
    ξ[ij, k[2]]
)

# 3rd order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 2}) where {NF} =
    ifelse(
    u > 0, (2ξ[ij, k[1]] + 5ξ[ij, k[2]] - ξ[ij, k[3]]) * 1 // 6,
    (2ξ[ij, k[4]] + 5ξ[ij, k[3]] - ξ[ij, k[2]]) * 1 // 6
)

# 5th order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 3}) where {NF} =
    ifelse(
    u > 0, (2ξ[ij, k[1]] - 13ξ[ij, k[2]] + 47ξ[ij, k[3]] + 27ξ[ij, k[4]] - 3ξ[ij, k[5]]) * 1 // 60,
    (2ξ[ij, k[6]] - 13ξ[ij, k[5]] + 47ξ[ij, k[4]] + 27ξ[ij, k[3]] - 3ξ[ij, k[2]]) * 1 // 60
)

# 2nd order centered
@inline reconstructed_at_face(ξ, ij, k, u, ::CenteredVerticalAdvection{NF, 1}) where {NF} =
    (ξ[ij, k[1]] + ξ[ij, k[2]]) * 1 // 2

# 4th order centered
@inline reconstructed_at_face(ξ, ij, k, u, ::CenteredVerticalAdvection{NF, 2}) where {NF} =
    (-ξ[ij, k[1]] + 7ξ[ij, k[2]] + 7ξ[ij, k[3]] - ξ[ij, k[4]]) * 1 // 12

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

@inline function reconstructed_at_face(ξ, ij, k, u, ::WENOVerticalAdvection)
    if u > 0
        S₀ = (ξ[ij, k[3]], ξ[ij, k[4]], ξ[ij, k[5]])
        S₁ = (ξ[ij, k[2]], ξ[ij, k[3]], ξ[ij, k[4]])
        S₂ = (ξ[ij, k[1]], ξ[ij, k[2]], ξ[ij, k[3]])
    else
        S₀ = (ξ[ij, k[4]], ξ[ij, k[3]], ξ[ij, k[2]])
        S₁ = (ξ[ij, k[5]], ξ[ij, k[4]], ξ[ij, k[3]])
        S₂ = (ξ[ij, k[6]], ξ[ij, k[5]], ξ[ij, k[4]])
    end
    return weno_reconstruction(S₀, S₁, S₂)
end
