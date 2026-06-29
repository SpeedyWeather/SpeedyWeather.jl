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

    # scratch buffer (npoints) for the CPU flux-form loop below; reused across u, v,
    # temperature, humidity and tracers as it is write-before-read within each call
    face = vars.scratch.grid.a_2D

    # unrolled over compile-time variable names (instead of a loop over runtime symbols)
    # to avoid Union-typed variables which Enzyme cannot differentiate
    vertical_advection!(Val(:u), vars, w, Δσ, advection_scheme, model, face)
    vertical_advection!(Val(:v), vars, w, Δσ, advection_scheme, model, face)
    vertical_advection!(Val(:temperature), vars, w, Δσ, advection_scheme, model, face)
    vertical_advection!(Val(:humidity), vars, w, Δσ, advection_scheme, model, face)

    for (name, tracer) in model.tracers
        if tracer.active
            ξ_tend = vars.tendencies.grid_tracers[name]
            ξ = vars.grid.tracers[name]
            s_tend = which_tendency_step(ξ_tend, model.time_stepping, advection_scheme)
            s_prog = which_prognostic_step(ξ, model.time_stepping, advection_scheme, model)
            _vertical_advection!(ξ_tend, s_tend, w, ξ, s_prog, Δσ, advection_scheme, face)
        end
    end
    return nothing
end

# var is a compile-time constant so that haskey and getproperty constant-fold to
# concrete variables (type-stable, required for Enzyme differentiability)
@inline function vertical_advection!(::Val{var}, vars::Variables, w, Δσ, advection_scheme, model, face) where {var}
    haskey(vars.tendencies.grid, var) || return nothing
    # Pass the full step-dimensioned fields plus the step index, rather than a `get_*_step`
    # view. The stencil kernel reads ξ at many vertical offsets per point, and indexing a
    # (`SubArray`-backed) view there is ~2x slower than indexing the contiguous parent array
    # (the step folds into the index).
    ξ_tend = vars.tendencies.grid[var]
    ξ = vars.grid[var]
    s_tend = which_tendency_step(ξ_tend, model.time_stepping, advection_scheme)
    s_prog = which_prognostic_step(ξ, model.time_stepping, advection_scheme)
    return _vertical_advection!(ξ_tend, s_tend, w, ξ, s_prog, Δσ, advection_scheme, face)
end

function _vertical_advection!(
        ξ_tend::AbstractField,      # tendency of quantity ξ (full, with step dimension)
        s_tend::Integer,            # step of ξ_tend to write into
        w::AbstractField,           # vertical velocity at k+1/2
        ξ::AbstractField,           # ξ (full, with step dimension)
        s_prog::Integer,            # step of ξ to advect
        Δσ,                         # layer thickness on σ levels
        adv::VerticalAdvection,     # vertical advection scheme of order B
        face::AbstractField,        # scratch buffer (npoints) for the CPU flux-form loop
    )
    grids_match(ξ_tend, w, ξ) || throw(DimensionMismatch(ξ_tend, w, ξ))

    nlayers = size(ξ, 2)
    arch = architecture(ξ_tend)
    return _vertical_advection!(arch, ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv, face)
end

# launching a kernel for the per-(ij,k) iteration space is unreasonably slow on CPU (same
# issue as vertical_integration!), so the CPU path instead loops over layers k (outer) and
# points ij (inner). The "+" face of layer k (tail of its stencil) is bit-identical to the
# "-" face of layer k+1 (front of its stencil) — proved by the face-consistency test — so
# each interior face is reconstructed once and carried over in `face` instead of twice.
function _vertical_advection!(::CPU, ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv, face::AbstractField)
    npoints = size(ξ_tend, 1)

    # boundary face at k=1/2 (no neighbouring layer to carry it over from)
    k_stencil₁ = retrieve_stencil(1, nlayers, adv)
    @inbounds for ij in 1:npoints
        face[ij] = reconstruct_face(gather_stencil_values(ξ, ij, s_prog, Base.front(k_stencil₁)), w[ij, 1], adv)
    end

    for k in 1:nlayers
        k_stencil = retrieve_stencil(k, nlayers, adv)
        Δσₖ⁻¹ = inv(Δσ[k])
        @inbounds for ij in 1:npoints
            w⁺ = w[ij, k]
            w⁻ = w[ij, max(1, k - 1)]
            ξᶠ⁺ = reconstruct_face(gather_stencil_values(ξ, ij, s_prog, Base.tail(k_stencil)), w⁺, adv)
            ξᶠ⁻ = face[ij]              # = "+" face of layer k-1, carried over

            # -= as the tendencies already contain the parameterizations
            ξ_tend[ij, k, s_tend] -= Δσₖ⁻¹ * (w⁺ * ξᶠ⁺ - w⁻ * ξᶠ⁻ - ξ[ij, k, s_prog] * (w⁺ - w⁻))

            face[ij] = ξᶠ⁺
        end
    end
    return nothing
end

function _vertical_advection!(arch::GPU, ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv, face::AbstractField)
    # worksize is the horizontal × vertical iteration space (skip the step dimension)
    launch!(
        arch, RingGridWorkOrder, (size(ξ_tend, 1), size(ξ_tend, 2)),
        vertical_advection_kernel!,
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv
    )
    return nothing
end

@kernel inbounds = true function vertical_advection_kernel!(
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv
    )
    ij, k = @index(Global, NTuple)

    Δσₖ⁻¹ = inv(Δσ[k])

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
    S⁺ = gather_stencil_values(ξ, ij, s_prog, Base.tail(k_stencil))
    S⁻ = gather_stencil_values(ξ, ij, s_prog, Base.front(k_stencil))

    ξᶠ⁺ = reconstruct_face(S⁺, w⁺, adv)
    ξᶠ⁻ = reconstruct_face(S⁻, w⁻, adv)

    # -= as the tendencies already contain the parameterizations
    ξ_tend[ij, k,s_tend] -= Δσₖ⁻¹ * (w⁺ * ξᶠ⁺ - w⁻ * ξᶠ⁻ - ξ[ij, k, s_prog] * (w⁺ - w⁻))
end

# gathers the stencil values ξ[ij, k[i], s] into a register tuple, once, so that
# `reconstruct_face` below operates purely on values: this lets the same reconstruction
# code be reused by memory-access patterns other than per-(ij,k) array indexing (e.g. a
# CPU column loop keeping a sliding window, or a GPU per-column kernel).
@inline function gather_stencil_values(ξ, ij, s, k::NTuple{N, Integer}) where {N}
    return ntuple(i -> ξ[ij, k[i], s], N)
end

# reconstruct_face operates on `S`, the gathered stencil values (one half of `k_stencil`,
# in the same index order as before): `S[1] == ξ[ij, k[1], s]` etc.

# 1st order upwind
@inline reconstruct_face(S, u, ::UpwindVerticalAdvection{NF, 1}) where {NF} =
    ifelse(u > 0, S[1], S[2])

# 3rd order upwind
@inline reconstruct_face(S, u, ::UpwindVerticalAdvection{NF, 2}) where {NF} =
    ifelse(
    u > 0, (2S[1] + 5S[2] - S[3]) * 1 // 6,
    (2S[4] + 5S[3] - S[2]) * 1 // 6
)

# 5th order upwind
@inline reconstruct_face(S, u, ::UpwindVerticalAdvection{NF, 3}) where {NF} =
    ifelse(
    u > 0, (2S[1] - 13S[2] + 47S[3] + 27S[4] - 3S[5]) * 1 // 60,
    (2S[6] - 13S[5] + 47S[4] + 27S[3] - 3S[2]) * 1 // 60
)

# 2nd order centered
@inline reconstruct_face(S, u, ::CenteredVerticalAdvection{NF, 1}) where {NF} =
    (S[1] + S[2]) * 1 // 2

# 4th order centered
@inline reconstruct_face(S, u, ::CenteredVerticalAdvection{NF, 2}) where {NF} =
    (-S[1] + 7S[2] + 7S[3] - S[4]) * 1 // 12

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

@inline function reconstruct_face(S, u, ::WENOVerticalAdvection)
    if u > 0
        S₀ = (S[3], S[4], S[5])
        S₁ = (S[2], S[3], S[4])
        S₂ = (S[1], S[2], S[3])
    else
        S₀ = (S[4], S[3], S[2])
        S₁ = (S[5], S[4], S[3])
        S₂ = (S[6], S[5], S[4])
    end
    return weno_reconstruction(S₀, S₁, S₂)
end
