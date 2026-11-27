export SubcycledVerticalAdvection

"""
Column-local subcycled vertical advection. The other vertical advection schemes are explicit with one
evaluation per time step and become unstable for vertical Courant numbers |σ̇| Δt / Δσ ≳ 1. With many thin
layers (e.g. L100 with Δσ ≈ 0.001-0.01 near the surface) this is exceeded routinely in ascent, e.g. in deep
convection or tropical cyclones, leading to a 2Δz checkerboard blow-up.

This scheme splits the interval of the current time step (`time_step(time_stepping, clock)`, i.e. Δt/2, Δt
and then 2Δt for leapfrog) into `n` substeps per column. Each substep is a 3-stage strong-stability-preserving
Runge-Kutta step (Shu and Osher, 1988) of the flux-form operator of the inner `scheme` (default
[`WENOVerticalAdvection`](@ref)). The tendency is the mean over the interval, `(ξ_end - ξ₀) / interval`,
added to the tendencies.

- With `courant_target::Real` (default 1), `n` is the smallest number of substeps for which the Courant number
  of every substep is at most `courant_target`, clamped to `1 ≤ n ≤ substeps`. Columns with weak vertical
  motion take `n = 1`. SSP-RK3 is linearly stable up to a Courant number of ≈1.43 with WENO5 (its linear,
  fifth-order upwind operator) and √3 ≈ 1.73 with centred operators.
- With `courant_target = nothing`, every column takes exactly `substeps` substeps. The runtime is then
  independent of the flow.

Like all `DiffusiveVerticalAdvection` schemes, it advects ξ of the *old* time step (`which_prognostic_step`,
t-Δt for leapfrog) over the whole interval, which is consistent with how leapfrog uses diffusive tendencies.
Vertical advection is therefore effectively operator-split from the other tendencies.

Uses the scratch arrays `vars.scratch.grid.a`, `.b` and `.c` for the Runge-Kutta stages.
Fields are $(TYPEDFIELDS)"""
struct SubcycledVerticalAdvection{NF, B, S <: DiffusiveVerticalAdvection{NF, B}, C <: Union{NF, Nothing}} <:
    DiffusiveVerticalAdvection{NF, B}
    "[OPTION] Vertical advection scheme of each Runge-Kutta stage"
    scheme::S

    "[OPTION] Maximum Courant number per substep, or `nothing` for a fixed number of `substeps`"
    courant_target::C

    "[OPTION] Maximum number of substeps per column and time step, or the fixed number for `courant_target = nothing`"
    substeps::Int
end

function SubcycledVerticalAdvection(
        spectral_grid::SpectralGrid;
        scheme::DiffusiveVerticalAdvection = WENOVerticalAdvection(spectral_grid),
        courant_target = 1,
        substeps::Integer = 16,
    )
    (; NF) = spectral_grid
    substeps ≥ 1 || throw(ArgumentError("substeps = $substeps has to be ≥ 1"))
    C = isnothing(courant_target) ? nothing : convert(NF, courant_target)
    B = stencil_halfwidth(scheme)
    return SubcycledVerticalAdvection{NF, B, typeof(scheme), typeof(C)}(scheme, C, substeps)
end

stencil_halfwidth(::VerticalAdvection{NF, B}) where {NF, B} = B

function Base.show(io::IO, adv::SubcycledVerticalAdvection{NF}) where {NF}
    return print(
        io, "SubcycledVerticalAdvection{$NF}(scheme = $(nameof(typeof(adv.scheme))), ",
        "courant_target = $(adv.courant_target), substeps = $(adv.substeps))"
    )
end

# Runge-Kutta stage arrays, shared with other components that use the same scratch variables
variables(::SubcycledVerticalAdvection) = (
    ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
    ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
    ScratchVariable(:c, GridXYZ(), desc = "Scratch array", namespace = :grid),
)

# integration interval [s/radius] of this time step (Δt/2, Δt, then 2Δt for leapfrog),
# w is the radius-scaled σ̇, so the interval is divided by the radius.
# Other schemes don't need it (`nothing`), see `vertical_advection!(vars, model)`
advection_interval(::SubcycledVerticalAdvection, vars, model) =
    time_step(model.time_stepping, vars.prognostic.clock) / model.planet.radius

function _vertical_advection!(
        ξ_tend::AbstractField,      # tendency of quantity ξ (full, with step dimension)
        s_tend::Integer,            # step of ξ_tend to write into
        w::AbstractField,           # vertical velocity at k+1/2
        ξ::AbstractField,           # ξ (full, with step dimension)
        s_prog::Integer,            # step of ξ to advect
        Δσ,                         # layer thickness on σ levels
        adv::SubcycledVerticalAdvection,
        interval,                   # integration interval [s/radius] from `advection_interval`
        scratch,                    # scratch variables, uses scratch.grid.a, .b, .c
    )
    grids_match(ξ_tend, w, ξ) || throw(DimensionMismatch(ξ_tend, w, ξ))
    (; a, b, c) = scratch.grid
    npoints, nlayers = size(ξ, 1), size(ξ, 2)
    launch!(
        architecture(ξ_tend), LinearWorkOrder, (npoints,), subcycled_vertical_advection_kernel!,
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, a.data, b.data, c.data,
        convert(eltype(a), interval), adv.courant_target, adv.substeps, adv.scheme,
    )
    return nothing
end

# flux-form advection operator L(X) at layer k, same discretisation as `vertical_advection_kernel!`
@inline function column_advection(X, ij, k, w, Δσ, nlayers, scheme)
    w⁻ = w[ij, max(1, k - 1)]
    w⁻ = ifelse(k == 1, zero(w⁻), w⁻)              # σ̇ = 0 at the model top
    w⁺ = w[ij, k]                                  # σ̇ = 0 at the surface (w[:, nlayers] = 0)
    stencil = retrieve_stencil(k, nlayers, scheme)
    Xᶠ⁺ = reconstructed_at_face(X, ij, 1, Base.tail(stencil), w⁺, scheme)
    Xᶠ⁻ = reconstructed_at_face(X, ij, 1, Base.front(stencil), w⁻, scheme)
    return -(w⁺ * Xᶠ⁺ - w⁻ * Xᶠ⁻ - X[ij, k] * (w⁺ - w⁻)) / Δσ[k]
end

# number of substeps: fixed
@inline number_of_substeps(ij, w, Δσ, nlayers, interval, ::Nothing, substeps) = substeps

# number of substeps: from the column's maximum Courant number over the whole interval
@inline function number_of_substeps(ij, w, Δσ, nlayers, interval, courant_target, substeps)
    NF = typeof(courant_target)
    C = zero(NF)
    for k in 1:(nlayers - 1)
        C = max(C, abs(w[ij, k]) * interval / min(Δσ[k], Δσ[k + 1]))
    end
    # clamp in floating point first, ceil(Int, NaN) would throw (on GPUs an opaque KernelException);
    # a non-finite column is left to propagate so that the NaN detection reports it
    r = C / courant_target
    r = ifelse(isfinite(r), r, NF(substeps))
    return unsafe_trunc(Int, clamp(ceil(r), one(NF), NF(substeps)))
end

@kernel inbounds = true function subcycled_vertical_advection_kernel!(
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, S, A, B, interval, courant_target, substeps, scheme,
    )
    ij = @index(Global, Linear)

    n = number_of_substeps(ij, w, Δσ, nlayers, interval, courant_target, substeps)
    h = interval / n

    for k in 1:nlayers
        S[ij, k] = ξ[ij, k, s_prog]
    end

    for _ in 1:n
        # SSP-RK3 (Shu-Osher): A = S + hL(S); B = 3/4 S + 1/4 (A + hL(A)); S = 1/3 S + 2/3 (B + hL(B))
        for k in 1:nlayers
            A[ij, k] = S[ij, k] + h * column_advection(S, ij, k, w, Δσ, nlayers, scheme)
        end
        for k in 1:nlayers
            B[ij, k] = 3 * S[ij, k] / 4 + (A[ij, k] + h * column_advection(A, ij, k, w, Δσ, nlayers, scheme)) / 4
        end
        for k in 1:nlayers
            S[ij, k] = S[ij, k] / 3 + 2 * (B[ij, k] + h * column_advection(B, ij, k, w, Δσ, nlayers, scheme)) / 3
        end
    end

    # mean tendency over the interval, added to the (parameterization) tendencies
    for k in 1:nlayers
        ξ_tend[ij, k, s_tend] += (S[ij, k] - ξ[ij, k, s_prog]) / interval
    end
end
