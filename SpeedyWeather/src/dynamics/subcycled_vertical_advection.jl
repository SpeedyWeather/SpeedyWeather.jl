export SubcycledVerticalAdvection

"""
Courant-adaptive, column-local subcycled vertical advection. The other vertical advection schemes
are explicit with one evaluation per time step and become unstable for vertical Courant numbers
|σ̇| Δt / Δσ ≳ 1, which with many thin layers (L100: Δσ ≈ 0.001-0.01) is easily exceeded in strong
ascent, e.g. in deep convection or tropical cyclones, leading to a 2Δz checkerboard blow-up.

This scheme advects from the same time step as the other diffusive schemes over the interval of
the current time step (`time_step(time_stepping, clock)`, i.e. 2Δt for leapfrog), but splits that
interval per column into `n` substeps such that the Courant number of every substep is at most
`courant_target`. Each substep is a 3-stage strong-stability-preserving Runge-Kutta step
(Shu and Osher, 1988) of the WENO5 flux-form operator of [`WENOVerticalAdvection`](@ref).
The tendency is the mean over the interval, `(ξ_end - ξ₀) / interval`, added to the tendencies.
Columns with weak vertical motion take `n = 1` (3 operator evaluations), `n` is at most `max_substeps`.

Work arrays of size (npoints, nlayers) are allocated three times at construction.
Fields are $(TYPEDFIELDS)"""
struct SubcycledVerticalAdvection{NF, A, R} <: DiffusiveVerticalAdvection{NF, 3}
    "[OPTION] Maximum Courant number per substep (SSP-RK3 with WENO5 is stable to about 1.4)"
    courant_target::NF

    "[OPTION] Maximum number of substeps per column and time step"
    max_substeps::Int

    "Work arrays for the Runge-Kutta stages"
    S::A
    A::A
    B::A

    "Integration interval [s/radius], set every time step"
    interval::R
end

function SubcycledVerticalAdvection(spectral_grid::SpectralGrid; courant_target = 1, max_substeps = 16)
    (; NF, npoints, nlayers, architecture) = spectral_grid
    work() = on_architecture(architecture, zeros(NF, npoints, nlayers))
    S, A, B = work(), work(), work()
    return SubcycledVerticalAdvection{NF, typeof(S), Base.RefValue{NF}}(
        NF(courant_target), max_substeps, S, A, B, Ref(zero(NF))
    )
end

function Base.show(io::IO, adv::SubcycledVerticalAdvection{NF}) where {NF}
    return print(io, "SubcycledVerticalAdvection{$NF}(courant_target = $(adv.courant_target), max_substeps = $(adv.max_substeps))")
end

# per-variable entry point: set the integration interval of this time step (Δt/2, Δt, then 2Δt for leapfrog),
# w is radius-scaled σ̇, so the interval is divided by the radius
@inline function vertical_advection!(
        ::Val{var}, vars::Variables, w, Δσ, adv::SubcycledVerticalAdvection, model,
    ) where {var}
    haskey(vars.tendencies.grid, var) || return nothing
    adv.interval[] = time_step(model.time_stepping, vars.prognostic.clock) / model.planet.radius
    ξ_tend = vars.tendencies.grid[var]
    ξ = vars.grid[var]
    s_tend = which_tendency_step(ξ_tend, model.time_stepping, adv)
    s_prog = which_prognostic_step(ξ, model.time_stepping, adv)
    return _vertical_advection!(ξ_tend, s_tend, w, ξ, s_prog, Δσ, adv)
end

function _vertical_advection!(
        ξ_tend::AbstractField,      # tendency of quantity ξ (full, with step dimension)
        s_tend::Integer,            # step of ξ_tend to write into
        w::AbstractField,           # vertical velocity at k+1/2
        ξ::AbstractField,           # ξ (full, with step dimension)
        s_prog::Integer,            # step of ξ to advect
        Δσ,                         # layer thickness on σ levels
        adv::SubcycledVerticalAdvection,
    )
    grids_match(ξ_tend, w, ξ) || throw(DimensionMismatch(ξ_tend, w, ξ))
    npoints, nlayers = size(ξ, 1), size(ξ, 2)
    launch!(
        architecture(ξ_tend), LinearWorkOrder, (npoints,), subcycled_vertical_advection_kernel!,
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, adv.S, adv.A, adv.B,
        adv.interval[], adv.courant_target, adv.max_substeps, WENOVerticalAdvection{eltype(adv.S)}(),
    )
    return nothing
end

# flux-form advection operator L(X) at layer k, same discretisation as `vertical_advection_kernel!`
@inline function column_advection(X, ij, k, w, Δσ, nlayers, weno)
    w⁻ = w[ij, max(1, k - 1)]
    w⁻ = ifelse(k == 1, zero(w⁻), w⁻)              # σ̇ = 0 at the model top
    w⁺ = w[ij, k]                                  # σ̇ = 0 at the surface (w[:, nlayers] = 0)
    stencil = retrieve_stencil(k, nlayers, weno)
    Xᶠ⁺ = reconstructed_at_face(X, ij, 1, Base.tail(stencil), w⁺, weno)
    Xᶠ⁻ = reconstructed_at_face(X, ij, 1, Base.front(stencil), w⁻, weno)
    return -(w⁺ * Xᶠ⁺ - w⁻ * Xᶠ⁻ - X[ij, k] * (w⁺ - w⁻)) / Δσ[k]
end

@kernel inbounds = true function subcycled_vertical_advection_kernel!(
        ξ_tend, s_tend, w, ξ, s_prog, Δσ, nlayers, S, A, B, interval, courant_target, max_substeps, weno,
    )
    ij = @index(Global, Linear)
    NF = eltype(S)

    # column Courant number over the whole interval → number of substeps
    C = zero(NF)
    for k in 1:(nlayers - 1)
        C = max(C, abs(w[ij, k]) * interval / min(Δσ[k], Δσ[k + 1]))
    end
    # clamp in floating point first, ceil(Int, NaN) would throw (on GPUs an opaque KernelException);
    # a non-finite column is left to propagate so that the NaN detection reports it
    r = C / courant_target
    r = ifelse(isfinite(r), r, NF(max_substeps))
    n = unsafe_trunc(Int, clamp(ceil(r), one(NF), NF(max_substeps)))
    h = interval / n

    for k in 1:nlayers
        S[ij, k] = ξ[ij, k, s_prog]
    end

    for _ in 1:n
        # SSP-RK3 (Shu-Osher): A = S + hL(S); B = 3/4 S + 1/4 (A + hL(A)); S = 1/3 S + 2/3 (B + hL(B))
        for k in 1:nlayers
            A[ij, k] = S[ij, k] + h * column_advection(S, ij, k, w, Δσ, nlayers, weno)
        end
        for k in 1:nlayers
            B[ij, k] = 3 * S[ij, k] / 4 + (A[ij, k] + h * column_advection(A, ij, k, w, Δσ, nlayers, weno)) / 4
        end
        for k in 1:nlayers
            S[ij, k] = S[ij, k] / 3 + 2 * (B[ij, k] + h * column_advection(B, ij, k, w, Δσ, nlayers, weno)) / 3
        end
    end

    # mean tendency over the interval, added to the (parameterization) tendencies
    for k in 1:nlayers
        ξ_tend[ij, k, s_tend] += (S[ij, k] - ξ[ij, k, s_prog]) / interval
    end
end
