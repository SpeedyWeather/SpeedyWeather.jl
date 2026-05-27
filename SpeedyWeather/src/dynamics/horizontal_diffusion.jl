abstract type AbstractHorizontalDiffusion <: AbstractModelComponent end

export HyperDiffusion

"""
Horizontal hyper diffusion of vorticity, div, temp, humid; implicitly in spectral space
with a `power` of the Laplacian (default = 4) and the strength controlled by
`time_scale` (default = 1 hour). For vorticity and divergence, by default,
the `time_scale` (=1/strength of diffusion) is reduced with increasing resolution
through `resolution_scaling` and the power is linearly decreased in the vertical
above the `tapering_σ` sigma level to `power_stratosphere` (default 2).

For the BarotropicModel and ShallowWaterModel no tapering or scaling is applied.
Fields and options are
$(TYPEDFIELDS)"""
@kwdef mutable struct HyperDiffusion{
        NF,
        MatrixType,
        IntType,
        S,
    } <: AbstractHorizontalDiffusion

    # DIMENSIONS
    "spectral resolution"
    trunc::IntType

    "number of vertical levels"
    nlayers::IntType

    # PARAMETERS
    "[OPTION] power of Laplacian"
    power::NF = 4

    "[OPTION] diffusion time scale"
    time_scale::S = Hour(4)

    "[OPTION] diffusion time scale for divergence"
    time_scale_div::S = Hour(1)

    "[OPTION] stronger diffusion with resolution? 0: constant with trunc, 1: (inverse) linear with trunc, etc"
    resolution_scaling::NF = 1

    # incrased diffusion in stratosphere
    "[OPTION] different power for tropopause/stratosphere"
    power_stratosphere::NF = 2

    "[OPTION] linearly scale towards power_stratosphere above this σ"
    tapering_σ::NF = 0.2

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    expl::MatrixType = zeros(NF, trunc + 2, nlayers)      # explicit part
    impl::MatrixType = ones(NF, trunc + 2, nlayers)       # implicit part

    # ARRAYS using time_scale_div
    expl_div::MatrixType = zeros(NF, trunc + 2, nlayers)    # explicit part
    impl_div::MatrixType = ones(NF, trunc + 2, nlayers)     # implicit part

    # Pre-tiled damping matrices matched to the :prognostic fuse parent slot layout
    # `(vor, div, T, pres, [humid])`. Sized for PrimitiveWet (`4·nlayers + 1`); PrimitiveDry
    # uses the leading `3·nlayers + 1` slots. Slot order matches `vars.fused.prognostic.slot_map`:
    #
    #   slots 1..L            : vorticity   — `expl, impl`
    #   slots L+1..2L         : divergence  — `expl_div, impl_div`  (stronger damping)
    #   slots 2L+1..3L        : temperature — `expl, impl`
    #   slot  3L+1            : pressure    — no-op damping (0 / 1)
    #   slots 3L+2..4L+1      : humidity    — `expl, impl`          (wet only)
    #
    # Used by the batched `horizontal_diffusion!(::Variables, …, ::PrimitiveEquation, lf)`
    # kernel which reads from the prognostic fuse parent at the corresponding slots.
    expl_all::MatrixType = zeros(NF, trunc + 2, 4 * nlayers + 1)
    impl_all::MatrixType = ones(NF, trunc + 2, 4 * nlayers + 1)
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function HyperDiffusion(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlayers, ArrayType) = spectral_grid        # take resolution parameters from spectral_grid
    MatrixType = ArrayType{NF, 2}
    return HyperDiffusion{NF, MatrixType, typeof(trunc), Dates.Second}(; trunc, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms in `diffusion` based on the
model time step, and possibly with a changing strength/power in
the vertical."""
function initialize!(
        diffusion::HyperDiffusion,
        model::AbstractModel
    )
    return initialize!(diffusion, model.geometry, model.time_stepping)
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms for all layers based on the
model time step in `L`, the vertical level sigma level in `G`."""
function initialize!(
        diffusion::HyperDiffusion,
        G::AbstractGeometry,
        L::AbstractTimeStepper,
    )
    (; trunc, nlayers, resolution_scaling) = diffusion
    (; power, power_stratosphere, tapering_σ) = diffusion
    (; Δt, radius) = L

    # Reduce diffusion time scale (=increase diffusion, always in seconds) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = Second(diffusion.time_scale).value / radius * (32 / (trunc + 1))^resolution_scaling
    time_scale_div = Second(diffusion.time_scale_div).value / radius * (32 / (trunc + 1))^resolution_scaling

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -trunc * (trunc + 1)

    # Get architecture and arrays
    ∇²ⁿ = diffusion.expl
    ∇²ⁿ_implicit = diffusion.impl
    ∇²ⁿ_div = diffusion.expl_div
    ∇²ⁿ_div_implicit = diffusion.impl_div
    σ_levels_full = G.σ_levels_full

    # Launch kernel
    arch = architecture(∇²ⁿ)
    worksize = (trunc + 2, nlayers)
    launch!(
        arch, Array3DWorkOrder, worksize, _initialize_hyperdiffusion_kernel!,
        ∇²ⁿ, ∇²ⁿ_implicit, ∇²ⁿ_div, ∇²ⁿ_div_implicit, σ_levels_full,
        trunc, power, power_stratosphere, tapering_σ,
        time_scale, time_scale_div, Δt, largest_eigenvalue
    )

    _tile_diffusion_all!(diffusion)
    return nothing
end

@kernel inbounds = true function _initialize_hyperdiffusion_kernel!(
        ∇²ⁿ,
        ∇²ⁿ_implicit,
        ∇²ⁿ_div,
        ∇²ⁿ_div_implicit,
        @Const(σ_levels_full),
        trunc,
        power,
        power_stratosphere,
        tapering_σ,
        time_scale,
        time_scale_div,
        Δt,
        largest_eigenvalue
    )
    l_plus_1, k = @index(Global, NTuple)  # l+1 index (1-based), layer index

    l = l_plus_1 - 1  # actual degree l (0-based)

    # last degree is only used by vector quantities; set to zero for implicit and explicit
    # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
    if l_plus_1 == trunc + 2
        ∇²ⁿ[l_plus_1, k] = 0
        ∇²ⁿ_implicit[l_plus_1, k] = 0
        ∇²ⁿ_div[l_plus_1, k] = 0
        ∇²ⁿ_div_implicit[l_plus_1, k] = 0
    else
        # VERTICAL TAPERING for the stratosphere
        # go from 1 to 0 between σ=0 and tapering_σ
        σ = σ_levels_full[k]
        tapering = max(0, (tapering_σ - σ) / tapering_σ)  # ∈ [0, 1]
        p = power + tapering * (power_stratosphere - power)

        # Normalized eigenvalue
        eigenvalue_norm = -l * (l + 1) / largest_eigenvalue

        # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
        ∇²ⁿ[l_plus_1, k] = -eigenvalue_norm^power / time_scale
        ∇²ⁿ_div[l_plus_1, k] = -eigenvalue_norm^p / time_scale_div

        # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
        ∇²ⁿ_implicit[l_plus_1, k] = 1 / (1 - 2Δt * ∇²ⁿ[l_plus_1, k])
        ∇²ⁿ_div_implicit[l_plus_1, k] = 1 / (1 - 2Δt * ∇²ⁿ_div[l_plus_1, k])
    end
end

"""
$(TYPEDSIGNATURES)
Tile `diffusion.expl_all`, `diffusion.impl_all` (shape `(lmax, 4·nlayers + 1)`) from the
per-variable matrices so that they line up slot-for-slot with the `:prognostic` fuse
parent. The actual declaration order on PrimitiveDry/Wet places members as
`(vor, div, T, pres, [humid])` (pressure before humidity because pressure is declared
in PrimitiveDry and humidity is appended by PrimitiveWet), so slot ranges are:

- `1..L`            = vorticity     → `expl, impl`
- `L+1..2L`         = divergence    → `expl_div, impl_div`  (stronger damping)
- `2L+1..3L`        = temperature   → `expl, impl`
- `3L+1`            = pressure      → no-op damping (0 / 1) — diffusion never touches pressure
- `3L+2..4L+1`      = humidity      → `expl, impl`  (only consumed by PrimitiveWet)
"""
function _tile_diffusion_all!(diffusion::AbstractHorizontalDiffusion)
    L = diffusion.nlayers
    expl_all = diffusion.expl_all
    impl_all = diffusion.impl_all
    fill!(expl_all, 0)
    fill!(impl_all, 1)
    # vorticity (slots 1..L)
    copyto!(view(expl_all, :, 1:L),                diffusion.expl)
    copyto!(view(impl_all, :, 1:L),                diffusion.impl)
    # divergence (slots L+1..2L) — stronger damping
    copyto!(view(expl_all, :, (L + 1):(2L)),       diffusion.expl_div)
    copyto!(view(impl_all, :, (L + 1):(2L)),       diffusion.impl_div)
    # temperature (slots 2L+1..3L)
    copyto!(view(expl_all, :, (2L + 1):(3L)),      diffusion.expl)
    copyto!(view(impl_all, :, (2L + 1):(3L)),      diffusion.impl)
    # pressure (slot 3L+1) — stays no-op (expl=0, impl=1) from the fill above
    # humidity (slots 3L+2..4L+1) — only used by PrimitiveWet, harmless for dry
    copyto!(view(expl_all, :, (3L + 2):(4L + 1)),  diffusion.expl)
    copyto!(view(impl_all, :, (3L + 2):(4L + 1)),  diffusion.impl)
    return nothing
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to a 2D field `var` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `expl` and an implicit part `impl`, such that both can be calculated
in a single forward step by using `var` as well as its tendency `tendency`."""
function horizontal_diffusion!(
        tendency::LowerTriangularArray,     # tendency of a
        var::LowerTriangularArray,          # spectral horizontal field to diffuse
        expl::AbstractMatrix,               # explicit spectral damping (lmax x nlayers matrix)
        impl::AbstractMatrix,               # implicit spectral damping (lmax x nlayers matrix)
    )
    lmax, mmax = size(tendency, OneBased, as = Matrix)
    nlayers = size(var, 2)

    @boundscheck size(tendency) == size(var) || throw(BoundsError(tendency))
    @boundscheck lmax <= size(expl, 1) == size(impl, 1) || throw(BoundsError(expl, lmax))
    @boundscheck nlayers <= size(expl, 2) == size(impl, 2) || throw(BoundsError(expl, nlayers))

    launch!(
        architecture(tendency), SpectralWorkOrder, size(tendency), _horizontal_diffusion_kernel!,
        tendency, var, expl, impl, var.spectrum.l_indices
    )
    return nothing
end

@kernel inbounds = true function _horizontal_diffusion_kernel!(
        tendency, var, expl, impl, l_indices
    )

    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(var) == 1 ? 1 : I[2]

    # Get the degree l for this coefficient
    l = l_indices[lm]

    # Apply horizontal diffusion
    tendency[I] = (tendency[I] + expl[l, k] * var[I]) * impl[l, k]
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity in the BarotropicModel."""
function horizontal_diffusion!(
        vars::Variables,
        diffusion::AbstractHorizontalDiffusion,
        model::Barotropic,
        lf::Integer = 1,    # leapfrog index used (2 is unstable)
    )
    (; expl, impl) = diffusion

    # Barotropic model diffuses vorticity (only variable)
    vor = get_step(vars.prognostic.vorticity, lf)                               # lta_view for leapfrog index
    vor_tend = vars.tendencies.vorticity
    horizontal_diffusion!(vor_tend, vor, expl, impl)

    for (name, tracer) in model.tracers
        tracer_var = get_step(vars.prognostic.tracers[name], lf)          # lta_view for leapfrog index
        tracer_tend = vars.tendencies.tracers[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity and divergence in the ShallowWaterModel."""
function horizontal_diffusion!(
        vars::Variables,
        diffusion::AbstractHorizontalDiffusion,
        model::ShallowWater,
        lf::Integer = 1,    # leapfrog index used (2 is unstable)
    )
    (; expl, impl, expl_div, impl_div) = diffusion

    # ShallowWater model diffuses vorticity and divergence
    vor = get_step(vars.prognostic.vorticity, lf)
    div = get_step(vars.prognostic.divergence, lf)
    vor_tend = vars.tendencies.vorticity
    div_tend = vars.tendencies.divergence
    horizontal_diffusion!(vor_tend, vor, expl, impl)
    horizontal_diffusion!(div_tend, div, expl_div, impl_div)

    for (name, tracer) in model.tracers
        tracer_var = get_step(vars.prognostic.tracers[name], lf)      # lta_view for leapfrog index
        tracer_tend = vars.tendencies.tracers[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity, divergence, temperature, and humidity
(PrimitiveWet only) in the PrimitiveEquation models.

Architecture-dispatched: on GPU, a single batched kernel reads from the `:prognostic`
fuse parent and writes to the scattered final tendency arrays in one launch, using the
pre-tiled `expl_all` / `impl_all` damping matrices that carry the right per-variable
damping at each slot. On CPU, the per-launch overhead is negligible while the wider
inner-loop memory footprint of the batched kernel measurably hurts at higher resolutions
(see `benchmark_horizontal_diffusion.jl`), so CPU stays on the per-variable code path.

Tracers stay on the single-variable code path regardless (one launch per active tracer)."""
function horizontal_diffusion!(
        vars::Variables,
        diffusion::AbstractHorizontalDiffusion,
        model::PrimitiveEquation,
        lf::Integer = 1,    # leapfrog index used (2 is unstable)
    )
    arch = architecture(vars.tendencies.vorticity)
    if arch isa Architectures.AbstractCPU
        _horizontal_diffusion_primeq_serial!(vars, diffusion, model, lf)
    else
        _horizontal_diffusion_primeq_batched!(vars, diffusion, model, lf)
    end

    # tracers all use the standard (expl, impl) damping; one launch per active tracer
    (; expl, impl) = diffusion
    for (name, tracer) in model.tracers
        tracer_var = get_step(vars.prognostic.tracers[name], lf)
        tracer_tend = vars.tendencies.tracers[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
CPU code path: one `horizontal_diffusion!(tendency, var, expl, impl)` launch per
prognostic variable. The per-launch cost on CPU is negligible while the smaller inner
loop (one read + one write per `(lm, k)`) vectorizes well — at T127–T255 this is
~1.4–2× faster than the batched code path."""
function _horizontal_diffusion_primeq_serial!(vars, diffusion, model, lf)
    (; expl, impl, expl_div, impl_div) = diffusion
    vor   = get_step(vars.prognostic.vorticity,   lf)
    div   = get_step(vars.prognostic.divergence,  lf)
    temp  = get_step(vars.prognostic.temperature, lf)
    horizontal_diffusion!(vars.tendencies.vorticity,   vor,  expl,     impl)
    horizontal_diffusion!(vars.tendencies.divergence,  div,  expl_div, impl_div)
    horizontal_diffusion!(vars.tendencies.temperature, temp, expl,     impl)
    if haskey(vars.tendencies, :humidity)
        humid = get_step(vars.prognostic.humidity, lf)
        horizontal_diffusion!(vars.tendencies.humidity, humid, expl, impl)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
GPU code path: a single batched kernel reads sources from the `:prognostic` fuse parent
at leapfrog step `lf` (one contiguous buffer covering `vor, div, T, pres, [q]`) and
writes to the scattered final tendency arrays. The pre-tiled `expl_all` / `impl_all`
matrices (slot-aligned with `:prognostic`) carry the right per-variable damping at each
slot, so the kernel uses uniform indexing. Collapses 3 (dry) / 4 (wet) launches into 1."""
function _horizontal_diffusion_primeq_batched!(vars, diffusion, model, lf)
    (; expl_all, impl_all) = diffusion
    prog_parent = parent(vars.fused.prognostic)
    prog_step = get_step(prog_parent, lf)               # 2D LTA, (lm, total_slots)
    slot_map = vars.fused.prognostic.slot_map

    vor_tend  = vars.tendencies.vorticity
    div_tend  = vars.tendencies.divergence
    temp_tend = vars.tendencies.temperature

    lmax, _ = size(vor_tend, OneBased, as = Matrix)
    @boundscheck size(vor_tend) == size(div_tend) == size(temp_tend) ||
        throw(BoundsError(vor_tend))
    @boundscheck lmax <= size(expl_all, 1) == size(impl_all, 1) ||
        throw(BoundsError(expl_all, lmax))

    vor_offset  = first(slot_map.vorticity)   - 1
    div_offset  = first(slot_map.divergence)  - 1
    temp_offset = first(slot_map.temperature) - 1
    l_indices = vor_tend.spectrum.l_indices
    arch = architecture(vor_tend)

    if haskey(vars.tendencies, :humidity)
        humid_tend = vars.tendencies.humidity
        humid_offset = first(slot_map.humidity) - 1
        @boundscheck size(humid_tend) == size(vor_tend) || throw(BoundsError(humid_tend))
        launch!(
            arch, SpectralWorkOrder, size(vor_tend), _horizontal_diffusion_primeq_wet_kernel!,
            vor_tend, div_tend, temp_tend, humid_tend,
            prog_step, expl_all, impl_all, l_indices,
            vor_offset, div_offset, temp_offset, humid_offset,
        )
    else
        launch!(
            arch, SpectralWorkOrder, size(vor_tend), _horizontal_diffusion_primeq_dry_kernel!,
            vor_tend, div_tend, temp_tend,
            prog_step, expl_all, impl_all, l_indices,
            vor_offset, div_offset, temp_offset,
        )
    end
    return nothing
end

@kernel inbounds = true function _horizontal_diffusion_primeq_dry_kernel!(
        vor_tend, div_tend, temp_tend,
        prog_step, expl_all, impl_all, l_indices,
        vor_offset, div_offset, temp_offset,
    )
    I = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]
    l  = l_indices[lm]

    s = vor_offset + k
    vor_tend[lm, k]  = (vor_tend[lm, k]  + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]

    s = div_offset + k
    div_tend[lm, k]  = (div_tend[lm, k]  + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]

    s = temp_offset + k
    temp_tend[lm, k] = (temp_tend[lm, k] + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]
end

@kernel inbounds = true function _horizontal_diffusion_primeq_wet_kernel!(
        vor_tend, div_tend, temp_tend, humid_tend,
        prog_step, expl_all, impl_all, l_indices,
        vor_offset, div_offset, temp_offset, humid_offset,
    )
    I = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]
    l  = l_indices[lm]

    s = vor_offset + k
    vor_tend[lm, k]   = (vor_tend[lm, k]   + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]

    s = div_offset + k
    div_tend[lm, k]   = (div_tend[lm, k]   + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]

    s = temp_offset + k
    temp_tend[lm, k]  = (temp_tend[lm, k]  + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]

    s = humid_offset + k
    humid_tend[lm, k] = (humid_tend[lm, k] + expl_all[l, s] * prog_step[lm, s]) * impl_all[l, s]
end

export SpectralFilter

"""Spectral filter for horizontal diffusion. Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SpectralFilter{
        NF,
        MatrixType,
        IntType,
    } <: AbstractHorizontalDiffusion

    # DIMENSIONS
    "spectral resolution"
    trunc::IntType

    "number of vertical levels"
    nlayers::IntType

    # PARAMETERS
    "[OPTION] shift diffusion to higher (positive shift) or lower (neg) wavenumbers, relative to trunc"
    shift::NF = 0

    "[OPTION] Scale-selectiveness, steepness of the sigmoid, higher is more selective"
    scale::NF = 0.05

    "[OPTION] diffusion time scale"
    time_scale::Second = Hour(4)

    "[OPTION] stronger diffusion time scale for divergence"
    time_scale_div::Second = Minute(30)

    "[OPTION] resolution scaling to shorten time_scale with trunc"
    resolution_scaling::NF = 1

    "[OPTION] power of the tanh function"
    power::NF = 4

    "[OPTION] power of the tanh function for divergence"
    power_div::NF = 4

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    expl::MatrixType = zeros(NF, trunc + 2, nlayers)  # explicit part
    impl::MatrixType = ones(NF, trunc + 2, nlayers)   # implicit part

    # ARRAYS using time_scale_div for divergence
    expl_div::MatrixType = zeros(NF, trunc + 2, nlayers)  # explicit part
    impl_div::MatrixType = ones(NF, trunc + 2, nlayers)   # implicit part

    # Pre-tiled damping aligned with the :prognostic fuse parent — see HyperDiffusion comment.
    expl_all::MatrixType = zeros(NF, trunc + 2, 4 * nlayers + 1)
    impl_all::MatrixType = ones(NF, trunc + 2, 4 * nlayers + 1)
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function SpectralFilter(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlayers, ArrayType) = spectral_grid        # take resolution parameters from spectral_grid
    MatrixType = ArrayType{NF, 2}
    return SpectralFilter{NF, MatrixType, typeof(trunc)}(; trunc, nlayers, kwargs...)
end

function initialize!(diffusion::SpectralFilter, model::AbstractModel)
    return initialize!(diffusion, model.time_stepping)
end

function initialize!(
        diffusion::SpectralFilter,
        L::AbstractTimeStepper,
    )
    (; trunc, nlayers) = diffusion
    (; expl, impl, expl_div, impl_div) = diffusion
    (; scale, shift, power, power_div, resolution_scaling) = diffusion
    (; Δt, radius) = L

    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = Second(diffusion.time_scale).value / radius * (32 / (trunc + 1))^resolution_scaling
    time_scale_div = Second(diffusion.time_scale_div).value / radius * (32 / (trunc + 1))^resolution_scaling

    for k in 1:nlayers
        for l in 0:trunc    # diffusion for every degree l, but indendent of order m
            # Explicit part for (tend + expl*var) * impl
            expl[l + 1, k] = -(1 + tanh(scale * (l - trunc - shift)))^power / time_scale
            expl_div[l + 1, k] = -(1 + tanh(scale * (l - trunc - shift)))^power_div / time_scale_div

            # and implicit part of the diffusion
            impl[l + 1, k] = 1 / (1 - 2Δt * expl[l + 1, k])
            impl_div[l + 1, k] = 1 / (1 - 2Δt * expl_div[l + 1, k])
        end

        # last degree is only used by vector quantities; set to zero for implicit and explicit
        # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
        expl[trunc + 2, k] = 0
        impl[trunc + 2, k] = 0
        expl_div[trunc + 2, k] = 0
        impl_div[trunc + 2, k] = 0
    end

    _tile_diffusion_all!(diffusion)
    return
end
