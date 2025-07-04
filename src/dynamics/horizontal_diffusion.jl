abstract type AbstractHorizontalDiffusion <: AbstractModelComponent end

export HyperDiffusion

"""
Horizontal hyper diffusion of vor, div, temp, humid; implicitly in spectral space
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
} <: AbstractHorizontalDiffusion

    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical levels"
    nlayers::Int

    # PARAMETERS
    "[OPTION] power of Laplacian"
    power::NF = 4

    "[OPTION] diffusion time scale"
    time_scale::Second = Hour(4)

    "[OPTION] diffusion time scale for temperature and humidity"
    time_scale_div::Second = Hour(1)

    "[OPTION] stronger diffusion with resolution? 0: constant with trunc, 1: (inverse) linear with trunc, etc"
    resolution_scaling::NF = 1

    # incrased diffusion in stratosphere
    "[OPTION] different power for tropopause/stratosphere"
    power_stratosphere::NF = 2

    "[OPTION] linearly scale towards power_stratosphere above this σ"
    tapering_σ::NF = 0.2

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    expl::MatrixType = zeros(NF, trunc+2, nlayers)      # explicit part
    impl::MatrixType = ones(NF, trunc+2, nlayers)       # implicit part

    # ARRAYS using time_scale_div
    expl_div::MatrixType = zeros(NF, trunc+2, nlayers)    # explicit part
    impl_div::MatrixType = ones(NF, trunc+2, nlayers)     # implicit part
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function HyperDiffusion(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlayers, ArrayType) = spectral_grid        # take resolution parameters from spectral_grid
    MatrixType = ArrayType{NF, 2}
    return HyperDiffusion{NF, MatrixType}(; trunc, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms in `diffusion` based on the
model time step, and possibly with a changing strength/power in
the vertical."""
function initialize!(   diffusion::HyperDiffusion,
                        model::AbstractModel)
    initialize!(diffusion, model.geometry, model.time_stepping)
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
    ∇²ⁿ = diffusion.expl
    ∇²ⁿ_implicit = diffusion.impl
    ∇²ⁿ_div = diffusion.expl_div
    ∇²ⁿ_div_implicit = diffusion.impl_div
    (; power, power_stratosphere, tapering_σ) = diffusion
    (; Δt, radius) = L

    # Reduce diffusion time scale (=increase diffusion, always in seconds) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = Second(diffusion.time_scale).value/radius * (32/(trunc+1))^resolution_scaling
    time_scale_div = Second(diffusion.time_scale_div).value/radius * (32/(trunc+1))^resolution_scaling

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -trunc*(trunc+1)

    for k in 1:nlayers
        # VERTICAL TAPERING for the stratosphere
        # go from 1 to 0 between σ=0 and tapering_σ
        σ = G.σ_levels_full[k]
        tapering = max(0, (tapering_σ-σ)/tapering_σ)         # ∈ [0, 1]
        p = power + tapering*(power_stratosphere - power)

        for l in 0:trunc    # diffusion for every degree l, but indendent of order m
            eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

            # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
            ∇²ⁿ[l+1, k] = -eigenvalue_norm^power/time_scale
            ∇²ⁿ_div[l+1, k] = -eigenvalue_norm^p/time_scale_div

            # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
            ∇²ⁿ_implicit[l+1, k] = 1/(1-2Δt*∇²ⁿ[l+1, k])
            ∇²ⁿ_div_implicit[l+1, k] = 1/(1-2Δt*∇²ⁿ_div[l+1, k])
        end

        # last degree is only used by vector quantities; set to zero for implicit and explicit
        # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
        ∇²ⁿ[trunc+2, k] = 0
        ∇²ⁿ_implicit[trunc+2, k] = 0
        ∇²ⁿ_div[trunc+2, k] = 0
        ∇²ⁿ_div_implicit[trunc+2, k] = 0
    end
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
    lmax, mmax = size(tendency, OneBased, as=Matrix)
    nlayers = size(var, 2)

    @boundscheck size(tendency) == size(var) || throw(BoundsError)
    @boundscheck lmax <= size(expl, 1) == size(impl, 1) || throw(BoundsError)
    @boundscheck nlayers <= size(expl, 2) == size(impl, 2) || throw(BoundsError)

    @inbounds for k in eachmatrix(tendency, var)
        lm = 0                      # running index
        for m in 1:mmax             # loops over all columns/order m
            for l in m:lmax
                lm += 1             # single index lm corresponding to harmonic l, m
                tendency[lm, k] = (tendency[lm, k] + expl[l, k]*var[lm, k]) * impl[l, k]
            end
        end
    end
end

# temp. function barrier for GPU version 
function horizontal_diffusion!(
    tendency::LowerTriangularArray{NF, N, AT, <:Spectrum{<:GPU}},     # tendency of a
    var::LowerTriangularArray{NF, N, AT, <:Spectrum{<:GPU}},          # spectral horizontal field to diffuse
    expl::AbstractMatrix,               # explicit spectral damping (lmax x nlayers matrix)
    impl::AbstractMatrix,               # implicit spectral damping (lmax x nlayers matrix)
) where {NF, N, AT <: AbstractArray}
    return _horizontal_diffusion_KA!(tendency, var, expl, impl)
end

"""$(TYPEDSIGNATURES)
Kernel-based implementation of horizontal diffusion for a 2D field `var` in spectral space.
Updates the tendency `tendency` with an implicitly calculated diffusion term.
The implicit diffusion of the next time step is split into an explicit part `expl` and 
an implicit part `impl`, such that both can be calculated in a single forward step."""
function _horizontal_diffusion_KA!(
    tendency::LowerTriangularArray,     # tendency of a
    var::LowerTriangularArray,          # spectral horizontal field to diffuse
    expl::AbstractMatrix,               # explicit spectral damping (lmax x nlayers matrix)
    impl::AbstractMatrix,               # implicit spectral damping (lmax x nlayers matrix)
)
    # Launch the kernel
    launch!(architecture(var), :lmk, size(var), _horizontal_diffusion_kernel!, 
            tendency, var, expl, impl)
    synchronize(architecture(var))
    
    return tendency
end

@kernel inbounds=true function _horizontal_diffusion_kernel!(
    tendency, var, @Const(expl), @Const(impl))

    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(var) == 1 ? 1 : I[2]
    
    # Get the degree l for this coefficient
    l = var.spectrum.l_indices[lm]
    
    # Apply horizontal diffusion
    tendency[I] = (tendency[I] + expl[l, k] * var[I]) * impl[l, k]
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity in the BarotropicModel."""
function horizontal_diffusion!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    diffusion::AbstractHorizontalDiffusion,
    model::Barotropic,
    lf::Integer = 1,    # leapfrog index used (2 is unstable)
)
    (; expl, impl) = diffusion

    # Barotropic model diffuses vorticity (only variable)
    vor = get_step(progn.vor, lf)                               # lta_view for leapfrog index
    (; vor_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, expl, impl)

    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)          # lta_view for leapfrog index
        tracer_tend = diagn.tendencies.tracers_tend[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity and divergence in the ShallowWaterModel."""
function horizontal_diffusion!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    diffusion::AbstractHorizontalDiffusion,
    model::ShallowWater,
    lf::Integer = 1,    # leapfrog index used (2 is unstable)
)
    (; expl, impl, expl_div, impl_div) = diffusion

    # ShallowWater model diffuses vorticity and divergence
    vor = get_step(progn.vor, lf)
    div = get_step(progn.div, lf)
    (; vor_tend, div_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, expl, impl)
    horizontal_diffusion!(div_tend, div, expl_div, impl_div)

    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)      # lta_view for leapfrog index
        tracer_tend = diagn.tendencies.tracers_tend[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion applied to vorticity, divergence, temperature, and
humidity (PrimitiveWet only) in the PrimitiveEquation models."""
function horizontal_diffusion!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    diffusion::AbstractHorizontalDiffusion,
    model::PrimitiveEquation,
    lf::Integer = 1,    # leapfrog index used (2 is unstable)
)
    # use stronger diffusion operators that taper and change power with height for divergence
    (; expl_div, impl_div) = diffusion

    # and those for all other variables
    (; expl, impl) = diffusion

    # Primitive equation models diffuse vor, divergence, temp (and humidity for wet core)
    vor   = get_step(progn.vor, lf)
    div   = get_step(progn.div, lf)
    temp  = get_step(progn.temp, lf)
    humid = get_step(progn.humid, lf)
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, expl, impl)
    horizontal_diffusion!(div_tend, div, expl_div, impl_div)
    horizontal_diffusion!(temp_tend, temp, expl, impl)
    model isa PrimitiveWet && horizontal_diffusion!(humid_tend, humid, expl, impl)

    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)      # lta_view for leapfrog index
        tracer_tend = diagn.tendencies.tracers_tend[name]
        tracer.active && horizontal_diffusion!(tracer_tend, tracer_var, expl, impl)
    end
end

export SpectralFilter

"""Spectral filter for horizontal diffusion. Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SpectralFilter{
    NF,
    MatrixType,
} <: AbstractHorizontalDiffusion

    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical levels"
    nlayers::Int

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
    expl::MatrixType = zeros(NF, trunc+2, nlayers)  # explicit part
    impl::MatrixType = ones(NF, trunc+2, nlayers)   # implicit part

    # ARRAYS using time_scale_div for divergence
    expl_div::MatrixType = zeros(NF, trunc+2, nlayers)  # explicit part
    impl_div::MatrixType = ones(NF, trunc+2, nlayers)   # implicit part
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function SpectralFilter(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlayers, ArrayType) = spectral_grid        # take resolution parameters from spectral_grid
    MatrixType = ArrayType{NF, 2}
    return SpectralFilter{NF, MatrixType}(; trunc, nlayers, kwargs...)
end

function initialize!(diffusion::SpectralFilter, model::AbstractModel)
    initialize!(diffusion, model.time_stepping)
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
    time_scale = Second(diffusion.time_scale).value/radius * (32/(trunc+1))^resolution_scaling
    time_scale_div = Second(diffusion.time_scale_div).value/radius * (32/(trunc+1))^resolution_scaling

    for k in 1:nlayers
        for l in 0:trunc    # diffusion for every degree l, but indendent of order m
            # Explicit part for (tend + expl*var) * impl
            expl[l+1, k] =  -(1 + tanh(scale*(l - trunc - shift)))^power / time_scale
            expl_div[l+1, k] =  -(1 + tanh(scale*(l - trunc - shift)))^power_div / time_scale_div

            # and implicit part of the diffusion
            impl[l+1, k] = 1/(1-2Δt*expl[l+1, k])
            impl_div[l+1, k] = 1/(1-2Δt*expl_div[l+1, k])
        end

        # last degree is only used by vector quantities; set to zero for implicit and explicit
        # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
        expl[trunc+2, k] = 0
        impl[trunc+2, k] = 0
        expl_div[trunc+2, k] = 0
        impl_div[trunc+2, k] = 0
    end
end