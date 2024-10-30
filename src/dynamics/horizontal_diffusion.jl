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
    time_scale::Second = Minute(60)

    "[OPTION] diffusion time scale for temperature and humidity"
    time_scale_temp_humid::Second = Minute(144)
    
    "[OPTION] stronger diffusion with resolution? 0: constant with trunc, 1: (inverse) linear with trunc, etc"
    resolution_scaling::NF = 0.5

    # incrased diffusion in stratosphere
    "[OPTION] different power for tropopause/stratosphere"
    power_stratosphere::NF = 2
    
    "[OPTION] linearly scale towards power_stratosphere above this σ"
    tapering_σ::NF = 0.2

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    ∇²ⁿ::MatrixType = zeros(NF, trunc+2, nlayers)           # explicit part
    ∇²ⁿ_implicit::MatrixType = ones(NF, trunc+2, nlayers)   # implicit part

    # ARRAYS but no scaling or tapering and using time_scale_temp_humid
    ∇²ⁿc::MatrixType = zeros(NF, trunc+2, nlayers)          # explicit part
    ∇²ⁿc_implicit::MatrixType = ones(NF, trunc+2, nlayers)  # implicit part
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
    (; ∇²ⁿ, ∇²ⁿ_implicit, ∇²ⁿc, ∇²ⁿc_implicit) = diffusion
    (; power, power_stratosphere, tapering_σ) = diffusion
    (; Δt, radius) = L

    # Reduce diffusion time scale (=increase diffusion, always in seconds) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = Second(diffusion.time_scale).value/radius * (32/(trunc+1))^resolution_scaling
    time_scale_constant = Second(diffusion.time_scale_temp_humid).value/radius

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
            ∇²ⁿ[l+1, k] = -eigenvalue_norm^p/time_scale
            ∇²ⁿc[l+1, k] = -eigenvalue_norm^power/time_scale_constant
            
            # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
            ∇²ⁿ_implicit[l+1, k] = 1/(1-2Δt*∇²ⁿ[l+1, k])           
            ∇²ⁿc_implicit[l+1, k] = 1/(1-2Δt*∇²ⁿc[l+1, k])           
        end
        
        # last degree is only used by vector quantities; set to zero for implicit and explicit
        # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
        ∇²ⁿ[trunc+2, k] = 0
        ∇²ⁿ_implicit[trunc+2, k] = 0
        ∇²ⁿc[trunc+2, k] = 0
        ∇²ⁿc_implicit[trunc+2, k] = 0
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `∇²ⁿ_expl` and an implicit part `∇²ⁿ_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( 
    tendency::LowerTriangularArray,     # tendency of a 
    var::LowerTriangularArray,          # spectral horizontal field to diffuse
    ∇²ⁿ_expl::AbstractMatrix,           # explicit spectral damping (lmax x nlayers matrix)
    ∇²ⁿ_impl::AbstractMatrix,           # implicit spectral damping (lmax x nlayers matrix)
)
    lmax, mmax = size(tendency, OneBased, as=Matrix)
    nlayers = size(var, 2)

    @boundscheck size(tendency) == size(var) || throw(BoundsError)
    @boundscheck lmax <= size(∇²ⁿ_expl, 1) == size(∇²ⁿ_impl, 1) || throw(BoundsError)
    @boundscheck nlayers <= size(∇²ⁿ_expl, 2) == size(∇²ⁿ_impl, 2) || throw(BoundsError)

    @inbounds for k in eachmatrix(tendency, var)
        lm = 0                      # running index
        for m in 1:mmax             # loops over all columns/order m
            for l in m:lmax
                lm += 1             # single index lm corresponding to harmonic l, m
                tendency[lm, k] = (tendency[lm, k] + ∇²ⁿ_expl[l, k]*var[lm, k]) * ∇²ⁿ_impl[l, k]
            end
        end
    end
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
    (; ∇²ⁿ, ∇²ⁿ_implicit) = model.horizontal_diffusion

    # Barotropic model diffuses vorticity (only variable)
    vor = progn.vor[lf]
    (; vor_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
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
    (; ∇²ⁿ, ∇²ⁿ_implicit) = model.horizontal_diffusion

    # ShallowWater model diffuses vorticity and divergence
    vor = progn.vor[lf]
    div = progn.div[lf]
    (; vor_tend, div_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion applied to vorticity, divergence, temperature, and
humidity (PrimitiveWet only) in the PrimitiveEquation models."""
function horizontal_diffusion!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    diffusion::HyperDiffusion,
    model::PrimitiveEquation,
    lf::Integer = 1,    # leapfrog index used (2 is unstable)
)
    # use weaker diffusion operators that don't taper or scale with resolution for temperature and humidity
    (; ∇²ⁿc, ∇²ⁿc_implicit) = model.horizontal_diffusion

    # and the ones that do for vorticity and divergence
    (; ∇²ⁿ, ∇²ⁿ_implicit) = model.horizontal_diffusion

    # Primitive equation models diffuse vor, divergence, temp (and humidity for wet core)
    vor = progn.vor[lf]
    div = progn.div[lf]
    temp = progn.temp[lf]
    humid = progn.humid[lf]
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(temp_tend, temp, ∇²ⁿc, ∇²ⁿc_implicit)
    model isa PrimitiveWet && horizontal_diffusion!(humid_tend, humid, ∇²ⁿc, ∇²ⁿc_implicit)
end

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
    "[OPTION] relative wavenumber"
    wavenumber::NF = 0

    "[OPTION] Scale-selectiveness"
    scale::NF = 0.05
    
    "[OPTION] diffusion time scale"
    time_scale::Second = Minute(144)

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    ∇²ⁿ::MatrixType = zeros(NF, trunc+2, nlayers)           # explicit part
    ∇²ⁿ_implicit::MatrixType = ones(NF, trunc+2, nlayers)   # implicit part
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function SpectralFilter(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlayers, ArrayType) = spectral_grid        # take resolution parameters from spectral_grid
    MatrixType = ArrayType{NF, 2}
    return SpectralFilter{NF, MatrixType}(; trunc, nlayers, kwargs...)
end