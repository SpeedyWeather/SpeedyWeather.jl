Base.@kwdef struct HyperDiffusion{NF} <: HorizontalDiffusion{NF}
    # DIMENSIONS
    trunc::Int                          # spectral resolution
    nlev::Int                           # number of vertical levels
    
    # PARAMETERS
    power::Float64 = 4.0                # Power of Laplacian
    time_scale::Float64 = 2.4           # Diffusion time scales [hrs]
    resolution_scaling::Float64 = 0.5   # 0: constant with trunc
                                        # 1: (inverse) linear with trunc
                                        # 2: (inverse) quadratic, etc

    # incrased diffusion in stratosphere
    power_stratosphere::Float64 = 2.0   # different power for stratosphere
    tapering_σ::Float64 = 0.2           # scale towards that power linearly above this σ

    # increase diffusion based on high vorticity levels
    adaptive::Bool = true               # swith on/off
    vor_max::Float64 = 1e-4             # [1/s] above this, diffusion is increased
    adaptive_strength::Float64 = 2.0    # increase strength above vor_max by this factor
                                        # times max(abs(vor))/vor_max

    # constant arrays to be initalised later
    # (Hyper) diffusion, precalculated for each spherical harm degree
    # and for each layer (to allow for varying orders/strength in the vertical)
    ∇²ⁿ::Vector{Vector{NF}} = [zeros(NF,trunc+2) for _ in 1:nlev]           # explicit part
    ∇²ⁿ_implicit::Vector{Vector{NF}} = [zeros(NF,trunc+2) for _ in 1:nlev]  # implicit part
end

function HyperDiffusion(spectral_grid::SpectralGrid,kwargs...)
    (;NF,trunc,nlev) = spectral_grid        # take resolution parameters from spectral_grid
    return HyperDiffusion{NF}(;trunc,nlev,kwargs...)
end

function initialize!(   scheme::HyperDiffusion,
                                P::Parameters,
                                C::DynamicsConstants,
                                G::Geometry,
                                S::SpectralTransform{NF}) where NF
    (;lmax) = S
    (;radius) = P.planet
    (;time_scale, resolution_scaling) = scheme
    (;nlev) = G

    # Reduce diffusion time scale (=increase diffusion) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = 1/radius*3600*scheme.time_scale * (32/(lmax+1))^resolution_scaling

    # PREALLOCATE as vector as only dependend on degree l
    # Damping coefficients for explicit part of the diffusion (=ν∇²ⁿ)
    ∇²ⁿ = [zeros(NF,lmax+2) for _ in 1:nlev]                # for temp, vor, div (explicit)
    ∇²ⁿ_implicit = [zeros(NF,lmax+2) for _ in 1:nlev]       # Implicit part (= 1/(1+2Δtν∇²ⁿ))
    HD = HorizontalDiffusion{NF}(lmax,time_scale,∇²ⁿ,∇²ⁿ_implicit)  # initialize struct

    # PRECOMPUTE diffusion operators based on non-adaptive diffusion
    vor_max = 0                                 # not increased diffusion initially
    adapt_diffusion!(HD,scheme,vor_max,G,C)     # fill diffusion arrays
    return HD
end

function adapt_diffusion!(  HD::HorizontalDiffusion,
                            scheme::HyperDiffusion,
                            vor_max::Real,
                            G::Geometry,
                            C::DynamicsConstants)
    nlev = length(HD.∇²ⁿ)
    for k in 1:nlev
        adapt_diffusion!(HD,scheme,vor_max,k,G,C)
    end
end

function adapt_diffusion!(  HD::HorizontalDiffusion,
                            scheme::HyperDiffusion,
                            vor_max::Real,
                            k::Int,
                            G::Geometry,
                            C::DynamicsConstants)

    (;lmax,time_scale,∇²ⁿ,∇²ⁿ_implicit) = HD
    (;power, power_stratosphere, tapering_σ) = scheme
    (;Δt) = C
    σ = G.σ_levels_full[k]

    # ADAPTIVE/FLOW AWARE
    # increase diffusion if maximum vorticity per layer is larger than scheme.vor_max
    if scheme.adaptive
        # /= as 1/time_scale*∇²ʰ below
        time_scale /= 1 + (scheme.adaptive_strength-1)*max(0,vor_max/scheme.vor_max - 1)
    end

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -lmax*(lmax+1)
    
    # VERTICAL TAPERING for the stratosphere
    # go from 1 to 0 between σ=0 and tapering_σ
    tapering = max(0,(tapering_σ-σ)/tapering_σ)         # ∈ [0,1]
    p = power + tapering*(power_stratosphere - power) 
        
    @inbounds for l in 0:lmax+1   # diffusion for every degree l, but indendent of order m
        eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

        # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
        ∇²ⁿ[k][l+1] = -eigenvalue_norm^p/time_scale
        
        # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
        ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
    end
end

"""
    horizontal_diffusion!(  tendency::LowerTriangularMatrix{Complex},
                            A::LowerTriangularMatrix{Complex},
                            ∇²ⁿ_expl::AbstractVector,
                            ∇²ⁿ_impl::AbstractVector)

Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `∇²ⁿ_expl` and an implicit part `∇²ⁿ_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( tendency::LowerTriangularMatrix{Complex{NF}},   # tendency of a 
                                A::LowerTriangularMatrix{Complex{NF}},          # spectral horizontal field
                                ∇²ⁿ_expl::AbstractVector{NF},                   # explicit spectral damping
                                ∇²ⁿ_impl::AbstractVector{NF}                    # implicit spectral damping
                                ) where {NF<:AbstractFloat}
    lmax,mmax = size(tendency)      # 1-based
    @boundscheck size(tendency) == size(A) || throw(BoundsError)
    @boundscheck lmax <= length(∇²ⁿ_expl) == length(∇²ⁿ_impl) || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax   # loops over all columns/order m
        for l in m:lmax-1       # but skips the lmax+2 degree (1-based)
            lm += 1             # single index lm corresponding to harmonic l,m
            tendency[lm] = (tendency[lm] + ∇²ⁿ_expl[l]*A[lm])*∇²ⁿ_impl[l]
        end
        lm += 1             # skip last row for scalar quantities
    end
end

# which variables to apply horizontal diffusion to
diffusion_vars(::Barotropic) = (:vor,)
diffusion_vars(::ShallowWater) = (:vor,:div)
diffusion_vars(::PrimitiveEquation) = (:vor,:div,:temp)

function horizontal_diffusion!( progn::PrognosticLayerTimesteps,
                                diagn::DiagnosticVariablesLayer,
                                model::ModelSetup,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    
    HD = model.horizontal_diffusion
    initialize!(HD,diagn,model)
    k = diagn.k                                 # current layer k
    ∇²ⁿ = HD.∇²ⁿ[k]                             # now pick operators at k
    ∇²ⁿ_implicit = HD.∇²ⁿ_implicit[k]

    for varname in diffusion_vars(model)
        var = getfield(progn.timesteps[lf],varname)
        var_tend = getfield(diagn.tendencies,Symbol(varname,:_tend))
        horizontal_diffusion!(var_tend,var,∇²ⁿ,∇²ⁿ_implicit)
    end
end

function initialize!(   scheme::HyperDiffusion,
                        diagn::DiagnosticVariablesLayer,
                        model::ModelSetup)
    scheme.adaptive || return nothing
    vor_min, vor_max = extrema(diagn.grid_variables.vor_grid)
    vor_abs_max = max(abs(vor_min), abs(vor_max))/model.geometry.radius
    initialize!(scheme,vor_abs_max,diagn.k,model)
end