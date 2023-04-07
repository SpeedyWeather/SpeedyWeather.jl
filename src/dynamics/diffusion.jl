"""
    horizontal_diffusion!(  tendency::LowerTriangularMatrix{Complex},
                            A::LowerTriangularMatrix{Complex},
                            ∇²ⁿ_expl::AbstractVector,
                            ∇²ⁿ_impl::AbstractVector)

Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `damp_expl` and an implicit part `damp_impl`, such that both can be calculated
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
    for m in 1:mmax         # loops over all columns/order m
        for l in m:lmax-1   # but skips the lmax+2 degree (1-based)
            lm += 1         # single index lm corresponding to harmonic l,m
            tendency[lm] = (tendency[lm] + ∇²ⁿ_expl[l]*A[lm])*∇²ⁿ_impl[l]
        end
        lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
    end
end

# which variables to apply horizontal diffusion to
diffusion_vars(::Barotropic) = (:vor,)
diffusion_vars(::ShallowWater) = (:vor,:div)
diffusion_vars(::PrimitiveEquation) = (:vor,:div,:temp)

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                model::ModelSetup,
                                lf::Int=1)      # leapfrog index used (2 is unstable)

    k = diagn.k     # pick precalculated hyperdiffusion operator for layer k
    ∇²ⁿ = model.horizontal_diffusion.∇²ⁿ[k]
    ∇²ⁿ_implicit = model.horizontal_diffusion.∇²ⁿ_implicit[k]

    for varname in diffusion_vars(model)
        var = getfield(progn.leapfrog[lf],varname)
        var_tend = getfield(diagn.tendencies,Symbol(varname,:_tend))
        horizontal_diffusion!(var_tend,var,∇²ⁿ,∇²ⁿ_implicit)
    end
end


"""
    stratospheric_zonal_drag!(  tendency::AbstractArray{Complex{NF},3}, # tendency of
                                A::AbstractArray{Complex{NF},3},        # spectral vorticity or divergence
                                drag::Real                              # drag coefficient [1/s]
                                ) where {NF<:AbstractFloat}             # number format NF

Zonal drag in the uppermost layer of the stratosphere of 3D spectral field `A` (vorticity or divergence).
Drag is applied explicitly to the time step in `A` and its tendency `tendency` is changed in-place.
`drag` is the drag coefficient of unit 1/s.
"""
function stratospheric_zonal_drag!( tendency::AbstractArray{Complex{NF},3}, # tendency of
                                    A::AbstractArray{Complex{NF},3},        # spectral vorticity or divergence
                                    drag::Real                              # drag coefficient [1/s]
                                    ) where {NF<:AbstractFloat}             # number format NF
    
    lmax,mmax,nlev = size(A)    # spherical harmonic degree l, order m, number of vertical levels nlev
    lmax -= 1                   # convert to 0-based l,m
    mmax -= 1
    @boundscheck size(A) == size(tendency) || throw(BoundsError())

    drag_NF = convert(NF,drag)

    @inbounds for l in 1:lmax+1     # loop over degree l, but 1-based
        # size(A) = lmax x mmax x nlev, nlev = 1 is uppermost model level
        # apply drag only to largest zonal wavenumber (m = 0) and in the uppermost model level (k=1)
        tendency[l,1,1] = tendency[l,1,1] - drag_NF*A[l,1,1]
    end
end

"""Orographic temperature correction for absolute temperature to be applied before the horizontal diffusion."""
function orographic_correction!(A_corrected::AbstractArray{Complex{NF},3},  # correction of 
                                A::AbstractArray{Complex{NF},3},            # 3D spectral temperature or humidity
                                correction_horizontal::Matrix{Complex{NF}}, # horizontal correction matrix
                                correction_vertical::Vector{NF},            # vertical correction vector
                                ) where NF
    
    lmax,mmax,nlev = size(A)    # degree l, order m of the spherical harmonics
    lmax -= 1                   # convert to 0-based
    mmax -= 1

    @boundscheck size(A) == size(A_corrected) || throw(BoundsError())
    @boundscheck (lmax+1,mmax+1) == size(correction_horizontal) || throw(BoundsError())
    @boundscheck (nlev,) == size(correction_vertical) || throw(BoundsError())

    @inbounds for k in 1:nlev       # vertical levels
        for m in 1:mmax+1           # order of spherical harmonics
            for l in m:lmax+1       # degree of spherical harmonics
                A_corrected[l,m,k] = A[l,m,k] + hori_correction[l,m]*vert_correction[k]
            end
        end
    end
end