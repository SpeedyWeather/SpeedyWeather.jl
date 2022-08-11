"""
    horizontal_diffusion!(  tendency::AbstractMatrix{Complex{NF}}, # tendency of a 
                            A::AbstractMatrix{Complex{NF}},        # spectral horizontal field
                            damp_expl::AbstractMatrix{NF},         # explicit spectral damping
                            damp_impl::AbstractMatrix{NF}          # implicit spectral damping
                            ) where {NF<:AbstractFloat}

Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `damp_expl` and an implicit part `damp_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( tendency::LowerTriangularMatrix{Complex{NF}},   # tendency of a 
                                A::LowerTriangularMatrix{Complex{NF}},          # spectral horizontal field
                                damp_expl::LowerTriangularMatrix{NF},           # explicit spectral damping
                                damp_impl::LowerTriangularMatrix{NF}            # implicit spectral damping
                                ) where {NF<:AbstractFloat}
    
    for lm in eachharmonic(tendency,A,damp_expl,damp_impl)
        @inbounds tendency[lm] = (tendency[lm] - damp_expl[lm]*A[lm])*damp_impl[lm]
    end
end

@kernel function horizontal_diffusion_kernel!(tendency, @Const(A), @Const(damp_expl), @Const(damp_impl))
    lm = @index(Global, Linear)
   
    tendency[lm] = (tendency[lm] - damp_expl[lm]*A[lm])*damp_impl[lm]
end

function horizontal_diffusion!( tendency::LowerTriangularMatrix{Complex{NF}},   # tendency of a 
                                A::LowerTriangularMatrix{Complex{NF}},          # spectral horizontal field
                                damp_expl::LowerTriangularMatrix{NF},           # explicit spectral damping
                                damp_impl::LowerTriangularMatrix{NF},           # implicit spectral damping
                                device_setup::DeviceSetup,        # device the function is executed on, we have to add that
                                ) where {NF<:AbstractFloat}

    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    @boundscheck size(A) == size(damp_expl) || throw(BoundsError())
    @boundscheck size(A) == size(damp_impl) || throw(BoundsError())

    ev = launch_kernel!(device_setup, horizontal_diffusion_kernel!, length(tendency), tendency, A, damp_expl, damp_impl)
    wait(ev)
end 

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::BarotropicModel,
                                lf::Int=1,                          # leapfrog index used (2 is unstable)
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor = progn.leapfrog[lf]
    @unpack vor_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend,vor,damping,damping_impl)    # diffusion of vorticity
end

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::ShallowWaterModel,
                                lf::Int=1,                          # leapfrog index used (2 is unstable)
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor, div = progn.leapfrog[lf]
    @unpack vor_tend,div_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend,vor,damping,damping_impl)        # diffusion of vorticity
    horizontal_diffusion!(div_tend,div,damping,damping_impl)        # diffusion of divergence
end

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::PrimitiveEquationModel,
                                lf::Int=1,                          # leapfrog index used (2 is unstable)
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor,div,temp = progn.leapfrog[lf]
    @unpack vor_tend,div_tend,temp_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend, vor, damping,damping_impl)      # diffusion of vorticity
    horizontal_diffusion!(div_tend, div, damping,damping_impl)      # diffusion of divergence
    horizontal_diffusion!(temp_tend,temp,damping,damping_impl)      # diffusion of temperature
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