"""
    timestep!()

Time step of a 2D field `input` (a `Lon x Lat x 2` array)

Perform one time step starting from F(1) and F(2) and using the following scheme:
Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
F(2) = Fnew
Input:
If j1 == 1, j2 == 1 : forward time step (eps = 0)
If j1 == 1, j2 == 2 : initial leapfrog time step (eps = 0)
If j1 == 2, j2 == 2 : leapfrog time step with time filter (eps = ROB)
dt = time step

"""
function leapfrog!( A::AbstractArray{NF,3},             # a prognostic variable
                    tendency::AbstractArray{NF,2},      # its tendency
                    l1::Int,                            # leapfrog index for time filtering
                    G::GeoSpectral{NF},                 # struct with precomputed arrays for spectral transform
                    C::Constants{NF}                    # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}

    nlon,nlat,nleapfrog = size(A)                       # longitude, latitude, 2 leapfrog steps

    @boundscheck (nlon,nlat) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError())     # last dim is 2 for leapfrog
    @boundscheck l1 in [1,2] || throw(BoundsError())        # index l1 is calls leapfrog dim
    
    # get coefficients for the Robert and Williams' filter for 3rd order accuracy in time stepping
    @unpack robert_filter, williams_filter, Δt = C
    eps = l1 == 1 ? zero(NF) : robert_filter
    two = convert(NF,2)
    eps = one(NF) - two*eps

    # truncate the tendency to the spectral resolution
    spectral_truncation!(tendency,G)

    # LEAP FROG time step with
    w1 = williams_filter*eps                # Robert time filter to compress computational mode
    w2 = (one(NF) - williams_filter)*eps    # and Williams' filter for 3rd order accuracy
    @inbounds for j in 1:nlat
        for i in 1:nlon
            Anew = A[i,j,1] + Δt*tendency[i,j]                          # forward step
            A[i,j,1] = A[i,j,l1] + w1*(A[i,j,1] - two*A[i,j,l1] + Anew) # Robert's filter
            A[i,j,2] = Anew - w2*(A[i,j,1] - two*A[i,j,l1] + Anew)      # Williams' filter
        end
    end
end

"""3D version that loops over all vertical layers."""
function leapfrog!( A::AbstractArray{NF,4},             # a prognostic variable
                    tendency::AbstractArray{NF,3},      # its tendency
                    l1::Int,                            # index for time filtering
                    C::Constants{NF}                    # struct containing all constants used at runtime
                    ) where {NF<:AbstractFloat}

    _,_,nlev,_ = size(A)        # A is of size nlon x nlat x nlev x 2

    for k in 1:nlev
        # extract vertical layers as views to not allocate any memory
        A_layer = view(A,:,:,k,:)
        tendency_layer = view(tendency,:,:,k)
        
        # make a timestep forward for each layer
        leapfrog!(A_layer,tendency_layer,l1,C)
    end
end

# Call initialization of semi-implicit scheme and perform initial time step
function first_step()
    initialize_implicit(half*Δt)

    step(1, 1, half*Δt)

    initialize_implicit(Δt)

    step(1, 2, Δt)

    initialize_implicit(2*Δt)
end

function timestep!( Prog::PrognosticVariables{NF},
                    Diag::PrognosticVariables{NF},
                    C::Constants{NF},
                    G::GeoSpectral{NF},
                    HD::HorizontalDiffusion{NF}
                    ) where {NF<:AbstractFloat}
                        
    @unpack vor,div,ps,temp = Prog
    @unpack vor_tend,div_tend,ps_tend,temp_tend = Diag.Tendencies
    @unpack dmp,dmpd,dmps,dmp1,dmp1d,dmp1s = HD
    @unpack sdrag = C

    # Compute tendencies of prognostic variables
    # =========================================================================

    get_tendencies!(vorU_tend, divU_tend, tem_tend, pₛ_tend, tr_tend, j2)

    # =========================================================================
    # Horizontal diffusion
    # =========================================================================

    # Diffusion of wind and temperature
    horizontal_diffusion!(vorU[:,:,:,1], vorU_tend, dmp,  dmp1)
    horizontal_diffusion!(divU[:,:,:,1], divU_tend, dmpd, dmp1d)

    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = tem[m,n,k,1] + tcorh[m,n]*tcorv[k]
            end
        end
    end

    horizontal_diffusion!(ctmp, tem_tend, dmp, dmp1)

    stratospheric_zonal_drag!(vor,vor_tend,sdrag)
    stratospheric_zonal_drag!(div,div_tend,sdrag)

    horizontal_diffusion!(vorU[:,:,:,1],  vorU_tend, dmps, dmp1s)
    horizontal_diffusion!(divU[:,:,:,1],  divU_tend, dmps, dmp1s)
    horizontal_diffusion!(ctmp, tem_tend,   dmps, dmp1s)

    # Diffusion of tracers
    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = tr[m,n,k,1,1] + qcorh[m,n]*qcorv[k]
            end
        end
    end

    horizontal_diffusion!(ctmp, @view(tr_tend[:,:,:,1]), dmpd, dmp1d)

    if ntracers > 1
        for i in 2:ntracers
            tracer          = view(tracers,:,:,:,1,i)                   # the i-th tracer, leapfrog index 1
            tracer_tendency = view(tracers_tendencies,:,:,:,i)          # its tendency
            horizontal_diffusion!(tracer,tracer_tendency,dmp,dmp1)
            leapfrog!(tracer,tracer_tendency,j1,C)
        end
    end

    # Time integration
    leapfrog!(ps,ps_tend,j1,G,C)
    leapfrog!(vor,vor_tend,j1,G,C)
    leapfrog!(div,div_tend,j1,G,C)
    leapfrog!(temp,temp_tend,j1,G,C)

    #TODO swap leapfrog indices here?
end

