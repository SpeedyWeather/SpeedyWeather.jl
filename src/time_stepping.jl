"""
Perform one leapfrog time step with Robert's or Robert+William's filter (see William (2009),
Montly Weather Review, Eq. 7-9)
"""
function leapfrog!( A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},2}, # its tendency (dynamics+physics)
                    dt::NF,                                 # time step (=2Δt, but for init steps =Δt,Δt/2)
                    l1::Int,                                # leapfrog index to dis/enable William's filter
                    C::Constants{NF}                        # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}             # number format NF

    mx,nx,nleapfrog = size(A)                               # spectral mx x nx, 2 leapfrog steps

    @boundscheck (mx,nx) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError())     # last dim is 2 for leapfrog
    @boundscheck l1 in [1,2] || throw(BoundsError())        # index l1 calls leapfrog dim
    
    # get coefficients for the Robert and Williams' filter for 3rd order accuracy in time stepping
    @unpack robert_filter, williams_filter = C
    two = convert(NF,2)

    # LEAP FROG time step with Robert or Robert+William's filter
    # see William (2009), Eq. 7-9
    # for l1 == 1 (initial time step) William's is disabled (α=1)
    # for l1 == 2 (later steps) William's is included (α=0.53 by default)
    α = l1 == 1 ? one(NF) : williams_filter # α=1 means no William's filter (=> w2=0)
    w1 = robert_filter*α/two                # Robert time filter to compress computational mode
    w2 = robert_filter*(one(NF) - α)/two    # and Williams' filter for 3rd order accuracy                     

    @inbounds for j in 1:nx
        for i in 1:mx
            Aold = A[i,j,1]                         # double filtered value from previous time step (t-Δt)
            Anew = Aold + dt*tendency[i,j]          # Leapfrog/Euler step depending on dt=Δt,2Δt (unfiltered at t+Δt)
            Aupdate = Aold - two*A[i,j,l1] + Anew   # Eq. 8&9 in William (2009), calculate only once
            A[i,j,1] = A[i,j,l1] + w1*Aupdate       # Robert's filter: A[i,j,1] becomes the double filtered value at t
            A[i,j,2] = Anew - w2*Aupdate            # Williams' filter: A[i,j,2] becomes the single filtred value at t+Δt
        end
    end
end

"""Leapfrog! for 3D arrays that loops over all vertical layers."""
function leapfrog!( A::AbstractArray{Complex{NF},4},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},3}, # its tendency (dynamics + physics)
                    dt::NF,                                 # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    l1::Int,                                # leapfrog index to dis/enable William's filter
                    C::Constants{NF}                        # struct containing all constants used at runtime
                    ) where {NF<:AbstractFloat}             # number format NF

    _,_,_,nlev = size(A)        # A is of size mx x nx x 2 x nlev

    for k in 1:nlev
        # extract vertical layers as views to not allocate any memory
        A_layer = view(A,:,:,:,k)
        tendency_layer = view(tendency,:,:,k)
        
        # make a timestep forward for each layer
        leapfrog!(A_layer,tendency_layer,dt,l1,C)
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

    # TODO apply to all tendencies
    # truncate the tendency to the spectral resolution
    spectral_truncation!(tendency,G)

    # Time integration
    leapfrog!(ps,ps_tend,j1,G,C)
    leapfrog!(vor,vor_tend,j1,G,C)
    leapfrog!(div,div_tend,j1,G,C)
    leapfrog!(temp,temp_tend,j1,G,C)

    #TODO swap leapfrog indices here?
end

