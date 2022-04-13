"""
    leapfrog!(  A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                tendency::AbstractMatrix{Complex{NF}},  # tendency (dynamics+physics) of A
                dt::NF,                                 # time step (=2Δt, but for init steps =Δt,Δt/2)
                lf::Int,                                # leapfrog index to dis/enable William's filter
                C::Constants{NF}                        # struct with constants used at runtime
                ) where {NF<:AbstractFloat}             # number format NF

Performs one leapfrog time step with or without Robert+William's filter (see William (2009),
Montly Weather Review, Eq. 7-9).
"""
function leapfrog!( A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                    tendency::AbstractMatrix{Complex{NF}},  # tendency (dynamics+physics) of A
                    dt::NF,                                 # time step (=2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                                # leapfrog index to dis/enable William's filter
                    C::Constants{NF}                        # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}             # number format NF

    lmax,mmax,nleapfrog = size(A)                           # 1-based max degree l, order m of sph harmonics
                                                            # 2 leapfrog steps

    @boundscheck (lmax,mmax) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError()) # last dim is 2 for leapfrog
    @boundscheck lf in [1,2] || throw(BoundsError())    # index l1 calls leapfrog dim
    
    @unpack robert_filter, williams_filter = C          # coefficients for the Robert and William's filter
    two = convert(NF,2)                                 # 2 in number format NF

    # LEAP FROG time step with or without Robert+William's filter
    # Robert time filter to compress computational mode, Williams' filter for 3rd order accuracy
    # see William (2009), Eq. 7-9
    # for lf == 1 (initial time step) no filter applied (w1=w2=0)
    # for lf == 2 (later steps) Robert+William's filter is applied
    w1 = lf == 1 ? zero(NF) : robert_filter*williams_filter/two         # = ν*α/2 in William (2009, Eq. 8)
    w2 = lf == 1 ? zero(NF) : robert_filter*(1-williams_filter)/two     # = ν(1-α)/2 in William (2009, Eq. 9)

    @inbounds for m in 1:mmax                       # 1-based max degree l, order m 
        for l in m:lmax                             # skip upper right triangle as they should be zero anyway
            Aold = A[l,m,1]                         # double filtered value from previous time step (t-Δt)
            Anew = Aold + dt*tendency[l,m]          # Leapfrog/Euler step depending on dt=Δt,2Δt (unfiltered at t+Δt)
            Aupdate = Aold - two*A[l,m,lf] + Anew   # Eq. 8&9 in William (2009), calculate only once
            A[l,m,1] = A[l,m,lf] + w1*Aupdate       # Robert's filter: A[l,m,1] becomes 2xfiltered value at t
            A[l,m,2] = Anew - w2*Aupdate            # Williams' filter: A[l,m,2] becomes 1xfiltred value at t+Δt
        end
    end
end

"""Leapfrog! for 3D arrays that loops over all vertical layers."""
function leapfrog!( A::AbstractArray{Complex{NF},4},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},3}, # tendency (dynamics + physics) of A
                    dt::NF,                                 # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                                # leapfrog index to dis/enable William's filter
                    C::Constants{NF}                        # struct containing all constants used at runtime
                    ) where {NF<:AbstractFloat}             # number format NF

    _,_,_,nlev = size(A)        # A is of size lmax x mmax x nlev x 2

    for k in 1:nlev
        A_layer = view(A,:,:,k,:)                   # extract vertical layers as views to not allocate any memory
        tendency_layer = view(tendency,:,:,k)
        leapfrog!(A_layer,tendency_layer,dt,lf,C)   # make a timestep forward for each layer
    end
end

"""TODO write a leapfrog! function that loops over all prognostic variables."""
function leapfrog!( Prog::PrognosticVariables,
                    Tend::Tendencies)
    
    # @distributed 
    for (var,tend) in zip((pres_surf,),(pres_surf_tend,))
        leapfrog!(var,tend,)
    end 
end    

"""Call initialization of semi-implicit scheme and perform initial time step."""
function first_timestep!(   Prog::PrognosticVariables{NF},  # all prognostic variables
                            Diag::PrognosticVariables{NF},  # all pre-allocated diagnostic variables
                            C::Constants{NF},               # struct containing constants
                            G::GeoSpectral{NF},             # struct containing geometry and spectral transform constants
                            HD::HorizontalDiffusion{NF}     # struct containing horizontal diffusion constants
                            ) where {NF<:AbstractFloat}
    
    @unpack Δt = C

    # FIRST TIME STEP (l1=l2=1, dt=Δt/2)
    IMP = initialize_implicit(half*Δt)
    timestep!(Prog,Diag,1,1,Δt/2,C,G,HD)

    # SECOND TIME STEP (l1=1,l2=2,dt=Δt)
    IMP = initialize_implicit(Δt)
    timestep!(Prog,Diag,1,2,Δt,C,G,HD)

    # Initialize implicit arrays for further time steps (dt=2Δt)
    IMP = initialize_implicit(2Δt)
    return IMP
end

"""Calculate a single time step for SpeedyWeather.jl"""
function timestep!( Prog::PrognosticVariables{NF},  # all prognostic variables
                    Diag::PrognosticVariables{NF},  # all pre-allocated diagnostic variables
                    l1::Int,                        # leapfrog index 1 (en/disables Robert+William's filter)
                    l2::Int,                        # leapfrog index 2 (time step used for tendencies)
                    dt::NF,                         # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    C::Constants{NF},               # struct containing constants
                    G::GeoSpectral{NF},             # struct containing geometry and spectral transform constants
                    HD::HorizontalDiffusion{NF}     # struct containing horizontal diffusion constants
                    ) where {NF<:AbstractFloat}
                        
    @unpack vor,div,ps,temp = Prog
    @unpack vor_tend,div_tend,ps_tend,temp_tend = Diag.Tendencies
    @unpack dmp,dmpd,dmps,dmp1,dmp1d,dmp1s = HD
    @unpack tcorh,tcorv,qcorh,qcorv = HD
    @unpack sdrag = C

    # COMPUTE TENDENCIES OF PROGNOSTIC VARIABLES
    get_tendencies!(Prog,Diag,l2,C)                   

    # DIFFUSION FOR WIND
    vor_l1 = view(vor,:,:,1,:)                                  # array view for leapfrog index
    div_l1 = view(div,:,:,1,:)                                  # TODO l1/l2 dependent?
    horizontal_diffusion!(vor_l1,vor_tend,dmp,dmp1)             # diffusion of vorticity
    horizontal_diffusion!(div_l1,div_tend,dmpd,dmp1d)           # diffusion of divergence

    # DIFFUSION FOR TEMPERATURE
    orographic_correction!(temp_corrected,temp,1,tcorh,tcorv)   # orographic correction for temperature
    horizontal_diffusion!(temp_corrected,temp_tend,dmp,dmp1)    # diffusion for corrected abs temperature

    # DISSIPATION IN THE STRATOSPHERE
    stratospheric_zonal_drag!(vor,vor_tend,sdrag)               # zonal drag for wind
    stratospheric_zonal_drag!(div,div_tend,sdrag)

    horizontal_diffusion!(vor_l1,vor_tend,dmps,dmp1s)           # stratospheric diffusion for wind
    horizontal_diffusion!(div_l1,div_tend,dmps,dmp1s)           
    horizontal_diffusion!(temp_corrected,temp_tend,dmps,dmp1s)  # stratospheric diffusion for temperature

    # DIFFUSION OF HUMIDITY
    orographic_correction!(humid_corrected,humid,1,qcorh,qcorv) # orographic correction for humidity
    horizontal_diffusion!(humid_corrected,humid_tend,dmp,dmp1)  # diffusion for corrected humidity

    # if ntracers > 1
    #     for i in 2:ntracers
    #         tracer          = view(tracers,:,:,:,1,i)                   # the i-th tracer, leapfrog index 1
    #         tracer_tendency = view(tracers_tendencies,:,:,:,i)          # its tendency
    #         horizontal_diffusion!(tracer,tracer_tendency,dmp,dmp1)
    #         leapfrog!(tracer,tracer_tendency,j1,C)
    #     end
    # end

    # SPECTRAL TRUNCATION of all tendencies to the spectral resolution
    for tendency in (logp0_tend,vor_tend,div_tend,temp_tend,humid_tend)
        spectral_truncation!(tendency,G)
    end

    # Time integration via leapfrog step forward (filtered with Robert+William's depending on l1)
    leapfrog!(logp0,logp0_tend,dt,l1,C)
    leapfrog!(vor,vor_tend,dt,l1,C)
    leapfrog!(div,div_tend,dt,l1,C)
    leapfrog!(temp,temp_tend,dt,l1,C)
    leapfrog!(humid,humid_tend,dt,l1,C)
end

"""Calculate a single time step for SpeedyWeather.jl"""
function time_stepping!(prog::PrognosticVariables{NF},  # all prognostic variables
                        diag::DiagnosticVariables{NF},  # all pre-allocated diagnostic variables
                        M::ModelSetup{NF}               # all precalculated structs
                        ) where {NF<:AbstractFloat}     # number format NF
    
    @unpack n_timesteps, Δt = M.constants
    @unpack output = M.parameters

    # FEEDBACK, OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS
    feedback = initialize_feedback(M)
    # netcdf_files = initialize_output(diag,feedback,M)

    # first_timestep!(prog,diag,C,G,HD)

    for i in 1:n_timesteps
        # timestep!(prog,diag,2,2,2Δt,C,G,HD)

        # FEEDBACK AND OUTPUT
        feedback!(feedback,i)
        # output_netcdf!(i,netcdf_files,diag,M)
    end

    feedback_end!(feedback)

    return prog
end
    