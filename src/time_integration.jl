"""
Perform one leapfrog time step with or without Robert+William's filter (see William (2009),
Montly Weather Review, Eq. 7-9)
"""
function leapfrog!( A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},2}, # its tendency (dynamics+physics)
                    dt::NF,                                 # time step (=2Δt, but for init steps =Δt,Δt/2)
                    l1::Int,                                # leapfrog index to dis/enable William's filter
                    C::Constants{NF}                        # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}             # number format NF

    mx,nx,nleapfrog = size(A)                           # spectral mx x nx, 2 leapfrog steps

    @boundscheck (mx,nx) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError()) # last dim is 2 for leapfrog
    @boundscheck l1 in [1,2] || throw(BoundsError())    # index l1 calls leapfrog dim
    
    @unpack robert_filter, williams_filter = C          # coefficients for the Robert and William's filter
    two = convert(NF,2)

    # LEAP FROG time step with or without Robert+William's filter
    # Robert time filter to compress computational mode, Williams' filter for 3rd order accuracy
    # see William (2009), Eq. 7-9
    # for l1 == 1 (initial time step) no filter applied (w1=w2=0)
    # for l1 == 2 (later steps) Robert+William's filter is applied
    w1 = l1 == 1 ? zero(NF) : robert_filter*williams_filter/two         # = ν*α/2 in William (2009, Eq. 8)
    w2 = l1 == 1 ? zero(NF) : robert_filter*(1-williams_filter)/two     # = ν(1-α)/2 in William (2009, Eq. 9)

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
    get_tendencies!(Prog,Diag,l2)                   

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
function time_stepping!(Prog::PrognosticVariables{NF},  # all prognostic variables
                        Diag::PrognosticVariables{NF},  # all pre-allocated diagnostic variables
                        C::Constants{NF},               # struct containing constants
                        G::GeoSpectral{NF},             # struct containing geometry and spectral transform constants
                        HD::HorizontalDiffusion{NF}     # struct containing horizontal diffusion constants
                        P::Params                       # struct containing all model parameters
                        ) where {NF<:AbstractFloat}
    
    @unpack n_timesteps, Δt = C
    @unpack output = P

    # FEEDBACK, OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS
    # feedback = feedback_initialise(S)
    # ncfile = output_initialise(feedback,S)

    first_timestep!(Prog,Diag,C,G,HD)

    for i in 1:n_timesteps
        timestep!(Prog,Diag,2,2,2Δt,C,G,HD)

        # FEEDBACK AND OUTPUT
        # feedback.i = i
        # feedback!(Prog,feedback,S)
        # output_nc!(i,netCDFfiles,Prog,Diag,S)
    end

    return Prog
end
    