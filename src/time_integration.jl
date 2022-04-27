"""
    leapfrog!(  A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                tendency::AbstractMatrix{Complex{NF}},  # tendency (dynamics+physics) of A
                dt::NF,                                 # time step (=2Δt, but for init steps =Δt,Δt/2)
                C::Constants{NF},                       # struct with constants used at runtime
                lf::Int=2                               # leapfrog index to dis/enable(default) William's filter
                ) where {NF<:AbstractFloat}             # number format NF

Performs one leapfrog time step with or without Robert+William's filter (see William (2009),
Montly Weather Review, Eq. 7-9).
"""
function leapfrog!( A::AbstractArray{Complex{NF},3},        # a prognostic variable (spectral)
                    tendency::AbstractMatrix{Complex{NF}},  # tendency (dynamics+physics) of A
                    dt::Real,                               # time step (=2Δt, but for init steps =Δt,Δt/2)
                    C::Constants{NF},                       # struct with constants used at runtime
                    lf::Int=2,                              # leapfrog index to dis/enable William's filter
                    ) where {NF<:AbstractFloat}             # number format NF

    lmax,mmax,nleapfrog = size(A)                           # 1-based max degree l, order m of sph harmonics
                                                            # 2 leapfrog steps

    @boundscheck (lmax,mmax) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError()) # last dim is 2 for leapfrog
    @boundscheck lf in [1,2] || throw(BoundsError())    # index l1 calls leapfrog dim
    
    @unpack robert_filter, williams_filter = C          # coefficients for the Robert and William's filter
    two = convert(NF,2)                                 # 2 in number format NF
    dt_NF = convert(NF,dt)                              # time step dt in number format NF

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
            Anew = Aold + dt_NF*tendency[l,m]       # Leapfrog/Euler step depending on dt=Δt,2Δt (unfiltered at t+Δt)
            Aupdate = Aold - two*A[l,m,lf] + Anew   # Eq. 8&9 in William (2009), calculate only once
            A[l,m,1] = A[l,m,lf] + w1*Aupdate       # Robert's filter: A[l,m,1] becomes 2xfiltered value at t
            A[l,m,2] = Anew - w2*Aupdate            # Williams' filter: A[l,m,2] becomes 1xfiltred value at t+Δt
        end
    end
end

"""Leapfrog! for 3D arrays that loops over all vertical layers."""
function leapfrog!( A::AbstractArray{Complex{NF},4},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},3}, # tendency (dynamics + physics) of A
                    dt::Real,                               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    C::Constants{NF},                       # struct containing all constants used at runtime
                    lf::Int=2                               # leapfrog index to dis/enable(default) William's filter
                    ) where {NF<:AbstractFloat}             # number format NF

    for k in 1:size(A)[end]                         # loop over vertical levels (last dimension)
        A_layer = view(A,:,:,:,k)                   # extract vertical layers as views to not allocate any memory
        tendency_layer = view(tendency,:,:,k)
        leapfrog!(A_layer,tendency_layer,dt,C,lf)   # make a timestep forward for each layer
    end
end

"""TODO write a leapfrog! function that loops over all prognostic variables."""
function leapfrog!( progn::PrognosticVariables{NF},         # all prognostic variables
                    tend::Tendencies{NF},                   # all tendencies
                    dt::Real,                               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    C::Constants{NF},                       # struct containing all constants used at runtime
                    lf::Int) where NF                       # leapfrog index to dis/enable William's filter        
    
    # convert fields in progn, tend to a generator tuple for for-loop
    all_progn_variables = (getproperty(progn,prop) for prop in propertynames(progn))
    all_tend_variables  = (getproperty(tend, prop) for prop in propertynames(tend))

    for (var,var_tend) in zip(all_progn_variables,all_tend_variables)
        leapfrog!(var,var_tend,dt,C,lf)
    end 
end    

"""Call initialization of semi-implicit scheme and perform initial time step."""
function first_timesteps!(  progn::PrognosticVariables{NF}, # all prognostic variables
                            diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                            M::ModelSetup                   # everything that is constant at runtime
                            ) where NF
    
    @unpack Δt,Δt_sec = M.constants
    time_sec = 0

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    # IMP = initialize_implicit(half*Δt)
    lf1 = 1     # without Robert+William's filter
    lf2 = 1     # evaluates all tendencies at t=0, the first leapfrog index (=>Euler forward)
    timestep!(progn,diagn,Δt/2,M,lf1,lf2)
    time_sec += Δt_sec÷2

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt)
    # IMP = initialize_implicit(Δt)
    lf1 = 1     # without Robert+William's filter
    lf2 = 2     # evaluate all tendencies at t=dt/2, the 2nd leapfrog index (=>Leapfrog)
    timestep!(progn,diagn,Δt,M,lf1,lf2)
    time_sec += Δt_sec÷2

    # Initialize implicit arrays for further time steps (dt=2Δt)
    # IMP = initialize_implicit(2Δt)
    # return IMP
    return time_sec
end

"""Calculate a single time step for SpeedyWeather.jl"""
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    M::ModelSetup,                  # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                    lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                    ) where {NF<:AbstractFloat}
    
    @unpack vor = progn
    @unpack vor_tend = diagn.tendencies
    @unpack damping, damping_impl = M.horizontal_diffusion

    # set all tendencies to zero
    fill!(vor_tend,zero(Complex{NF}))

    # @unpack tcorh,tcorv,qcorh,qcorv = HD
    # @unpack sdrag = C

    # PROPAGATE THE SPECTRAL STATE INTO THE DIAGNOSTIC VARIABLES
    gridded!(diagn,progn,M,lf2)

    # COMPUTE TENDENCIES OF PROGNOSTIC VARIABLES
    get_tendencies!(diagn,progn,M,lf2)                   

    # DIFFUSION FOR WIND
    vor_lf = view(vor,:,:,1,:)                                      # array view for leapfrog index
    # div_l1 = view(div,:,:,1,:)                                    # TODO l1/l2 dependent?
    horizontal_diffusion!(vor_tend,vor_lf,damping,damping_impl)     # diffusion of vorticity
    # horizontal_diffusion!(div_l1,div_tend,dmpd,dmp1d)             # diffusion of divergence

    # # DIFFUSION FOR TEMPERATURE
    # orographic_correction!(temp_corrected,temp,1,tcorh,tcorv)   # orographic correction for temperature
    # horizontal_diffusion!(temp_corrected,temp_tend,dmp,dmp1)    # diffusion for corrected abs temperature

    # # DISSIPATION IN THE STRATOSPHERE
    # stratospheric_zonal_drag!(vor,vor_tend,sdrag)               # zonal drag for wind
    # stratospheric_zonal_drag!(div,div_tend,sdrag)

    # horizontal_diffusion!(vor_l1,vor_tend,dmps,dmp1s)           # stratospheric diffusion for wind
    # horizontal_diffusion!(div_l1,div_tend,dmps,dmp1s)           
    # horizontal_diffusion!(temp_corrected,temp_tend,dmps,dmp1s)  # stratospheric diffusion for temperature

    # # DIFFUSION OF HUMIDITY
    # orographic_correction!(humid_corrected,humid,1,qcorh,qcorv) # orographic correction for humidity
    # horizontal_diffusion!(humid_corrected,humid_tend,dmp,dmp1)  # diffusion for corrected humidity

    # if ntracers > 1
    #     for i in 2:ntracers
    #         tracer          = view(tracers,:,:,:,1,i)                   # the i-th tracer, leapfrog index 1
    #         tracer_tendency = view(tracers_tendencies,:,:,:,i)          # its tendency
    #         horizontal_diffusion!(tracer,tracer_tendency,dmp,dmp1)
    #         leapfrog!(tracer,tracer_tendency,j1,C)
    #     end
    # end

    # # SPECTRAL TRUNCATION of all tendencies to the spectral resolution
    # for tendency in (logp0_tend,vor_tend,div_tend,temp_tend,humid_tend)
    #     spectral_truncation!(tendency,G)
    # end

    # time integration via leapfrog step forward (filtered with Robert+William's depending on l1)
    # leapfrog!(progn,diagn,dt,C,lf1)
    leapfrog!(vor,vor_tend,dt,M.constants,lf1)
end

"""Calculate a single time step for SpeedyWeather.jl"""
function time_stepping!(progn::PrognosticVariables{NF}, # all prognostic variables
                        diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                        M::ModelSetup{NF}               # all precalculated structs
                        ) where {NF<:AbstractFloat}     # number format NF
    
    @unpack n_timesteps, Δt, Δt_sec = M.constants
    @unpack output = M.parameters

    progn.vor[4,3,1,:] .= 5e-6
    lmax, mmax = 15,15
    progn.vor[1:lmax,1:mmax,1,:] .+= 1e-7*randn(Complex{NF},lmax,mmax,M.parameters.nlev)
    spectral_truncation!(progn.vor[:,:,1,:],M.parameters.trunc)
    
    gridded!(diagn,progn,M,1)

    # FEEDBACK, OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS
    feedback = initialize_feedback(M)
    netcdf_file = initialize_netcdf_output(diagn,feedback,M)

    # FIRST TIMESTEP: EULER FORWARD THEN LEAPFROG IN MAIN LOOP
    time_sec = first_timesteps!(progn,diagn,M)

    for i in 1:n_timesteps
        time_sec += Δt_sec
        timestep!(progn,diagn,2Δt,M)

        # FEEDBACK AND OUTPUT
        feedback!(feedback,i)
        write_netcdf_output!(netcdf_file,feedback,i,time_sec,diagn,M)
    end

    feedback_end!(feedback)

    return progn
end
    