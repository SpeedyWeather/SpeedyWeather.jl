"""
    leapfrog!(  A_old::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t
                A_new::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t+dt
                tendency::LowerTriangularMatrix{Complex{NF}},   # tendency (dynamics+physics) of A
                dt::Real,                                       # time step (=2Δt, but for init steps =Δt,Δt/2)
                lf::Int=2,                                      # leapfrog index to dis/enable William's filter
                C::Constants{NF},                               # struct with constants used at runtime
                ) where {NF<:AbstractFloat}                     # number format NF

Performs one leapfrog time step with (`lf=2`) or without (`lf=1`) Robert+William's filter
(see William (2009), Montly Weather Review, Eq. 7-9).
"""
function leapfrog!( A_old::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t
                    A_new::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t+dt
                    tendency::LowerTriangularMatrix{Complex{NF}},   # tendency (dynamics+physics) of A
                    dt::Real,                                       # time step (=2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                                        # leapfrog index to dis/enable William's filter
                    C::Constants{NF},                               # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}                     # number format NF

    @boundscheck lf in [1,2] || throw(BoundsError())    # index lf picks leapfrog dim
    
    A_lf = lf == 1 ? A_old : A_new                      # view on either t or t+dt to dis/enable William's filter
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

    @inbounds for lm in eachharmonic(A_old,A_new,A_lf,tendency)
        a_old = A_old[lm]                       # double filtered value from previous time step (t-Δt)
        a_new = a_old + dt_NF*tendency[lm]      # Leapfrog/Euler step depending on dt=Δt,2Δt (unfiltered at t+Δt)
        a_update = a_old - two*A_lf[lm] + a_new # Eq. 8&9 in William (2009), calculate only once
        A_old[lm] = A_lf[lm] + w1*a_update      # Robert's filter: A_old[lm] becomes 2xfiltered value at t
        A_new[lm] = a_new - w2*a_update         # Williams' filter: A_new[lm] becomes 1xfiltred value at t+Δt
    end
end

function leapfrog!( progn::PrognosticVariablesLeapfrog,
                    diagn::DiagnosticVariablesLayer,
                    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                # leapfrog index to dis/enable William's filter
                    M::BarotropicModel)
                    
    vor_old  = progn.leapfrog[1].vor        # unpack vorticity and call leapfrog!
    vor_new  = progn.leapfrog[2].vor
    @unpack vor_tend = diagn.tendencies
    leapfrog!(vor_old,vor_new,vor_tend,dt,lf,M.constants)
end

function leapfrog!( progn::PrognosticVariablesLeapfrog,
                    diagn::DiagnosticVariablesLayer,
                    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                # leapfrog index to dis/enable William's filter
                    M::ShallowWaterModel)
                    
    vor_old  = progn.leapfrog[1].vor        # vorticity
    vor_new  = progn.leapfrog[2].vor
    @unpack vor_tend = diagn.tendencies
    leapfrog!(vor_old,vor_new,vor_tend,dt,lf,M.constants)

    div_old  = progn.leapfrog[1].div        # divergence
    div_new  = progn.leapfrog[2].div
    @unpack div_tend = diagn.tendencies
    leapfrog!(div_old,div_new,div_tend,dt,lf,M.constants)
end

"""
    first_timesteps!(   progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        M::ModelSetup,              # everything that is constant at runtime
                        feedback::Feedback          # feedback struct
                        )

Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
prognostic variables with two time steps (t=0,Δt) that can then be used in the normal leap frogging."""
function first_timesteps!(  progn::PrognosticVariables, # all prognostic variables
                            diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                            M::ModelSetup,              # everything that is constant at runtime
                            feedback::Feedback          # feedback struct
                            )
    
    @unpack Δt,Δt_sec = M.constants
    time_sec = 0    # overall time counter in seconds using Int64 for error free accumulation

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    initialize_implicit!(Δt/2,M)        # update precomputed implicit terms with time step Δt/2
    lf1 = 1                             # without Robert+William's filter
    lf2 = 1                             # evaluates all tendencies at t=0,
                                        # the first leapfrog index (=>Euler forward)
    timestep!(progn,diagn,Δt/2,M,lf1,lf2)
    time_sec += Δt_sec÷2
    progress!(feedback)

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
    initialize_implicit!(Δt,M)          # update precomputed implicit terms with time step Δt
    lf1 = 1                             # without Robert+William's filter
    lf2 = 2                             # evaluate all tendencies at t=dt/2,
                                        # the 2nd leapfrog index (=>Leapfrog)
    timestep!(progn,diagn,Δt,M,lf1,lf2)
    time_sec += Δt_sec÷2                # update by half the leapfrog time step Δt used here
    progress!(feedback)

    # update precomputed implicit terms with time step 2Δt for further time steps
    initialize_implicit!(2Δt,M)
    return time_sec
end

"""
    timestep!(  progn::PrognosticVariables{NF}, # all prognostic variables
                diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                M::PrimitiveEquationModel,      # everything that's constant at runtime
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                ) where {NF<:AbstractFloat}

Calculate a single time step for the primitive equation model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    M::PrimitiveEquationModel,      # everything that's constant at runtime
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
    # always use first leapfrog index for diffusion (otherwise unstable)
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

"""
    timestep!(  progn::PrognosticVariables{NF}, # all prognostic variables
                diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                M::ShallowWaterModel,           # everything that's constant at runtime
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                ) where {NF<:AbstractFloat}

Calculate a single time step for the shallow water model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    M::ShallowWaterModel,           # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                    lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                    ) where {NF<:AbstractFloat}

    # IMPLICIT, DIFFUSION, LEAPFROGGING AND PROPAGATE STATE TO GRID
    progn_layer = progn.layers[1]                       # only calculate tendencies for the first layer
    diagn_layer = diagn.layers[1]
    diagn_surface = diagn.surface
    @unpack pres = progn

    get_tendencies!(diagn_layer,diagn_surface,M)        # tendency of vorticity, divergence, pressure
    implicit_correction!(diagn_layer,progn_layer,diagn_surface,pres,M)     # correct divergence
    horizontal_diffusion!(progn_layer,diagn_layer,M)    # diffusion for vorticity and divergence
    leapfrog!(progn_layer,diagn_layer,dt,lf1,M)         # leapfrog vorticity forward
    gridded!(diagn_layer,progn_layer,lf2,M)             # propagate spectral state to grid

    # SURFACE LAYER (pressure)
    @unpack pres_tend = diagn_surface
    pres_old,pres_new = pres.leapfrog
    leapfrog!(pres_old,pres_new,pres_tend,dt,lf1,M.constants)
    gridded!(diagn.surface.pres_grid,progn.pres.leapfrog[lf2],M.spectral_transform)
end

"""
    timestep!(  progn::PrognosticVariables,     # all prognostic variables
                diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2,                     # leapfrog index 2 (time step used for tendencies)
                M::BarotropicModel,             # everything that's constant at runtime
                )

Calculate a single time step for the barotropic vorticity equation model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables,     # all prognostic variables
                    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    M::BarotropicModel,             # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                    lf2::Int=2,                     # leapfrog index 2 (time step used for tendencies)
                    )

    # LOOP OVER LAYERS FOR DIFFUSION, LEAPFROGGING AND PROPAGATE STATE TO GRID
    for (progn_layer,diagn_layer) in zip(progn.layers,diagn.layers)
        get_tendencies!(diagn_layer,M)                      # tendency of vorticity
        horizontal_diffusion!(progn_layer,diagn_layer,M)    # diffusion for vorticity
        leapfrog!(progn_layer,diagn_layer,dt,lf1,M)         # leapfrog vorticity forward
        gridded!(diagn_layer,progn_layer,lf2,M)             # propagate spectral state to grid
    end
end

"""
    time_stepping!( progn::PrognosticVariables,     # all prognostic variables
                    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                    M::ModelSetup                   # all precalculated structs
                    )

Main time loop that that initialises output and feedback, loops over all time steps
and calls the output and feedback functions."""
function time_stepping!(progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        M::ModelSetup,              # all precalculated structs
                        )
    
    @unpack n_timesteps, Δt, Δt_sec = M.constants
    @unpack output = M.parameters
    
    # propagate spectral state to grid variables for initial condition output
    lf = 1                              # use first leapfrog index
    gridded!(diagn,progn,lf,M)

    # FEEDBACK, OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS
    feedback = initialize_feedback(M)
    netcdf_file = initialize_netcdf_output(progn,diagn,feedback,M)

    # FIRST TIMESTEP: EULER FORWARD THEN LEAPFROG IN MAIN LOOP
    time_sec = first_timesteps!(progn,diagn,M,feedback)

    # MAIN LOOP
    for i in 2:n_timesteps              # 2: as first Δt time step in first_timesteps!
        time_sec += Δt_sec
        timestep!(progn,diagn,2Δt,M)

        # FEEDBACK AND OUTPUT
        progress!(feedback)             # updates the progress meter bar
        write_netcdf_output!(netcdf_file,feedback,time_sec,progn,diagn,M)
    end

    write_restart_file(time_sec,progn,feedback,M)
    progress_finish!(feedback)          # finishes the progress meter bar

    return progn
end
    