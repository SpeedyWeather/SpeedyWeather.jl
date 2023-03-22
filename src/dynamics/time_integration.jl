"""
    leapfrog!(  A_old::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t
                A_new::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t+dt
                tendency::LowerTriangularMatrix{Complex{NF}},   # tendency (dynamics+physics) of A
                dt::Real,                                       # time step (=2Δt, but for init steps =Δt,Δt/2)
                lf::Int=2,                                      # leapfrog index to dis/enable William's filter
                C::DynamicsConstants{NF},                       # struct with constants used at runtime
                ) where {NF<:AbstractFloat}                     # number format NF

Performs one leapfrog time step with (`lf=2`) or without (`lf=1`) Robert+William's filter
(see William (2009), Montly Weather Review, Eq. 7-9).
"""
function leapfrog!( A_old::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t
                    A_new::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t+dt
                    tendency::LowerTriangularMatrix{Complex{NF}},   # tendency (dynamics+physics) of A
                    dt::Real,                                       # time step (=2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                                        # leapfrog index to dis/enable William's filter
                    C::DynamicsConstants{NF},                               # struct with constants used at runtime
                    ) where {NF<:AbstractFloat}                     # number format NF

    @boundscheck lf == 1 || lf == 2 || throw(BoundsError())         # index lf picks leapfrog dim
    
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

# variables that are leapfrogged in the respective models that are on layers (so excl surface pressure)
leapfrog_layer_vars(::Barotropic) = (:vor,)
leapfrog_layer_vars(::ShallowWater) = (:vor, :div)
leapfrog_layer_vars(::PrimitiveDryCore) = (:vor, :div, :temp)
leapfrog_layer_vars(::PrimitiveWetCore) = (:vor, :div, :temp, :humid)

function leapfrog!( progn::PrognosticVariablesLeapfrog,
                    diagn::DiagnosticVariablesLayer,
                    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                # leapfrog index to dis/enable William's filter
                    model::ModelSetup)
               
    for var in leapfrog_layer_vars(model)
        var_old = getproperty(progn.leapfrog[1],var)
        var_new = getproperty(progn.leapfrog[2],var)
        var_tend = getproperty(diagn.tendencies,Symbol(var,:_tend))
        leapfrog!(var_old,var_new,var_tend,dt,lf,model.constants)
    end
end

"""
    first_timesteps!(   progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        time::DateTime,             # time at timestep
                        M::ModelSetup,              # everything that is constant at runtime
                        feedback::AbstractFeedback  # feedback struct
                        )

Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
prognostic variables with two time steps (t=0,Δt) that can then be used in the normal leap frogging."""
function first_timesteps!(  progn::PrognosticVariables, # all prognostic variables
                            diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                            time::DateTime,             # time at timestep
                            model::ModelSetup,          # everything that is constant at runtime
                            feedback::AbstractFeedback, # feedback struct
                            outputter::AbstractOutput
                            )
    
    @unpack n_timesteps, Δt, Δt_sec = model.constants
    n_timesteps == 0 && return time     # exit immediately for no time steps

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    initialize_implicit!(Δt/2,model)    # update precomputed implicit terms with time step Δt/2
    lf1 = 1                             # without Robert+William's filter
    lf2 = 1                             # evaluates all tendencies at t=0,
                                        # the first leapfrog index (=>Euler forward)
    timestep!(progn,diagn,time,Δt/2,model,lf1,lf2)
    time += Dates.Second(Δt_sec÷2)      # update by half the leapfrog time step Δt used here
    progress!(feedback,progn)

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
    initialize_implicit!(Δt,model)      # update precomputed implicit terms with time step Δt
    lf1 = 1                             # without Robert+William's filter
    lf2 = 2                             # evaluate all tendencies at t=dt/2,
                                        # the 2nd leapfrog index (=>Leapfrog)
    timestep!(progn,diagn,time,Δt,model,lf1,lf2)
    time += Dates.Second(Δt_sec÷2)      # now 2nd leapfrog step is at t=Δt
    progress!(feedback,progn)
    write_netcdf_output!(outputter,time,diagn,model)

    return time
end

"""
    timestep!(  progn::PrognosticVariables,     # all prognostic variables
                diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                time::DateTime,                 # time at timestep
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2,                     # leapfrog index 2 (time step used for tendencies)
                M::BarotropicModel,             # everything that's constant at runtime
                )

Calculate a single time step for the barotropic vorticity equation model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables,     # all prognostic variables
                    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                    time::DateTime,                 # time at time step 
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
    timestep!(  progn::PrognosticVariables{NF}, # all prognostic variables
                diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                time::DateTime,                 # time at timestep
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                M::ShallowWaterModel,           # everything that's constant at runtime
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                ) where {NF<:AbstractFloat}

Calculate a single time step for the shallow water model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    time::DateTime,                 # time at timestep
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
    pres_lf = pres.leapfrog[lf2]

    get_tendencies!(diagn_layer,diagn_surface,pres_lf,time,M)               # tendency of vor, div, pres
    implicit_correction!(diagn_layer,progn_layer,diagn_surface,pres,M)      # dampen gravity waves
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
    timestep!(  progn::PrognosticVariables{NF}, # all prognostic variables
                diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                time::DateTime,                 # time at timestep
                dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                M::PrimitiveEquation,      # everything that's constant at runtime
                lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                ) where {NF<:AbstractFloat}

Calculate a single time step for the primitive equation model of SpeedyWeather.jl """
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    time::DateTime,                 # time at timestep
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    model::PrimitiveEquation,       # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+William's filter)
                    lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                    ) where {NF<:AbstractFloat}

    get_tendencies!(diagn,progn,time,model,lf2)
    implicit_correction!(diagn,progn,model)

    # LOOP OVER ALL LAYERS for diffusion, leapfrog time integration
    # and progn state from spectral to grid for next time step
    @threads for k in 1:diagn.nlev+1
        if k <= diagn.nlev
            diagn_layer = diagn.layers[k]
            progn_layer = progn.layers[k]

            horizontal_diffusion!(progn_layer,diagn_layer,model)    # implicit diffusion of vor, div, temp
            leapfrog!(progn_layer,diagn_layer,dt,lf1,model)         # time step forward for vor, div, temp
            gridded!(diagn_layer,progn_layer,lf2,model)             # propagate spectral state to grid
        else
            # SURFACE LAYER (log of surface pressure)
            @unpack pres_tend = diagn.surface
            pres_old,pres_new = progn.pres.leapfrog
            leapfrog!(pres_old,pres_new,pres_tend,dt,lf1,model.constants)
            gridded!(diagn.surface.pres_grid,progn.pres.leapfrog[lf2],model.spectral_transform)
        end
    end
end

"""
    time_stepping!( progn::PrognosticVariables,     # all prognostic variables
                    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                    model::ModelSetup)              # all precalculated structs

Main time loop that that initialises output and feedback, loops over all time steps
and calls the output and feedback functions."""
function time_stepping!(progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        model::ModelSetup)          # all precalculated structs
    
    @unpack n_timesteps, Δt, Δt_sec = model.constants
    time = model.parameters.startdate

    # SCALING: we use vorticity*radius,divergence*radius in the dynamical core
    scale!(progn,model)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    lf = 1                                  # use first leapfrog index
    gridded!(diagn,progn,lf,model)
    outputter = initialize_netcdf_output(diagn,model)
    feedback = initialize_feedback(outputter,model)

    # FIRST TIMESTEPS: EULER FORWARD THEN 1x LEAPFROG
    time = first_timesteps!(progn,diagn,time,model,feedback,outputter)
    initialize_implicit!(2Δt,model)         # from now on precomputed implicit terms with 2Δt

    # MAIN LOOP
    for i in 2:n_timesteps                  # start at 2 as first Δt in first_timesteps!
        timestep!(progn,diagn,time,2Δt,model)   # calculate tendencies and leapfrog forward
        time += Dates.Second(Δt_sec)        # time of lf=2 and diagn after timestep!

        progress!(feedback,progn)           # updates the progress meter bar
        write_netcdf_output!(outputter,time,diagn,model)
    end

    unscale!(progn,model)                   # unscale radius-scaling from the dynamical core
    write_restart_file(time,progn,outputter,model)
    progress_finish!(feedback)              # finishes the progress meter bar

    return progn
end