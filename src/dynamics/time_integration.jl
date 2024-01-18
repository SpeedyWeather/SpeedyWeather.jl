"""
Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Leapfrog{NF} <: TimeStepper{NF}

    # DIMENSIONS
    "spectral resolution (max degree of spherical harmonics)"
    trunc::Int                      

    # OPTIONS
    "time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::Second = Minute(30)

    "radius of sphere [m], used for scaling"
    radius::NF = DEFAULT_RADIUS

    "adjust Δt_at_T31 with the output_dt to reach output_dt exactly in integer time steps"
    adjust_with_output::Bool = true

    # NUMERICS
    "Robert (1966) time filter coefficeint to suppress comput. mode"
    robert_filter::NF = 0.05

    "Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF = 0.53

    # DERIVED FROM OPTIONS    
    "time step Δt [ms] at specified resolution"
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value/1000

    "time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec/radius  
end

"""
$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T31` is additionally 
adjusted to the closest divisor of `output_dt` so that the output time axis is keeping
`output_dt` exactly.
"""
function get_Δt_millisec(
    Δt_at_T31::Dates.TimePeriod,
    trunc,
    radius,
    adjust_with_output::Bool,
    output_dt::Dates.TimePeriod = DEFAULT_OUTPUT_DT,
)
    # linearly scale Δt with trunc+1 (which are often powers of two)
    resolution_factor = (DEFAULT_TRUNC+1)/(trunc+1)

    # radius also affects grid spacing, scale proportionally
    radius_factor = radius/DEFAULT_RADIUS

    # maybe rename to _at_trunc_and_radius?
    Δt_at_trunc = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (output_dt > Millisecond(0))
        k = round(Int,Second(output_dt).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(output_dt).value)
        sort!(divisors)
        i = findfirst(x -> x>=k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(output_dt).value/k_new))

        # provide info when time step is significantly shortened or lengthened
        Δt_millisec_unadjusted = round(Int,1000*Δt_at_trunc)
        Δt_ratio = Δt_millisec.value/Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05     # only when +-5% changes
            p = round(Int,(Δt_ratio - 1)*100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else 
        Δt_millisec = Millisecond(round(Int, 1000*Δt_at_trunc))
    end

    return Δt_millisec
end 

"""
$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(spectral_grid::SpectralGrid;kwargs...)
    (;NF,trunc,radius) = spectral_grid
    return Leapfrog{NF}(;trunc,radius,kwargs...)
end

function Base.show(io::IO,L::Leapfrog)
    println(io,"$(typeof(L)) <: TimeStepper")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

"""
$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the timestep given the output time step
`output_dt` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog,model::ModelSetup)
    (;output_dt) = model.output

    if L.adjust_with_output
        # take actual output dt from model.output and recalculate timestep
        L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, L.radius, L.adjust_with_output, output_dt)
        L.Δt_sec = L.Δt_millisec.value/1000
        L.Δt = L.Δt_sec/L.radius
    end

    # check how time stepping time step and output time step align
    n = round(Int,Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @info "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms (=$(nΔt.value/1000)s), but output_dt = $(output_dt.value)s"
    end
end

"""
$(TYPEDSIGNATURES)
Performs one leapfrog time step with (`lf=2`) or without (`lf=1`) Robert+Williams filter
(see Williams (2009), Montly Weather Review, Eq. 7-9)."""
function leapfrog!( A_old::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t
                    A_new::LowerTriangularMatrix{Complex{NF}},      # prognostic variable at t+dt
                    tendency::LowerTriangularMatrix{Complex{NF}},   # tendency (dynamics+physics) of A
                    dt::Real,                                       # time step (=2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                                        # leapfrog index to dis/enable Williams filter
                    L::Leapfrog{NF},                                # struct with constants
                    ) where {NF<:AbstractFloat}                     # number format NF

    @boundscheck lf == 1 || lf == 2 || throw(BoundsError())         # index lf picks leapfrog dim
    
    A_lf = lf == 1 ? A_old : A_new                      # view on either t or t+dt to dis/enable Williams filter        
    (;robert_filter, williams_filter) = L               # coefficients for the Robert and Williams filter
    dt_NF = convert(NF,dt)                              # time step dt in number format NF

    # LEAP FROG time step with or without Robert+Williams filter
    # Robert time filter to compress computational mode, Williams filter for 3rd order accuracy
    # see Williams (2009), Eq. 7-9
    # for lf == 1 (initial time step) no filter applied (w1=w2=0)
    # for lf == 2 (later steps) Robert+Williams filter is applied
    w1 = lf == 1 ? zero(NF) : robert_filter*williams_filter/2       # = ν*α/2 in Williams (2009, Eq. 8)
    w2 = lf == 1 ? zero(NF) : robert_filter*(1-williams_filter)/2   # = ν(1-α)/2 in Williams (2009, Eq. 9)

    @inbounds for lm in eachharmonic(A_old,A_new,A_lf,tendency)
        a_old = A_old[lm]                       # double filtered value from previous time step (t-Δt)
        a_new = a_old + dt_NF*tendency[lm]      # Leapfrog/Euler step depending on dt=Δt,2Δt (unfiltered at t+Δt)
        a_update = a_old - 2A_lf[lm] + a_new # Eq. 8&9 in Williams (2009), calculate only once
        A_old[lm] = A_lf[lm] + w1*a_update      # Robert's filter: A_old[lm] becomes 2xfiltered value at t
        A_new[lm] = a_new - w2*a_update         # Williams filter: A_new[lm] becomes 1xfiltered value at t+Δt
    end
end

# variables that are leapfrogged in the respective models that are on layers (so excl surface pressure)
leapfrog_layer_vars(::Barotropic) = (:vor,)
leapfrog_layer_vars(::ShallowWater) = (:vor, :div)
leapfrog_layer_vars(::PrimitiveDry) = (:vor, :div, :temp)
leapfrog_layer_vars(::PrimitiveWet) = (:vor, :div, :temp, :humid)

function leapfrog!( progn::PrognosticLayerTimesteps,
                    diagn::DiagnosticVariablesLayer,
                    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                # leapfrog index to dis/enable Williams filter
                    model::ModelSetup)
               
    for var in leapfrog_layer_vars(model)
        var_old = getproperty(progn.timesteps[1],var)
        var_new = getproperty(progn.timesteps[2],var)
        var_tend = getproperty(diagn.tendencies,Symbol(var,:_tend))
        spectral_truncation!(var_tend)      # set lmax+1 mode to zero
        leapfrog!(var_old,var_new,var_tend,dt,lf,model.time_stepping)
    end
end

function leapfrog!( progn::PrognosticSurfaceTimesteps,
                    diagn::SurfaceVariables,
                    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    lf::Int,                # leapfrog index to dis/enable Williams filter
                    model::ModelSetup)
               
    (;pres_tend) = diagn
    pres_old = progn.timesteps[1].pres
    pres_new = progn.timesteps[2].pres
    spectral_truncation!(pres_tend)         # set lmax+1 mode to zero
    leapfrog!(pres_old,pres_new,pres_tend,dt,lf,model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
prognostic variables with two time steps (t=0,Δt) that can then be used in the normal leap frogging."""
function first_timesteps!(  
    progn::PrognosticVariables,         # all prognostic variables
    diagn::DiagnosticVariables,         # all pre-allocated diagnostic variables
    model::ModelSetup,                  # everything that is constant at runtime
    output::AbstractOutputWriter,
)
    
    (;clock) = progn
    clock.n_timesteps == 0 && return nothing    # exit immediately for no time steps
    
    (;implicit) = model
    (;Δt, Δt_millisec) = model.time_stepping
    Δt_millisec_half = Dates.Millisecond(Δt_millisec.value÷2)   # this might be 1ms off

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    lf1 = 1                             # without Robert+Williams filter
    lf2 = 1                             # evaluates all tendencies at t=0,
                                        # the first leapfrog index (=>Euler forward)
    initialize!(implicit,Δt/2,diagn,model)  # update precomputed implicit terms with time step Δt/2
    timestep!(progn,diagn,Δt/2,model,lf1,lf2)
    clock.time += Δt_millisec_half      # update by half the leapfrog time step Δt used here

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
    initialize!(implicit,Δt,diagn,model)    # update precomputed implicit terms with time step Δt
    lf1 = 1                             # without Robert+Williams filter
    lf2 = 2                             # evaluate all tendencies at t=dt/2,
                                        # the 2nd leapfrog index (=>Leapfrog)
    timestep!(progn,diagn,Δt,model,lf1,lf2)
    clock.time -= Δt_millisec_half      # remove prev Δt/2 in case not even ms
    clock.time += Δt_millisec           # otherwise time is off by 1ms
    write_output!(output,clock.time,diagn)

    # from now on precomputed implicit terms with 2Δt
    initialize!(implicit,2Δt,diagn,model) 

    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model <: Barotropic`."""
function timestep!( progn::PrognosticVariables,     # all prognostic variables
                    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    model::Barotropic,              # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+Williams filter)
                    lf2::Int=2)                     # leapfrog index 2 (time step used for tendencies)

    model.feedback.nars_detected && return nothing  # exit immediately if NaRs already present
    (;time) = progn.clock                           # current time
    zero_tendencies!(diagn)                         # start with zero for accumulation 

    # LOOP OVER LAYERS FOR TENDENCIES, DIFFUSION, LEAPFROGGING AND PROPAGATE STATE TO GRID
    for (progn_layer,diagn_layer) in zip(progn.layers,diagn.layers)
        progn_lf = progn_layer.timesteps[lf2]       # pick the leapfrog time step lf2 for tendencies
        dynamics_tendencies!(diagn_layer,progn_lf,time,model)
        horizontal_diffusion!(diagn_layer,progn_layer,model)
        leapfrog!(progn_layer,diagn_layer,dt,lf1,model)
        gridded!(diagn_layer,progn_lf,model)
    end
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model <: ShallowWater`."""
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    model::ShallowWater,            # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+Williams filter)
                    lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                    ) where {NF<:AbstractFloat}

    model.feedback.nars_detected && return nothing  # exit immediately if NaRs already present
    (;time) = progn.clock                           # current time

    progn_layer = progn.layers[1]                   # only calculate tendencies for the first layer
    diagn_layer = diagn.layers[1]                   # multi-layer shallow water not supported

    progn_lf = progn_layer.timesteps[lf2]           # pick the leapfrog time step lf2 for tendencies
    (;pres) = progn.surface.timesteps[lf2]
    (;implicit, spectral_transform) = model

    zero_tendencies!(diagn)
    
    # GET TENDENCIES, CORRECT THEM FOR SEMI-IMPLICIT INTEGRATION
    dynamics_tendencies!(diagn_layer,progn_lf,diagn.surface,pres,time,model)
    implicit_correction!(diagn_layer,progn_layer,diagn.surface,progn.surface,implicit)
    
    # APPLY DIFFUSION, STEP FORWARD IN TIME, AND TRANSFORM NEW TIME STEP TO GRID
    horizontal_diffusion!(progn_layer,diagn_layer,model)
    leapfrog!(progn_layer,diagn_layer,dt,lf1,model)
    gridded!(diagn_layer,progn_lf,model)

    # SURFACE LAYER (pressure), no diffusion though
    (;pres_grid) = diagn.surface
    leapfrog!(progn.surface,diagn.surface,dt,lf1,model)
    gridded!(pres_grid,pres,spectral_transform)
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model<:PrimitiveEquation`"""
function timestep!( progn::PrognosticVariables{NF}, # all prognostic variables
                    diagn::DiagnosticVariables{NF}, # all pre-allocated diagnostic variables
                    dt::Real,                       # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    model::PrimitiveEquation,       # everything that's constant at runtime
                    lf1::Int=2,                     # leapfrog index 1 (dis/enables Robert+Williams filter)
                    lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                    ) where {NF<:AbstractFloat}

    model.feedback.nars_detected && return nothing  # exit immediately if NaRs already present
    (;time) = progn.clock                           # current time

    if model.physics                                # switch on/off all physics parameterizations
        # time step ocean (temperature and TODO sea ice) and land (temperature and soil moisture)
        ocean_timestep!(progn.ocean,time,model)
        land_timestep!(progn.land,time,model)
        soil_moisture_availability!(diagn.surface,progn.land,model)

        # calculate all parameterizations
        parameterization_tendencies!(diagn,progn,time,model)
    else                                            # set tendencies to zero otherwise for accumulators
        zero_tendencies!(diagn)
    end

    if model.dynamics                               # switch on/off all dynamics
        dynamics_tendencies!(diagn,progn,model,lf2)         # dynamical core
        implicit_correction!(diagn,model.implicit,progn)    # semi-implicit time stepping corrections
    else    # just transform physics tendencies to spectral space
        for k in 1:diagn.nlev
            diagn_layer = diagn.layers[k]
            tendencies_physics_only!(diagn_layer,model)
        end
    end

    # LOOP OVER ALL LAYERS for diffusion, leapfrog time integration
    # and progn state from spectral to grid for next time step
    @floop for k in 1:diagn.nlev+1
        if k <= diagn.nlev                  # model levels
            diagn_layer = diagn.layers[k]
            progn_layer = progn.layers[k]
            progn_layer_lf = progn_layer.timesteps[lf2]

            horizontal_diffusion!(progn_layer,diagn_layer,model)    # implicit diffusion of vor, div, temp
            leapfrog!(progn_layer,diagn_layer,dt,lf1,model)         # time step forward for vor, div, temp
            gridded!(diagn_layer,progn_layer_lf,model)              # propagate spectral state to grid
        else                                                        # surface level
            leapfrog!(progn.surface,diagn.surface,dt,lf1,model)
            (;pres_grid) = diagn.surface
            pres_lf = progn.surface.timesteps[lf2].pres
            gridded!(pres_grid,pres_lf,model.spectral_transform)
        end
    end
end

"""
$(TYPEDSIGNATURES)
Main time loop that that initializes output and feedback, loops over all time steps
and calls the output and feedback functions."""
function time_stepping!(
    progn::PrognosticVariables,     # all prognostic variables
    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
    model::ModelSetup,              # all precalculated structs
)          
    
    (;clock) = progn
    (;Δt, Δt_millisec) = model.time_stepping
    (;time_stepping) = model

    # SCALING: we use vorticity*radius,divergence*radius in the dynamical core
    scale!(progn,model.spectral_grid.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    (;output,feedback) = model
    lf = 1                                  # use first leapfrog index
    gridded!(diagn,progn,lf,model)
    initialize!(output,feedback,time_stepping,clock,diagn,model)
    initialize!(feedback,clock,model)

    # FIRST TIMESTEPS: EULER FORWARD THEN 1x LEAPFROG
    first_timesteps!(progn,diagn,model,output)

    # MAIN LOOP
    for i in 2:clock.n_timesteps            # start at 2 as first Δt in first_timesteps!
        
        # calculate tendencies and leapfrog forward
        timestep!(progn,diagn,2Δt,model)   
        clock.time += Δt_millisec           # time of lf=2 and diagn after timestep!

        progress!(feedback,progn)           # updates the progress meter bar
        write_output!(output,clock.time,diagn)
    end

    # UNSCALE, CLOSE, FINISH
    unscale!(progn)                         # undo radius-scaling for vor,div from the dynamical core
    close(output)                           # close netCDF file
    write_restart_file(progn,output)        # as JLD2 
    progress_finish!(feedback)              # finishes the progress meter bar

    return progn
end