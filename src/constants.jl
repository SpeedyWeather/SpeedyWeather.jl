"""
Struct holding the parameters needed at runtime in number format NF.
"""
@with_kw struct Constants{NF<:AbstractFloat}

    # PHYSICAL CONSTANTS
    radius_earth::NF # Radius of Earth
    Ω::NF            # Angular frequency of Earth's rotation
    gravity::NF      # Gravitational acceleration
    akap::NF         # Ratio of gas constant to specific heat of dry air at constant pressure
    R_gas::NF        # Universal gas constant

    # TIME STEPPING
    Δt::NF                  # time step [s], use 2Δt for leapfrog
    Δt_sec::Int             # time step [s] but encoded as 64-bit integer for rounding error-free accumulation
    Δt_hrs::Float64         # time step [hrs]
    robert_filter::NF       # Robert (1966) time filter coefficient to suppress comput. mode
    williams_filter::NF     # Williams time filter (Amezcua 2011) coefficient for 3rd order acc
    n_timesteps::Int        # number of time steps to integrate for

    # OUTPUT TIME STEPPING
    output_every_n_steps::Int   # output every n time steps
    n_outputsteps::Int          # total number of output time steps

    # DIFFUSION AND DRAG
    drag_strat::NF               # drag [1/s] for zonal wind in the stratosphere
end

"""
Generator function for a Constants struct.
"""
function Constants(P::Parameters)

    # PHYSICAL CONSTANTS
    @unpack radius_earth, Ω, gravity, akap, R_gas = P
    
    # TIME INTEGRATION CONSTANTS
    @unpack robert_filter, williams_filter = P
    @unpack n_days, output_dt = P
    Δt      = P.Δt*60                           # convert time step Δt from minutes to seconds
    Δt_sec  = round(Int,Δt)                     # encode time step Δt [s] as integer
    Δt_hrs  = P.Δt/60                           # convert time step Δt from minutes to hours
    n_timesteps = ceil(Int,24*n_days/Δt_hrs)    # number of time steps to integrate for
    output_every_n_steps = max(1,floor(Int,output_dt/Δt_hrs))   # output every n time steps
    n_outputsteps = (n_timesteps ÷ output_every_n_steps)+1      # total number of output time steps

    # stratospheric drag [1/s] from damping_time_strat [hrs]
    @unpack damping_time_strat = P
    drag_strat = 1/(damping_time_strat*3600)

    # This implies conversion to NF
    return Constants{P.NF}( radius_earth,Ω,gravity,akap,R_gas,
                            Δt,Δt_sec,Δt_hrs,
                            robert_filter,williams_filter,n_timesteps,
                            output_every_n_steps, n_outputsteps,
                            drag_strat)
end  