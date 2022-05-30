"""
Struct holding the parameters needed at runtime in number format NF.
"""
@with_kw struct Constants{NF<:AbstractFloat}

    # PHYSICAL CONSTANTS
    radius_earth::NF        # Radius of Earth
    rotation_earth::NF      # Angular frequency of Earth's rotation
    gravity::NF             # Gravitational acceleration
    akap::NF                # Ratio of gas constant to specific heat of dry air at constant pressure
    R_gas::NF               # Universal gas constant

    # TIME STEPPING
    Δt::NF                  # time step [s/m], use 2Δt for leapfrog, scaled by Earth's radius
    Δt_unscaled::NF         # time step [s], as Δt but not scaled with Earth's radius
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

    # PARAMETRIZATIONS
    # Large-scale condensation
    humid_relax_time::NF    # Relaxation time for humidity (hours)
    RH_thresh_max ::NF      # Maximum relative humidity threshold (at σ = 1)
    ΔRH::NF                 # Vertical range of relative humidity threshold
    RH_thresh_boundary::NF  # Relative humidity threshold for boundary layer

end

"""
Generator function for a Constants struct.
"""
function Constants(P::Parameters)

    # PHYSICAL CONSTANTS
    @unpack radius_earth, rotation_earth, gravity, akap, R_gas = P

    # TIME INTEGRATION CONSTANTS
    @unpack robert_filter, williams_filter = P
    @unpack trunc, Δt_at_T85, n_days, output_dt = P

    # PARAMETRIZATION CONSTANTS
    @unpack humid_relax_time, RH_thresh_max, ΔRH, RH_thresh_boundary = P  # Large-scale condensation

    Δt_min_at_trunc = Δt_at_T85*(85/trunc)      # scale time step Δt to specified resolution
    Δt      = round(Δt_min_at_trunc*60)         # convert time step Δt from minutes to whole seconds
    Δt_sec  = convert(Int,Δt)                   # encode time step Δt [s] as integer
    Δt_hrs  = Δt/3600                           # convert time step Δt from minutes to hours
    n_timesteps = ceil(Int,24*n_days/Δt_hrs)    # number of time steps to integrate for
    output_every_n_steps = max(1,floor(Int,output_dt/Δt_hrs))   # output every n time steps
    n_outputsteps = (n_timesteps ÷ output_every_n_steps)+1      # total number of output time steps

    # stratospheric drag [1/s] from damping_time_strat [hrs]
    @unpack damping_time_strat = P
    drag_strat = 1/(damping_time_strat*3600)

    # SCALING
    Δt_unscaled = Δt        # [s] not scaled
    Δt /= radius_earth      # scale with Earth's radius

    # This implies conversion to NF
    return Constants{P.NF}( radius_earth,rotation_earth,gravity,akap,R_gas,
                            Δt,Δt_unscaled,Δt_sec,Δt_hrs,
                            robert_filter,williams_filter,n_timesteps,
                            output_every_n_steps, n_outputsteps,
                            drag_strat, humid_relax_time, RH_thresh_max, ΔRH,
                            RH_thresh_boundary)
end
