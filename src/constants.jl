"""
Struct holding the parameters needed at runtime in number format NF.
"""
@with_kw struct Constants{NF<:AbstractFloat}
    R_earth::NF      # Radius of Earth
    Ω::NF            # Angular frequency of Earth's rotation
    gravity::NF      # Gravitational acceleration
    akap::NF         # Ratio of gas constant to specific heat of dry air at constant pressure
    R_gas::NF        # Universal gas constant

    # TIME STEPPING
    Δt::NF                  # time step [s], use 2Δt for leapfrog
    robert_filter::NF       # Robert (1966) time filter coefficient to suppress comput. mode
    williams_filter::NF     # Williams time filter (Amezcua 2011) coefficient for 3rd order acc

    # DIFFUSION AND DRAG
    drag_strat::NF               # drag [1/s] for zonal wind in the stratosphere
end

"""
Generator function for a Constants struct.
"""
function Constants(P::Parameters)

    @unpack R_earth, Ω, gravity, akap, R_gas, Δt = P
    @unpack robert_filter, williams_filter = P
    @unpack damping_time_strat = P

    # stratospheric drag [1/s] from drag time timescale tdrs [hrs]
    drag_strat = 1/(damping_time_strat*3600)

    # This implies conversion to NF
    return Constants{P.NF}( R_earth,Ω,gravity,akap,R_gas,
                            Δt,robert_filter,williams_filter,
                            drag_strat)
end  