"""
Struct holding the parameters needed at runtime in number format NF.
"""
@with_kw struct Constants{NF<:AbstractFloat}
    R_earth::NF      # Radius of Earth
    Ω::NF            # Angular frequency of Earth's rotation
    g::NF            # Gravitational acceleration
    akap::NF         # Ratio of gas constant to specific heat of dry air at constant pressure
    R::NF            # Gas constant
    γ::NF            # Reference temperature lapse rate (-dT/dz in deg/km)
    hscale::NF       # Reference scale height for pressure (in km)
    hshum::NF        # Reference scale height for specific humidity (in km)
    rh_ref::NF       # Reference relative humidity of near-surface air

    # TIME STEPPING
    robert_filter::NF       # Robert (1966) time filter coefficient to suppress comput. mode
    williams_filter::NF     # Williams time filter (Amezcua 2011) coefficient for 3rd order acc
end

"""
Generator function for a Constants struct.
"""
function Constants{NF}(P::Params) where NF      # number format NF

    @unpack R_earth, Ω, g, akap, R, γ, hscale, hshum, rh_ref = P
    @unpack robert_filter, williams_filter = P

    # This implies conversion to NF
    return Constants{NF}(R_earth,Ω,g,akap,R,γ,hscale,hshum,rh_ref,
                            robert_filter,williams_filter)
end
