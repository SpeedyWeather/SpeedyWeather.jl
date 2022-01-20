"""
Struct containing relevant parameters and arrays nneded for humidity calculations. See humidity.jl
"""
struct HumidityConstants{NF<:AbstractFloat}     

#Set of constants used in the saturation humidity calculations
#Where do these come from? What do they represent?

e0 ::NF
c1 ::NF
c2 ::NF
t0 ::NF 
t1 ::NF 
t2 ::NF 

c3 ::NF
c4 ::NF 

end





"""
Struct holding the parameters needed at runtime in number format NF.
"""
@with_kw struct Constants{NF<:AbstractFloat}
    R_earth::NF      # Radius of Earth
    Ω::NF            # Angular frequency of Earth's rotation
    gravity::NF      # Gravitational acceleration
    akap::NF         # Ratio of gas constant to specific heat of dry air at constant pressure
    R::NF            # Gas constantσ

    # TIME STEPPING
    Δt::NF                  # time step [s] 
    robert_filter::NF       # Robert (1966) time filter coefficient to suppress comput. mode
    williams_filter::NF     # Williams time filter (Amezcua 2011) coefficient for 3rd order acc

    # DIFFUSION AND DRAG
    sdrag::NF               # drag [1/s] for zonal wind in the stratosphere
end

"""
Generator function for a Constants struct.
"""
function Constants{NF}(P::Parameters) where NF      # number format NF

    @unpack R_earth, Ω, gravity, akap, R, Δt = P
    @unpack robert_filter, williams_filter = P
    @unpack tdrs = P

    # stratospheric drag [1/s] from drag time timescale tdrs [hrs]
    sdrag = 1/(tdrs*3600)

    # This implies conversion to NF
    return Constants{NF}(   R_earth,Ω,gravity,akap,R,
                            Δt,robert_filter,williams_filter,
                            sdrag)
end  