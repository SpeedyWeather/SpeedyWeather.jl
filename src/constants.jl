"""
Struct containing relevant parameters and arrays nneded for humidity calculations. See humidity.jl
"""
struct HumidityConstants{NF<:AbstractFloat}     

#Set of constants used in the saturation humidity calculations
#Where do these come from? What do they represent?

e0 ::NF= 6.108e-3
c1 ::NF= 17.269
c2 ::NF= 21.875
t0 ::NF= 273.16 #Paxton/Chantry have a zero_c term here? 
t1 ::NF= 35.86 
t2 ::NF= 7.66 

c3 ::NF= 622.0
c4 ::NF=0.378 

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