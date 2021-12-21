struct Tendencies{NF<:AbstractFloat}
    vor_tend     ::Array{Complex{NF},3}      # Vorticity of horizontal wind field
    div_tend     ::Array{Complex{NF},3}      # Divergence of horizontal wind field
    Tabs_tend    ::Array{Complex{NF},3}      # Absolute temperature [K]
    logp0_tend   ::Array{Complex{NF},2}      # Log of surface pressure [log(Pa)]
    humid_tend   ::Array{Complex{NF},3}      # Specific humidity [g/kg]
end

struct GridVariables{NF<:AbstractFloat}
    u       ::Array{NF,3}       # Zonal velocity [m/s]
    v       ::Array{NF,3}       # Meridional velocity [m/s]
    Tabs    ::Array{NF,3}       # Absolute temperature [K]
    logp0   ::Array{NF,2}       # Log of surface pressure [log(Pa)]
    geopot  ::Array{NF,3}       # 
    humid   ::Array{NF,3}      # Specific humidity [g/kg]
end

"""Struct holding the diagnostic variables."""
struct DiagnosticVariables{NF<:AbstractFloat}
    tendencies  ::Tendencies{NF}
    gridvars    ::GridVariables{NF}
end