function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            model::BarotropicModel,
                            )
    
    # only (planetary) vorticity advection for the barotropic model
    G,S = model.geometry, model.spectral_transform
    vorticity_flux_divcurl!(diagn,G,S,curl=false)           # = -∇⋅(u(ζ+f),v(ζ+f))
end

function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            surface::SurfaceVariables,
                            pres::LowerTriangularMatrix,    # spectral pressure/η for geopotential
                            time::DateTime,                 # time to evaluate the tendencies at
                            model::ShallowWaterModel,       # struct containing all constants
                            )

    G,S,C = model.geometry, model.spectral_transform, model.constants

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    vorticity_flux_divcurl!(diagn,G,S)              # = -∇⋅(u(ζ+f),v(ζ+f)), tendency for vorticity
                                                    # and ∇×(u(ζ+f),v(ζ+f)), tendency for divergence
    geopotential!(diagn,pres,C)                     # Φ = gη in the shallow water model
    bernoulli_potential!(diagn,G,S)                 # = -∇²(E+Φ), tendency for divergence
    volume_flux_divergence!(diagn,surface,model)    # = -∇⋅(uh,vh), tendency pressure

    # interface forcing
    @unpack interface_relaxation = model.parameters
    interface_relaxation && interface_relaxation!(pres,surface,time,model)
end

function get_tendencies!(   diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            time::DateTime,
                            model::PrimitiveEquation,
                            lf::Int=2                   # leapfrog index to evaluate tendencies on
                            )

    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform
    @unpack surface = diagn

    # PARAMETERIZATIONS
    # parameterization_tendencies!(diagn,time,model)
 
    # DYNAMICS
    pressure_gradients!(diagn,progn,lf,S)               # calculate ∇ln(pₛ)

    for layer in diagn.layers
        thickness_weighted_divergence!(layer,surface,G)    # calculate Δσₖ[(uₖ,vₖ)⋅∇ln(pₛ) + ∇⋅(uₖ,vₖ)]
    end

    geopotential!(diagn,B,G)                        # from ∂Φ/∂ln(pₛ) = -RTᵥ
    vertical_averages!(diagn,progn,lf,G)            # get ū,v̄,D̄ and others
    surface_pressure_tendency!(surface,model)       # ∂ln(pₛ)/∂t = -(ū,v̄)⋅∇ln(pₛ) - ∇⋅(ū,v̄)

    for layer in diagn.layers
        vertical_velocity!(layer,surface,model)     # calculate σ̇ for the vertical mass flux M = pₛσ̇
    end

    vertical_advection!(diagn,model)                # use σ̇ for the vertical advection of u,v,T,q

    for layer in diagn.layers
        vordiv_tendencies!(layer,surface,model)     # vorticity advection
        temperature_tendency!(layer,model)          # hor. advection + adiabatic term
        humidity_tendency!(layer,model)             # horizontal advection of humid
        bernoulli_potential!(layer,G,S)             # add -∇²(E+ϕ) term to div tendency
    end
end