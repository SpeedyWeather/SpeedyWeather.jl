function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            M::BarotropicModel,
                            )
    
    # only (planetary) vorticity advection for the barotropic model
    G,S = M.geometry, M.spectral_transform
    vorticity_flux_divcurl!(diagn,G,S,curl=false)           # = -∇⋅(u(ζ+f),v(ζ+f))
end

function get_tendencies!(   pres::LowerTriangularMatrix,
                            diagn::DiagnosticVariablesLayer,
                            surface::SurfaceVariables,
                            time::DateTime,                 # time to evaluate the tendencies at
                            M::ShallowWaterModel,           # struct containing all constants
                            )

    G,S,B = M.geometry, M.spectral_transform, M.boundaries
    g,H₀  = M.constants.gravity, M.constants.layer_thickness

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    vorticity_flux_divcurl!(diagn,G,S)              # = -∇⋅(u(ζ+f),v(ζ+f)), tendency for vorticity
                                                    # and ∇×(u(ζ+f),v(ζ+f)), tendency for divergence
    bernoulli_potential!(diagn,surface,G,S,g)       # = -∇²(E+gη), tendency for divergence
    volume_flux_divergence!(diagn,surface,G,S,B,H₀) # = -∇⋅(uh,vh), tendency pressure

    # interface forcing
    @unpack interface_relaxation = M.parameters
    interface_relaxation && interface_relaxation!(pres,surface,time,M)
end

function get_tendencies!(   diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            time::DateTime,
                            model::PrimitiveEquationModel,
                            lf::Int=2                   # leapfrog index to evaluate tendencies on
                            )

    @unpack dry_core = model.parameters
    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform
    surf = diagn.surface
 
    pressure_gradients!(diagn,progn,lf,S)               # calculate ∇ln(pₛ)

    for layer in diagn.layers
        thickness_weighted_divergence!(layer,surf,G)    # calculate Δσₖ[(uₖ,vₖ)⋅∇ln(pₛ) + ∇⋅(uₖ,vₖ)]
    end

    geopotential!(diagn,B,G)                            # from ∂Φ/∂ln(pₛ) = -RTᵥ
    vertical_averages!(progn,diagn,lf,G)
    surface_pressure_tendency!(surf,model)              # 
    # vertical_velocity!(diagn,model)
    # vertical_advection!(diagn,model)

    for layer in diagn.layers
        vordiv_tendencies!(layer,diagn.surface,model)
        temperature_tendency!(layer,diagn.surface,model)
        dry_core || humidity_tendency!(layer,model)
        bernoulli_potential!(layer,G,S)
    end
end

function tendencies_not_zero(diagn::DiagnosticVariables)
    
    for k in 1:diagn.nlev
        any(diagn.layers[k].tendencies.vor_tend .!= 0) && println("ζ tendency not zero in layer $k")
        any(diagn.layers[k].tendencies.div_tend .!= 0) && println("D tendency not zero in layer $k")
        any(diagn.layers[k].tendencies.temp_tend .!= 0) && println("T tendency not zero in layer $k")
    end

    any(diagn.surface.pres_tend .!= 0) && println("ln(pₛ) tendency not zero")

    return nothing
end