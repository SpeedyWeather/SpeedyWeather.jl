function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            model::BarotropicModel,
                            )
    
    # only (absolute) vorticity advection for the barotropic model
    vorticity_flux_divcurl!(diagn,model,curl=false)         # = -∇⋅(u(ζ+f),v(ζ+f))
end

function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            surface::SurfaceVariables,
                            pres::LowerTriangularMatrix,    # spectral pressure/η for geopotential
                            time::DateTime,                 # time to evaluate the tendencies at
                            model::ShallowWaterModel,       # struct containing all constants
                            )

    S,C = model.spectral_transform, model.constants

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    vorticity_flux_divcurl!(diagn,model,curl=true)  # = -∇⋅(u(ζ+f),v(ζ+f)), tendency for vorticity
                                                    # and ∇×(u(ζ+f),v(ζ+f)), tendency for divergence
    geopotential!(diagn,pres,C)                     # Φ = gη in the shallow water model
    bernoulli_potential!(diagn,S)                   # = -∇²(E+Φ), tendency for divergence
    volume_flux_divergence!(diagn,surface,model)    # = -∇⋅(uh,vh), tendency pressure

    # interface forcing
    @unpack interface_relaxation = model.parameters
    interface_relaxation && interface_relaxation!(pres,surface,time,model)
end

function get_tendencies!(   diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            time::DateTime,
                            model::PrimitiveEquation,
                            lf::Int=2               # leapfrog index to evaluate tendencies on
                            )

    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform
    @unpack surface = diagn

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step) 
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = model.parameters.implicit_α == 0 ? lf : 1

    # PARAMETERIZATIONS
    # parameterization_tendencies!(diagn,time,model)
 
    # DYNAMICS
    pressure_gradients!(diagn,progn,lf,S)           # calculate ∇ln(pₛ)

    for (diagn_layer,progn_layer) in zip(diagn.layers,progn.layers)
        pressure_flux!(diagn_layer,surface)         # calculate (uₖ,vₖ)⋅∇ln(pₛ)

        # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
        linear_virtual_temperature!(diagn_layer,progn_layer,model,lf_implicit)
        temperature_anomaly!(diagn_layer,model)
    end

    geopotential!(diagn,B,G)                        # from ∂Φ/∂ln(pₛ) = -RTᵥ, used in bernoulli_potential!
    vertical_integration!(diagn,progn,lf_implicit,G)   # get ū,v̄,D̄ on grid; and
                                                    # and D̄ in spectral at prev time step i-1 via lf_linear
    surface_pressure_tendency!(surface,model)       # ∂ln(pₛ)/∂t = -(ū,v̄)⋅∇ln(pₛ) - ∇⋅(ū,v̄)

    for layer in diagn.layers
        vertical_velocity!(layer,surface,model)     # calculate σ̇ for the vertical mass flux M = pₛσ̇
        linear_pressure_gradient!(layer,progn,model,lf_implicit)
    end

    vertical_advection!(diagn,model)                # use σ̇ for the vertical advection of u,v,T,q

    for layer in diagn.layers
        vordiv_tendencies!(layer,surface,model)     # vorticity advection, pressure gradient term
        temperature_tendency!(layer,surface,model)  # hor. advection + adiabatic term
        humidity_tendency!(layer,model)             # horizontal advection of humid
        bernoulli_potential!(layer,S)               # add -∇²(E+ϕ+RTₖlnpₛ) term to div tendency
    end
end