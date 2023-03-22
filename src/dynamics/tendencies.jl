"""
    get_tendencies!(diagn,model)

Calculate all tendencies for the BarotropicModel."""
function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            model::BarotropicModel)
    
    # only (absolute) vorticity advection for the barotropic model
    vorticity_flux_divcurl!(diagn,model,curl=false)         # = -∇⋅(u(ζ+f),v(ζ+f))
end

"""
    get_tendencies!(diagn,surface,pres,time,model)

Calculate all tendencies for the ShallowWaterModel."""
function get_tendencies!(   diagn::DiagnosticVariablesLayer,
                            surface::SurfaceVariables,
                            pres::LowerTriangularMatrix,    # spectral pressure/η for geopotential
                            time::DateTime,                 # time to evaluate the tendencies at
                            model::ShallowWaterModel)       # struct containing all constants
    S,C = model.spectral_transform, model.constants

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    vorticity_flux_divcurl!(diagn,model,curl=true)  # = -∇⋅(u(ζ+f),v(ζ+f)), tendency for vorticity
                                                    # and ∇×(u(ζ+f),v(ζ+f)), tendency for divergence
    geopotential!(diagn,pres,C)                     # geopotential Φ = gη in the shallow water model
    bernoulli_potential!(diagn,S)                   # = -∇²(E+Φ), tendency for divergence
    volume_flux_divergence!(diagn,surface,model)    # = -∇⋅(uh,vh), tendency pressure

    # interface forcing
    @unpack interface_relaxation = model.parameters
    interface_relaxation && interface_relaxation!(pres,surface,time,model)
end

"""
    get_tendencies!(diagn,surface,pres,time,model)

Calculate all tendencies for the primitive equation model (wet or dry)."""
function get_tendencies!(   diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            time::DateTime,
                            model::PrimitiveEquation,
                            lf::Int=2               # leapfrog index for tendencies
                            )

    B, G, S = model.boundaries, model.geometry, model.spectral_transform
    @unpack surface = diagn

    # PARAMETERIZATIONS
    # parameterization_tendencies!(diagn,time,model)

    # DYNAMICS
    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step) 
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = model.parameters.implicit_α == 0 ? lf : 1

    pressure_gradient!(diagn,progn,lf,S)            # calculate ∇ln(pₛ)

    Threads.@threads for k in 1:diagn.nlev
        diagn_layer = diagn.layers[k]
        progn_layer = progn.layers[k]

        pressure_flux!(diagn_layer,surface)         # calculate (uₖ,vₖ)⋅∇ln(pₛ)

        # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
        linear_virtual_temperature!(diagn_layer,progn_layer,model,lf_implicit)
        temperature_anomaly!(diagn_layer,model)
    end

    geopotential!(diagn,B,G)                        # from ∂Φ/∂ln(pₛ) = -RTᵥ, used in bernoulli_potential!
    vertical_integration!(diagn,progn,lf_implicit,G)   # get ū,v̄,D̄ on grid; and and D̄ in spectral
    surface_pressure_tendency!(surface,model)       # ∂ln(pₛ)/∂t = -(ū,v̄)⋅∇ln(pₛ) - D̄

    Threads.@threads for layer in diagn.layers
        vertical_velocity!(layer,surface,model)     # calculate σ̇ for the vertical mass flux M = pₛσ̇
                                                    # add the RTₖlnpₛ term to geopotential
        linear_pressure_gradient!(layer,progn,model,lf_implicit)
    end

    vertical_advection!(diagn,model)                # use σ̇ for the vertical advection of u,v,T,q

    Threads.@threads for layer in diagn.layers
        vordiv_tendencies!(layer,surface,model)     # vorticity advection, pressure gradient term
        temperature_tendency!(layer,surface,model)  # hor. advection + adiabatic term
        humidity_tendency!(layer,model)             # horizontal advection of humidity (nothing for wetcore)
        bernoulli_potential!(layer,S)               # add -∇²(E+ϕ+RTₖlnpₛ) term to div tendency
    end
end