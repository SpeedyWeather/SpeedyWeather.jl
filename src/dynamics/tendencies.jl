"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the BarotropicModel."""
function dynamics_tendencies!(  diagn::DiagnosticVariablesLayer,
                                time::DateTime,
                                model::Barotropic)
    forcing!(diagn,model.forcing,time)      # = (Fᵤ, Fᵥ) forcing for u,v
    vorticity_flux!(diagn,model)            # = ∇×(v(ζ+f) + Fᵤ,v(ζ+f) + Fᵥ)
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the ShallowWaterModel."""
function dynamics_tendencies!(  diagn::DiagnosticVariablesLayer,
                                surface::SurfaceVariables,
                                pres::LowerTriangularMatrix,    # spectral pressure/η for geopotential
                                time::DateTime,                 # time to evaluate the tendencies at
                                model::ShallowWater)            # struct containing all constants

    S,C,G,O,F = model.spectral_transform, model.constants, model.geometry, model.orography, model.forcing

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    forcing!(diagn,F,time)                  # = (Fᵤ, Fᵥ, Fₙ) forcing for u,v,η
    vorticity_flux!(diagn,model)            # = ∇×(v(ζ+f) + Fᵤ,v(ζ+f) + Fᵥ), tendency for vorticity
                                            # = ∇⋅(v(ζ+f) + Fᵤ,v(ζ+f) + Fᵥ), tendency for divergence
    
    geopotential!(diagn,pres,C)                     # geopotential Φ = gη in the shallow water model
    bernoulli_potential!(diagn,S)                   # = -∇²(E+Φ), tendency for divergence
    volume_flux_divergence!(diagn,surface,O,C,G,S)  # = -∇⋅(uh,vh), tendency pressure
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the PrimitiveEquation model (wet or dry)."""
function dynamics_tendencies!(  diagn::DiagnosticVariables,
                                progn::PrognosticVariables,
                                model::PrimitiveEquation,
                                lf::Int=2)          # leapfrog index for tendencies

    O = model.orography
    G = model.geometry
    S = model.spectral_transform
    C = model.constants
    I = model.implicit
    (; surface ) = diagn

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step) 
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = model.implicit.α == 0 ? lf : 1

    pressure_gradient!(diagn,progn,lf,S)            # calculate ∇ln(pₛ)

    @floop for (diagn_layer, progn_layer) in zip(diagn.layers,progn.layers)
        pressure_flux!(diagn_layer,surface)         # calculate (uₖ,vₖ)⋅∇ln(pₛ)

        # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
        linear_virtual_temperature!(diagn_layer,progn_layer,model,lf_implicit)
        temperature_anomaly!(diagn_layer,I)           # temperature relative to profile
    end

    geopotential!(diagn,O,C)                        # from ∂Φ/∂ln(pₛ) = -RTᵥ, used in bernoulli_potential!
    vertical_integration!(diagn,progn,lf_implicit,G)   # get ū,v̄,D̄ on grid; and and D̄ in spectral
    surface_pressure_tendency!(surface,S)           # ∂ln(pₛ)/∂t = -(ū,v̄)⋅∇ln(pₛ) - D̄

    @floop for layer in diagn.layers
        vertical_velocity!(layer,surface,G)         # calculate σ̇ for the vertical mass flux M = pₛσ̇
                                                    # add the RTₖlnpₛ term to geopotential
        linear_pressure_gradient!(layer,progn.surface,lf_implicit,C,I)
    end                                             # wait all because vertical_velocity! needs to
                                                    # finish before vertical_advection!
    @floop for layer in diagn.layers
        vertical_advection!(layer,diagn,model)      # use σ̇ for the vertical advection of u,v,T,q

        vordiv_tendencies!(layer,surface,model)     # vorticity advection, pressure gradient term
        temperature_tendency!(layer,model)          # hor. advection + adiabatic term
        humidity_tendency!(layer,model)             # horizontal advection of humidity (nothing for wetcore)
        bernoulli_potential!(layer,S)               # add -∇²(E+ϕ+RTₖlnpₛ) term to div tendency
    end
end

"""
$(TYPEDSIGNATURES)
Set the tendencies in `diagn` to zero."""
function zero_tendencies!(diagn::DiagnosticVariables)
    for layer in diagn.layers
        fill!(layer.tendencies.u_tend_grid,0)
        fill!(layer.tendencies.v_tend_grid,0)
        fill!(layer.tendencies.temp_tend_grid,0)
        fill!(layer.tendencies.humid_tend_grid,0)
    end
    fill!(diagn.surface.pres_tend_grid,0)
    fill!(diagn.surface.pres_tend,0)
    return nothing
end
