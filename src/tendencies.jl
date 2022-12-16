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
 
    # # parametrization_tendencies!(diagn,time,model)
    # for layer in diagn.layers
    #     for fieldname in fieldnames(typeof(layer.tendencies))
    #         field = getfield(layer.tendencies,fieldname)
    #         fill!(field,0)
    #     end
    # end

    geopotential!(diagn,B,G)
    vertical_averages!(progn,diagn,lf,G)
    surface_pressure_tendency!(progn,diagn,lf,model)
    # vertical_velocity!(diagn,model)
    # vertical_advection!(diagn,model)

    for layer in diagn.layers
        vordiv_tendencies!(layer,diagn.surface,model)
        temperature_tendency!(layer,diagn.surface,model)
        dry_core || humidity_tendency!(layer,model)
        bernoulli_potential!(layer,G,S)
    end
end

"""
    add_tendencies!(tend::LowerTriangularMatrix{NF},    # tendency to accumulate into
                    term1::LowerTriangularMatrix{NF},   # with term1
                    term2::LowerTriangularMatrix{NF}    # and term2
                    ) where NF                          # number format real or complex

Accumulates three `LowerTriangularMatrix`s element-wise into the first.

    tend += term1 + term2."""
function add_tendencies!(   tend::LowerTriangularMatrix{NF},    # tendency to accumulate into
                            term1::LowerTriangularMatrix{NF},   # with term1
                            term2::LowerTriangularMatrix{NF}    # and term2
                            ) where NF                          # number format real or complex

    @inbounds for lm in eachharmonic(tend,term1,term2)
        tend[lm] += (term1[lm] + term2[lm])
    end
end

"""
    add_tendencies!(tend::LowerTriangularMatrix{NF},    # tendency to accumulate into
                    term::LowerTriangularMatrix{NF},    # with term
                    ) where NF                          # number format real or complex

Accumulates two `LowerTriangularMatrix`s element-wise into the first.

    tend += term."""
function add_tendencies!(   tend::LowerTriangularMatrix{NF},    # tendency to accumulate into
                            term::LowerTriangularMatrix{NF}     # with term
                            ) where NF                          # number format real or complex

    @inbounds for lm in eachharmonic(tend,term)
        tend[lm] += term[lm]
    end
end