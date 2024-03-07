abstract type AbstractDrag <: AbstractModelComponent end

## NO DRAG
export NoDrag
struct NoDrag <: AbstractDrag end
NoDrag(SG::SpectralGrid) = NoDrag()
initialize!(::NoDrag, ::ModelSetup) = nothing

function drag!(     diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLayer,
                    drag::NoDrag,
                    time::DateTime,
                    model::ModelSetup)
    return nothing
end

# Quadratic drag
export QuadraticDrag
Base.@kwdef mutable struct QuadraticDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1]"
    c_D::NF = 1e-5

    "drag coefficient divided by layer thickness H, scaled with radius R [1]"
    c::Base.RefValue{NF} = Ref(zero(NF))
end

QuadraticDrag(SG::SpectralGrid; kwargs...) = QuadraticDrag{SG.NF}(; kwargs...)

function initialize!(   drag::QuadraticDrag,
                        model::ModelSetup)
    # c = c_D / H * R
    drag.c[] = drag.c_D / model.atmosphere.layer_thickness * model.geometry.radius
end

# function barrier
function drag!(     diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLayer,
                    drag::QuadraticDrag,
                    time::DateTime,
                    model::ModelSetup)
    drag!(diagn, drag)
end

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u, v)|*(u, v)

with c_D the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`.
c_D and layer thickness `H` are precomputed in intialize!(::QuadraticDrag, ::ModelSetup)
and scaled by the radius as are the momentum equations."""
function drag!(     
    diagn::DiagnosticVariablesLayer,
    drag::QuadraticDrag{NF},
) where NF
    
    u = diagn.grid_variables.u_grid
    v = diagn.grid_variables.v_grid

    Fu = diagn.tendencies.u_tend_grid
    Fv = diagn.tendencies.v_tend_grid

    # total drag coefficient with radius scaling and /layer_thickness
    c = drag.c[]

    @inbounds for ij in eachgridpoint(u, v, Fu, Fv)
        speed = sqrt(u[ij]^2 + v[ij]^2)
        Fu[ij] -= c*speed*u[ij]     # -= as the tendencies already contain forcing
        Fv[ij] -= c*speed*v[ij]
    end
end