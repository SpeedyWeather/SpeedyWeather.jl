abstract type AbstractDrag <: AbstractModelComponent end

## NO DRAG
export NoDrag
struct NoDrag <: AbstractDrag end
NoDrag(SG::SpectralGrid) = NoDrag()
initialize!(::NoDrag, ::AbstractModel) = nothing

function drag!(     diagn::DiagnosticVariables,
                    progn::PrognosticVariables,
                    drag::NoDrag,
                    model::AbstractModel,
                    lf::Integer)
    return nothing
end

# Quadratic drag
export QuadraticDrag
Base.@kwdef mutable struct QuadraticDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1]"
    c_D::NF = 1e-5

    "drag coefficient divided by layer thickness H, scaled with radius R [1]"
    c::Base.RefValue{NF} = Ref(zero(NF))
end

QuadraticDrag(SG::SpectralGrid; kwargs...) = QuadraticDrag{SG.NF}(; kwargs...)

function initialize!(   drag::QuadraticDrag,
                        model::AbstractModel)
    # c = c_D / H * R
    drag.c[] = drag.c_D / model.atmosphere.layer_thickness * model.geometry.radius
end

# function barrier
function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::QuadraticDrag,
    model::AbstractModel,
    lf::Integer,
)
    drag!(diagn, drag)
end

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u, v)|*(u, v)

with c_D the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`.
c_D and layer thickness `H` are precomputed in initialize!(::QuadraticDrag, ::AbstractModel)
and scaled by the radius as are the momentum equations."""
function drag!(     
    diagn::DiagnosticVariables,
    drag::QuadraticDrag,
)
    u = diagn.grid.u_grid
    v = diagn.grid.v_grid

    Fu = diagn.tendencies.u_tend_grid
    Fv = diagn.tendencies.v_tend_grid

    # total drag coefficient with radius scaling and /layer_thickness
    c = drag.c[]

    k = diagn.nlayers   # only apply to surface layer 
    @inbounds for ij in eachgridpoint(u, v, Fu, Fv)
        speed = sqrt(u[ij, k]^2 + v[ij, k]^2)
        Fu[ij, k] -= c*speed*u[ij, k]     # -= as the tendencies already contain forcing
        Fv[ij, k] -= c*speed*v[ij, k]
    end
end