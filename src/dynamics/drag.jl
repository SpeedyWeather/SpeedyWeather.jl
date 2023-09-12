function Base.show(io::IO,F::AbstractDrag)
    println(io,"$(typeof(F)) <: AbstractDrag")
    keys = propertynames(F)
    print_fields(io,F,keys)
end

## NO DRAG
struct NoDrag{NF} <: AbstractDrag{NF} end
NoDrag(SG::SpectralGrid) = NoDrag{SG.NF}()

function initialize!(   drag::NoDrag,
                        model::ModelSetup)
    return nothing
end

function drag!(     diagn::DiagnosticVariablesLayer,
                    drag::NoDrag,
                    time::DateTime,
                    model::ModelSetup)
    return nothing
end

# Quadratic drag
Base.@kwdef struct QuadraticDrag{NF} <: AbstractDrag{NF}
    "drag coefficient [1]"
    c_D::Float64 = 1e-5
end

QuadraticDrag(SG::SpectralGrid;kwargs...) = QuadraticDrag{SG.NF}(;kwargs...)

function initialize!(   drag::QuadraticDrag,
                        model::ModelSetup)
    return nothing
end

# function barrier
function drag!(     diagn::DiagnosticVariablesLayer,
                    drag::QuadraticDrag,
                    time::DateTime,
                    model::ModelSetup)
    drag!(diagn,drag,model.constants)
end

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u,v)|*(u,v)

with c_D the non-dimensional drag coefficient
as defined in `drag::QuadraticDrag`. Layer thickness `H`
is taken from the `Atmosphere` via `DynamicsConstants`,
as defined for the `ShallowWaterModel`."""
function drag!(     diagn::DiagnosticVariablesLayer,
                    drag::QuadraticDrag{NF},
                    C::DynamicsConstants) where NF
    
    u = diagn.grid_variables.u_grid
    v = diagn.grid_variables.v_grid

    Fu = diagn.tendencies.u_tend_grid
    Fv = diagn.tendencies.v_tend_grid

    # total drag coefficient with radius scaling and /layer_thickness
    c = convert(NF,drag.c_D/C.layer_thickness*C.radius)

    @inbounds for ij in eachgridpoint(u,v,Fu,Fv)
        speed = sqrt(u[ij]^2 + v[ij]^2)
        Fu[ij] -= c*speed*u[ij]     # -= as the tendencies already contain forcing
        Fv[ij] -= c*speed*v[ij]
    end
end