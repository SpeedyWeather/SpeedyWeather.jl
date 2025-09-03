abstract type AbstractSnow <: AbstractParameterization end

export SnowModel    # maybe change for a more concise name later
@kwdef mutable struct SnowModel{NF} <: AbstractSnow
    parameter1::NF = 0
    parameter2::NF = 0
end

# generator function
SnowModel(SG::SpectralGrid; kwargs...) = SnowModel{SG.NF}(; kwargs...)

# initialize component
initialize!(snow::SnowModel, model::PrimitiveEquation) = nothing

# set initial conditions for snow depth in initial conditions
function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    snow::SnowModel,
    model::PrimitiveEquation,
)
    set!(progn, model.geometry, snow_depth=0)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    snow::SnowModel,
    model::PrimitiveEquation,
)
    (; snow_depth) = progn.land
    (; mask) = model.land_sea_mask
    p1 = snow.parameter1
    # your code here updating snow_depth, e.g.

    launch!(architecture(snow_depth), LinearWorkOrder, size(snow_depth),
        land_snow_kernel!, snow_depth, mask, p1)
end

@kernel inbounds=true function land_snow_kernel!(
    snow_depth, mask, @Const(p1),
)
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land
        snow_depth[ij] = 0                  # dummy operation, replace with real logic
    end
end