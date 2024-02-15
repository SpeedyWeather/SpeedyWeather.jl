abstract type AbstractLongwave{NF} <: AbstractRadiation{NF} end

struct NoLongwave{NF} <: AbstractLongwave{NF} end
NoLongwave(SG::SpectralGrid) = NoLongwave{SG.NF}()
initialize!(::NoLongwave,::PrimitiveEquation) = nothing

# function barrier for all AbstractLongwave
function longwave_radiation!(column::ColumnVariables,model::PrimitiveEquation)
    longwave_radiation!(column,model.longwave_radiation,model)
end

longwave_radiation!(::ColumnVariables,::NoLongwave,::PrimitiveEquation) = nothing

Base.@kwdef struct UniformCooling{NF} <: AbstractLongwave{NF}
    time_scale::Second = Hour(16)
end

UniformCooling(SG::SpectralGrid;kwargs...) = UniformCooling{SG.NF}(;kwargs...)
initialize!(scheme::UniformCooling,model::PrimitiveEquation) = nothing

function longwave_radiation!(
    column::ColumnVariables,
    scheme::UniformCooling,
    model::PrimitiveEquation,
)
    longwave_radiation!(column,scheme)
end

function longwave_radiation!(
    column::ColumnVariables{NF},
    scheme::UniformCooling,
) where NF
    (;temp_tend) = column
    cooling = -inv(convert(NF,scheme.time_scale.value))
    @inbounds for k in eachlayer(column)
        temp_tend[k] += cooling
    end
end