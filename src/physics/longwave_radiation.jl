abstract type AbstractLongwave <: AbstractRadiation end

export NoLongwave
struct NoLongwave <: AbstractLongwave end
NoLongwave(SG::SpectralGrid) = NoLongwave()
initialize!(::NoLongwave,::PrimitiveEquation) = nothing

# function barrier for all AbstractLongwave
function longwave_radiation!(column::ColumnVariables,model::PrimitiveEquation)
    longwave_radiation!(column,model.longwave_radiation,model)
end

longwave_radiation!(::ColumnVariables,::NoLongwave,::PrimitiveEquation) = nothing

export UniformCooling
Base.@kwdef struct UniformCooling{NF} <: AbstractLongwave
    time_scale::Second = Hour(16)
    temp_min::NF = 207.5
    temp_stratosphere::NF = 200
    time_scale_stratosphere::Second = Day(5)
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
    (;temp, temp_tend) = column
    (;temp_min, temp_stratosphere) = scheme
    
    cooling = -inv(convert(NF,scheme.time_scale.value))
    τ⁻¹ = inv(scheme.time_scale_stratosphere.value)

    @inbounds for k in eachlayer(column)
        temp_tend[k] += temp[k] > temp_min ? cooling : (temp_stratosphere - temp[k])*τ⁻¹
    end
end