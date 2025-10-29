using SpeedyWeather, Adapt 

@kwdef struct SimpleAlbedo{NF<:Number} <: SpeedyWeather.AbstractParameterization end
    land_albedo::NF=0.3
    seaice_albedo::NF=0.5
    ocean_albedo::NF=0.1
end 

Adapt.@adapt_structure SimpleAlbedo

SimpleAlbedo(SG::SpectralGrid; kwargs...) = SimpleAlbedo(; kwargs...)

SpeedyWeather.initialize!(::SimpleAlbedo, model::PrimitiveEquation) = nothing 

SpeedyWeather.variables(::SimpleAlbedo) = (
    DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo", units="1"),
)

function SpeedyWeather.parameterization!(ij, diagn::DiagnosticVariables, progn::PrognosticVariables, albedo::SimpleAlbedo, model_parameters) 
    
    (; land_sea_mask) = model_parameters
    (; sea_ice_concentration) = progn.ocean
    (; land_albedo, seaice_albedo, ocean_albedo) = albedo
    
    if land_sea_mask[ij] > 0.95 # if mostly land 
        diagn.physics.albedo[ij] = land_albedo
    else # if ocean
        diagn.physics.albedo[ij] = ocean_albedo + sea_ice_concentration[ij] * (seaice_albedo - ocean_albedo)
    end
end 

spectral_grid = SpectralGrid(trunc=31, nlayers=8)
model = PrimitiveWetModel(spectral_grid; albedo=SimpleAlbedo())
simulation = initialize!(model)  

