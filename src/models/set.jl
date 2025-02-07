"""$(TYPEDSIGNATURES)
Sets a boundary condition fields for `model`. The input can be a function, `RingGrid`, `LowerTriangularMatrix`, 
or scalar as for other `set!` functions. If the keyword `add==true` the input is added to the exisiting 
field instead."""
function set!(
    model::AbstractModel;
    orography = nothing,
    land_sea_mask = nothing,
    albedo = nothing,
    #TODO add vegetation?
    kwargs...
)
    # orography also needs spectral transform and gravity for corresponding geopot_surf
    isnothing(orography) || set!(model.orography, orography, model.geometry, model.spectral_transform;
        gravity=model.planet.gravity, kwargs...)

    isnothing(land_sea_mask) || set!(model.land_sea_mask, land_sea_mask, model.geometry; kwargs...)
    isnothing(albedo) || set!(model.albedo, albedo, model.geometry; kwargs...)
    return nothing
end