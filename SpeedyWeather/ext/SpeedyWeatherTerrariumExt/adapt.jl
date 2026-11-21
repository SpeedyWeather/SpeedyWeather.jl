# The Terrarium state lives in `vars.prognostic.land.terrarium` (see `TerrariumVars`) and is only used by the
# land component itself, which hands it to Terrarium directly. SpeedyWeather kernels never read it, but they
# receive (parts of) the `Variables`, which Adapt would otherwise convert including the whole Terrarium state.
# The state is therefore replaced by `nothing` when `Variables` are adapted, e.g. into kernel arguments.

# Replace the Terrarium state in the prognostic land namespace, if there is one; a no-op otherwise.
exclude_terrarium_state(prognostic::NamedTuple) =
    haskey(prognostic, :land) ? merge(prognostic, (land = exclude_terrarium_land_state(prognostic.land),)) : prognostic

# The land namespace is a `NamedTuple` too, so it needs a function of its own: a second
# `exclude_terrarium_state(::NamedTuple)` method would overwrite the one above.
exclude_terrarium_land_state(land::NamedTuple) = haskey(land, :terrarium) ? merge(land, (terrarium = nothing,)) : land
exclude_terrarium_land_state(land) = land

# More specific than SpeedyWeather's generic `Adapt.@adapt_structure Variables`, so it applies whenever this
# extension is loaded and otherwise adapts every group of variables as before.
function Adapt.adapt_structure(to, vars::Variables{<:NamedTuple})
    return Variables(
        Adapt.adapt(to, exclude_terrarium_state(vars.prognostic)),
        Adapt.adapt(to, vars.grid),
        Adapt.adapt(to, vars.tendencies),
        Adapt.adapt(to, vars.dynamics),
        Adapt.adapt(to, vars.parameterizations),
        Adapt.adapt(to, vars.particles),
        Adapt.adapt(to, vars.scratch),
        Adapt.adapt(to, vars.fused),
    )
end
