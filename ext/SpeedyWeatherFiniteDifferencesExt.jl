module SpeedyWeatherFiniteDifferencesExt 

using SpeedyWeather 
import FiniteDifferences
import FiniteDifferences: to_vec 

# FiniteDifferences needs to be able to convert data structures to Vectors and back 
# This doesn't work out of the box with our data types, so we'll define those 
# conversions here.
function FiniteDifferences.to_vec(x::Grid) where Grid <: AbstractGridArray
    x_vec, from_vec = FiniteDifferences.to_vec(Array(x))

    function GridArray_from_vec(x_vec)
        return Grid(reshape(from_vec(x_vec), size(x)), x.nlat_half)
    end 

    return x_vec, GridArray_from_vec
end 

function FiniteDifferences.to_vec(x::LTA) where LTA <: LowerTriangularArray
    x_vec, from_vec = FiniteDifferences.to_vec(x.data)

    function LowerTriangularArray_from_vec(x_vec)
        return LowerTriangularArray(reshape(from_vec(x_vec), size(x)), x.m, x.n)
    end 

    return x_vec, LowerTriangularArray_from_vec
end

# Vector{Particle} needs an extra modification because an empty vector yields Any[] with to_vec for Particle (which isn't the case for number types)
function FiniteDifferences.to_vec(x::Vector{Particle{NF}}) where NF 
    if isempty(x) 
        return NF[], identity
    else # the else statement is the unmodified to_vec(::DenseVector)
        x_vecs_and_backs = map(to_vec, x)
        x_vecs, backs = first.(x_vecs_and_backs), last.(x_vecs_and_backs)
        function Vector_from_vec(x_vec)
            sz = cumsum(map(length, x_vecs))
            x_Vec = [backs[n](x_vec[sz[n] - length(x_vecs[n]) + 1:sz[n]]) for n in eachindex(x)]
            return oftype(x, x_Vec)
        end
        # handle empty x
        x_vec = isempty(x_vecs) ? eltype(eltype(x_vecs))[] : reduce(vcat, x_vecs)
        return x_vec, Vector_from_vec
    end 
end

# A version of the generic fallback from FiniteDifferences that excludes some of the fields 
# that we don't want to be varied for our big data structures 
function FiniteDifferences.to_vec(x::T) where {T <: Union{PrognosticVariables, PrognosticVariablesOcean, PrognosticVariablesLand, DiagnosticVariables, Tendencies, GridVariables, DynamicsVariables, PhysicsVariables, ParticleVariables}}

    val_vecs_and_backs = map(name -> to_vec(getfield(x, name)), included_fields(T))
    vals = first.(val_vecs_and_backs)
    backs = last.(val_vecs_and_backs)

    vals_excluded_pre = map(name -> getfield(x, name), excluded_fields_pre(T))
    vals_excluded_post = map(name -> getfield(x, name), excluded_fields_post(T))

    v, vals_from_vec = to_vec(vals)
    function structtype_from_vec(v::Vector{<:Real})
        val_vecs = vals_from_vec(v)
        values = map((b, v) -> b(v), backs, val_vecs)
        
        T(vals_excluded_pre..., values..., vals_excluded_post...)
    end
    return v, structtype_from_vec
end

included_fields(::Type{<:PrognosticVariables}) = (:vor, :div, :temp, :humid, :pres, :random_pattern, :ocean, :land, :tracers, :particles)
excluded_fields_pre(::Type{<:PrognosticVariables}) = (:trunc, :nlat_half, :nlayers, :nparticles)
excluded_fields_post(::Type{<:PrognosticVariables}) = (:scale, :clock)

included_fields(::Type{<:PrognosticVariablesOcean}) = (:sea_surface_temperature, :sea_ice_concentration)
excluded_fields_pre(::Type{<:PrognosticVariablesOcean}) = (:nlat_half, )
excluded_fields_post(::Type{<:PrognosticVariablesOcean}) = ()

included_fields(::Type{<:PrognosticVariablesLand}) = (:land_surface_temperature, :snow_depth, :soil_moisture_layer1, :soil_moisture_layer2)
excluded_fields_pre(::Type{<:PrognosticVariablesLand}) = (:nlat_half, )
excluded_fields_post(::Type{<:PrognosticVariablesLand}) = ()

included_fields(::Type{<:DiagnosticVariables}) = (:tendencies, :grid, :dynamics, :physics, :particles, :column, :temp_average)
excluded_fields_pre(::Type{<:DiagnosticVariables}) = (:trunc, :nlat_half, :nlayers, :nparticles)
excluded_fields_post(::Type{<:DiagnosticVariables}) = (:scale,)

included_fields(::Type{<:Tendencies}) = (:vor_tend, :div_tend, :temp_tend, :humid_tend, :u_tend, :v_tend, :pres_tend, :tracers_tend, :u_tend_grid, :v_tend_grid, :temp_tend_grid, :humid_tend_grid, :pres_tend_grid, :tracers_tend_grid)
excluded_fields_pre(::Type{<:Tendencies}) = (:trunc, :nlat_half, :nlayers)
excluded_fields_post(::Type{<:Tendencies}) = ()

included_fields(::Type{<:GridVariables}) = (:vor_grid, :div_grid, :temp_grid, :temp_virt_grid, :humid_grid, :u_grid, :v_grid, :pres_grid, :tracers_grid, :random_pattern, :temp_grid_prev, :humid_grid_prev, :u_grid_prev, :v_grid_prev, :pres_grid_prev, :tracers_grid_prev)
excluded_fields_pre(::Type{<:GridVariables}) = (:nlat_half, :nlayers)
excluded_fields_post(::Type{<:GridVariables}) = ()

included_fields(::Type{<:DynamicsVariables}) = (:a, :b, :a_grid, :b_grid, :a_2D, :b_2D, :a_2D_grid, :b_2D_grid, :uv∇lnp, :uv∇lnp_sum_above, :div_sum_above, :temp_virt, :geopot, :σ_tend, :∇lnp_x, :∇lnp_y, :u_mean_grid, :v_mean_grid, :div_mean_grid, :div_mean)
excluded_fields_pre(::Type{<:DynamicsVariables}) = (:trunc, :nlat_half, :nlayers)
excluded_fields_post(::Type{<:DynamicsVariables}) = ()

included_fields(::Type{<:PhysicsVariables}) = (:precip_large_scale, :precip_convection, :precip_rate_large_scale, :precip_rate_convection, :cloud_top, :soil_moisture_availability, :sensible_heat_flux, :evaporative_flux, :surface_shortwave_up, :surface_shortwave_down, :surface_longwave_up, :surface_longwave_down, :outgoing_shortwave_radiation, :outgoing_longwave_radiation, :cos_zenith)
excluded_fields_pre(::Type{<:PhysicsVariables}) = (:nlat_half,)
excluded_fields_post(::Type{<:PhysicsVariables}) = ()

included_fields(::Type{<:ParticleVariables}) = (:locations, :u, :v, :σ_tend)
excluded_fields_pre(::Type{<:ParticleVariables}) = (:nparticles, :nlat_half)
excluded_fields_post(::Type{<:ParticleVariables}) = (:interpolator, )

end 