module SpeedyWeatherLuxExt

using SpeedyWeather, Lux, NPZ
import Random, Adapt

const z₀_SATURATION_LIMIT = 3.92f-3 # observed roughness saturation, from Curcic (2020)

function RingGrids.load_asset(path::String, name::String, ArrayType::Type{<:Array}, FileFormat::Type{<:Dict}, FillValue::Float64)
    data = NPZ.npzread(path)
    return data
end

function SpeedyWeather.LearnedSurfaceRoughness(
        SG::SpectralGrid;
        kwargs...
    )
    # Set up Lux NN
    land_nn = Lux.Chain(
        Lux.Dense(13 => 32, Lux.leakyrelu),
        Lux.Dense(32 => 64, Lux.leakyrelu),
        Lux.Dropout(0.2),
        Lux.Dense(64 => 64, Lux.leakyrelu),
        Lux.Dropout(0.1),
        Lux.Dense(64 => 32, Lux.leakyrelu),
        Lux.Dense(32 => 1)
    )

    rng = Random.default_rng()
    land_params, rand_states = Lux.setup(rng, land_nn)
    land_states = Lux.testmode(rand_states)

    return LearnedSurfaceRoughness{SG.NF, typeof(land_nn), typeof(land_params), typeof(land_states)}(;
        land_nn = land_nn,
        land_params = land_params,
        land_states = land_states,
        kwargs...
    )
end

function load_land_parameters!(params, weights::Dict{String, Array{Float32}})
    layer_map = (
        "embed_layer" => :layer_1,
        "layer_1" => :layer_2,
        "layer_2" => :layer_4,
        "layer_3" => :layer_6,
        "output_layer" => :layer_7,
    )
    for (py_name, lux_sym) in layer_map
        # Check layer names match
        if hasproperty(params, lux_sym)
            lux_layer_params = getproperty(params, lux_sym)
            lux_layer_params.weight .= weights[py_name * ".weight"]
            lux_layer_params.bias .= weights[py_name * ".bias"]
        else
            @warn "Layer $lux_sym not found in Lux params"
        end
    end
    return
end

Adapt.@adapt_structure SpeedyWeather.LearnedSurfaceRoughness
function SpeedyWeather.initialize!(surface_roughness::LearnedSurfaceRoughness, ::PrimitiveEquation)
    weights = RingGrids.get_asset(
        surface_roughness.path,
        from_assets = surface_roughness.from_assets,
        name = "land_weights",
        ArrayType = Array,
        FileFormat = Dict,
        version = surface_roughness.version
    )
    load_land_parameters!(surface_roughness.land_params, weights)
    return nothing
end

# Symbolic regression-derived expression for ice-free ocean surface roughness
@inline function ice_free_roughness(x1, x2, x3)
    return (x3 * ((((((-0.096541196 * x3) - 0.30319092) / exp(x3)) - -0.7495353) / exp(x3)) - -0.7495353)) - (-0.040766064 * (((-0.6404533 + x2) * x2) - (x1 - 0.69219255)))
end

@inline function sea_ice_roughness(ℵ)
    return max(1.0f-3, 0.93f-3 * (1.0f0 - ℵ) + (6.05f-3 * exp(-17.0f0 * (ℵ - 0.5f0)^2)))
end

@inline function normalise(a, m, s)
    return (a - m) / s
end

@inline function denormalise(a, m, s)
    return (a * s) + m
end

Base.@propagate_inbounds function surface_roughness_land(ij, vars, scheme::LearnedSurfaceRoughness)
    vₕ = vars.parameterizations.land.vegetation_high[ij]
    vₗ = vars.parameterizations.land.vegetation_low[ij]
    vᵦ = 1 - vₕ - vₗ  # bare soil
    g = vars.grid.geopotential[ij, end]
    sd = vars.prognostic.land.snow_depth[ij]
    soil_moisture = vars.prognostic.land.soil_moisture[ij, begin]  # currently top layer
    soil_temperature = vars.prognostic.land.soil_temperature[ij, end]  # currently bottom layer

    # Modelled interactions
    i₁ = soil_temperature * soil_moisture
    i₂ = soil_temperature * sd
    i₃ = soil_temperature * vₕ
    i₄ = soil_moisture * vₕ
    i₅ = sd * vᵦ
    i₆ = soil_temperature * vᵦ

    # Normalise inputs
    vₕ = normalise(vₕ, scheme.land_input_means[2], scheme.land_input_stds[2])
    vₗ = normalise(vₗ, scheme.land_input_means[3], scheme.land_input_stds[3])
    vᵦ = normalise(vᵦ, scheme.land_input_means[1], scheme.land_input_stds[1])
    g = normalise(g, scheme.land_input_means[4], scheme.land_input_stds[4])
    sd = normalise(sd, scheme.land_input_means[5], scheme.land_input_stds[5])
    soil_moisture = normalise(soil_moisture, scheme.land_input_means[12], scheme.land_input_stds[12])
    soil_temperature = normalise(soil_temperature, scheme.land_input_means[7], scheme.land_input_stds[7])

    # Normalise interaction terms
    i₁ = normalise(i₁, scheme.land_input_means[11], scheme.land_input_stds[11])
    i₂ = normalise(i₂, scheme.land_input_means[10], scheme.land_input_stds[10])
    i₃ = normalise(i₃, scheme.land_input_means[9], scheme.land_input_stds[9])
    i₄ = normalise(i₄, scheme.land_input_means[13], scheme.land_input_stds[13])
    i₅ = normalise(i₅, scheme.land_input_means[6], scheme.land_input_stds[6])
    i₆ = normalise(i₆, scheme.land_input_means[8], scheme.land_input_stds[8])

    scheme.land_input_buffer[:] .= (vᵦ, vₕ, vₗ, g, sd, i₅, soil_temperature, i₆, i₃, i₂, i₁, soil_moisture, i₄)

    prediction, _ = Lux.apply(scheme.land_nn, scheme.land_input_buffer, scheme.land_params, scheme.land_states)
    log_surface_roughness = (prediction[1] * scheme.land_output_std) + scheme.land_output_mean
    surface_roughness = exp(log_surface_roughness)
    return surface_roughness
end

Base.@propagate_inbounds function surface_roughness_ocean(ij, vars::Variables, scheme::LearnedSurfaceRoughness)
    ℵ = vars.prognostic.ocean.sea_ice_concentration[ij]
    Uₛ = vars.grid.u[ij, end]  # surface wind speed
    Vₛ = vars.grid.v[ij, end]
    UVₛ = sqrt(Uₛ^2 + Vₛ^2)

    Uₛ = normalise(Uₛ, scheme.ocean_input_means[1], scheme.ocean_input_stds[1])
    Vₛ = normalise(Vₛ, scheme.ocean_input_means[2], scheme.ocean_input_stds[2])
    UVₛ = normalise(UVₛ, scheme.ocean_input_means[3], scheme.ocean_input_stds[3])

    log_ocean_roughness = ice_free_roughness(Uₛ, Vₛ, UVₛ) * scheme.ocean_output_std + scheme.ocean_output_mean
    ocean_roughness = min(z₀_SATURATION_LIMIT, exp(log_ocean_roughness))
    ℵ_roughness = sea_ice_roughness(ℵ)  # From IFS documentation, CY49R1

    surface_roughness = ℵ * ℵ_roughness + (1 - ℵ) * ocean_roughness
    return surface_roughness
end

Base.@propagate_inbounds function SpeedyWeather.surface_roughness!(ij, vars::Variables, scheme::LearnedSurfaceRoughness, land_sea_mask::LandSeaMask)
    land_fraction = land_sea_mask.mask[ij]

    # Compute separate ocean and land surface roughness
    vars.parameterizations.land.surface_roughness[ij] = land_fraction > 0 ? surface_roughness_land(ij, vars, scheme) : zero(land_fraction)
    vars.parameterizations.ocean.surface_roughness[ij] = land_fraction < 1 ? surface_roughness_ocean(ij, vars, scheme) : zero(land_fraction)

    # Blend the two via arithmetic average
    vars.parameterizations.surface_roughness[ij] = land_fraction * vars.parameterizations.land.surface_roughness[ij] + (1 - land_fraction) * vars.parameterizations.ocean.surface_roughness[ij]
    return nothing
end


function PreNormResidualBlock(dim::Int, activation::Function, dropout_rate::Float32 = 0.05f0, expansion_factor::Int = 2)
    hidden_dim = dim * expansion_factor

    inner_path = Lux.Chain(;
        norm1 = Lux.BatchNorm(dim),
        linear1 = Lux.Dense(dim => hidden_dim),
        act1 = Lux.WrappedFunction(activation),
        linear2 = Lux.Dense(hidden_dim => dim),
        dropout = Lux.Dropout(dropout_rate)
    )
    return Lux.SkipConnection(inner_path, +)
end

# BRDFNet
function BRDFNet(; input_dim = 10, shared_dim = 32, head_dim = 32, activation = gelu, dropout_rate = 0.1)

    trunk_embed = Lux.Chain(;
        linear = Lux.Dense(input_dim => shared_dim),
        norm = Lux.BatchNorm(shared_dim),
        act = Lux.WrappedFunction(activation)
    )

    trunk_res_blocks = Lux.Chain(;
        block1 = PreNormResidualBlock(shared_dim, activation, dropout_rate),
        block2 = PreNormResidualBlock(shared_dim, activation, dropout_rate),
        block3 = PreNormResidualBlock(shared_dim, activation, dropout_rate)
    )

    build_head(in_dim, hidden_dim, out_features) = Lux.Chain(;
        linear1 = Lux.Dense(in_dim => hidden_dim),
        norm = Lux.BatchNorm(hidden_dim),
        act = Lux.WrappedFunction(activation),
        linear2 = Lux.Dense(hidden_dim => out_features)
    )

    heads = Lux.Parallel(
        tuple;
        head1 = build_head(shared_dim, head_dim, 2),
        head2 = build_head(shared_dim, head_dim, 2),
        head3 = build_head(shared_dim, head_dim, 2)
    )
    return Lux.Chain(; trunk_embed, trunk_res_blocks, heads)
end

function load_brdfnet_parameters!(ps, st, weights::Dict{String, Any})
    # --- Trunk Embed ---
    ps.trunk_embed.linear.weight .= weights["trunk_embed.0.weight"]
    ps.trunk_embed.linear.bias .= weights["trunk_embed.0.bias"]
    ps.trunk_embed.norm.scale .= weights["trunk_embed.1.weight"]
    ps.trunk_embed.norm.bias .= weights["trunk_embed.1.bias"]
    st.trunk_embed.norm.running_mean .= weights["trunk_embed.1.running_mean"]
    st.trunk_embed.norm.running_var .= weights["trunk_embed.1.running_var"]

    # --- Trunk Res Blocks ---
    b1_ps = ps.trunk_res_blocks.block1
    b1_st = st.trunk_res_blocks.block1
    b2_ps = ps.trunk_res_blocks.block2
    b2_st = st.trunk_res_blocks.block2
    b3_ps = ps.trunk_res_blocks.block3
    b3_st = st.trunk_res_blocks.block3

    # Block 1
    b1_ps.norm1.scale .= weights["trunk_res_blocks.0.norm1.weight"]
    b1_ps.norm1.bias .= weights["trunk_res_blocks.0.norm1.bias"]
    b1_st.norm1.running_mean .= weights["trunk_res_blocks.0.norm1.running_mean"]
    b1_st.norm1.running_var .= weights["trunk_res_blocks.0.norm1.running_var"]
    b1_ps.linear1.weight .= weights["trunk_res_blocks.0.linear1.weight"]
    b1_ps.linear1.bias .= weights["trunk_res_blocks.0.linear1.bias"]
    b1_ps.linear2.weight .= weights["trunk_res_blocks.0.linear2.weight"]
    b1_ps.linear2.bias .= weights["trunk_res_blocks.0.linear2.bias"]

    # Block 2
    b2_ps.norm1.scale .= weights["trunk_res_blocks.1.norm1.weight"]
    b2_ps.norm1.bias .= weights["trunk_res_blocks.1.norm1.bias"]
    b2_st.norm1.running_mean .= weights["trunk_res_blocks.1.norm1.running_mean"]
    b2_st.norm1.running_var .= weights["trunk_res_blocks.1.norm1.running_var"]
    b2_ps.linear1.weight .= weights["trunk_res_blocks.1.linear1.weight"]
    b2_ps.linear1.bias .= weights["trunk_res_blocks.1.linear1.bias"]
    b2_ps.linear2.weight .= weights["trunk_res_blocks.1.linear2.weight"]
    b2_ps.linear2.bias .= weights["trunk_res_blocks.1.linear2.bias"]

    # Block 3
    b3_ps.norm1.scale .= weights["trunk_res_blocks.2.norm1.weight"]
    b3_ps.norm1.bias .= weights["trunk_res_blocks.2.norm1.bias"]
    b3_st.norm1.running_mean .= weights["trunk_res_blocks.2.norm1.running_mean"]
    b3_st.norm1.running_var .= weights["trunk_res_blocks.2.norm1.running_var"]
    b3_ps.linear1.weight .= weights["trunk_res_blocks.2.linear1.weight"]
    b3_ps.linear1.bias .= weights["trunk_res_blocks.2.linear1.bias"]
    b3_ps.linear2.weight .= weights["trunk_res_blocks.2.linear2.weight"]
    b3_ps.linear2.bias .= weights["trunk_res_blocks.2.linear2.bias"]

    # Iso head
    ps.heads.head1.linear1.weight .= weights["iso_head.0.weight"]
    ps.heads.head1.linear1.bias .= weights["iso_head.0.bias"]
    ps.heads.head1.norm.scale .= weights["iso_head.1.weight"]
    ps.heads.head1.norm.bias .= weights["iso_head.1.bias"]
    st.heads.head1.norm.running_mean .= weights["iso_head.1.running_mean"]
    st.heads.head1.norm.running_var .= weights["iso_head.1.running_var"]
    ps.heads.head1.linear2.weight .= weights["iso_head.3.weight"]
    ps.heads.head1.linear2.bias .= weights["iso_head.3.bias"]

    # Vol head
    ps.heads.head2.linear1.weight .= weights["vol_head.0.weight"]
    ps.heads.head2.linear1.bias .= weights["vol_head.0.bias"]
    ps.heads.head2.norm.scale .= weights["vol_head.1.weight"]
    ps.heads.head2.norm.bias .= weights["vol_head.1.bias"]
    st.heads.head2.norm.running_mean .= weights["vol_head.1.running_mean"]
    st.heads.head2.norm.running_var .= weights["vol_head.1.running_var"]
    ps.heads.head2.linear2.weight .= weights["vol_head.3.weight"]
    ps.heads.head2.linear2.bias .= weights["vol_head.3.bias"]

    # Geo head
    ps.heads.head3.linear1.weight .= weights["geo_head.0.weight"]
    ps.heads.head3.linear1.bias .= weights["geo_head.0.bias"]
    ps.heads.head3.norm.scale .= weights["geo_head.1.weight"]
    ps.heads.head3.norm.bias .= weights["geo_head.1.bias"]
    st.heads.head3.norm.running_mean .= weights["geo_head.1.running_mean"]
    st.heads.head3.norm.running_var .= weights["geo_head.1.running_var"]
    ps.heads.head3.linear2.weight .= weights["geo_head.3.weight"]
    ps.heads.head3.linear2.bias .= weights["geo_head.3.bias"]
    return nothing
end

function SpeedyWeather.LearnedLandAlbedo(
        SG::SpectralGrid;
        snow_cover = SaturatingSnowCover(),
        input_dim::Int = 10,
        shared_dim::Int = 32,
        head_dim::Int = 32,
        dropout_rate::Float32 = 0.0f0,
        kwargs...
    )

    # Set up Lux NN
    brdf_nn = BRDFNet(; input_dim, shared_dim, head_dim, dropout_rate)

    rng = Random.default_rng()
    brdf_params, rand_states = Lux.setup(rng, brdf_nn)
    brdf_states = Lux.testmode(rand_states)

    input_buffer = zeros(SG.NF, 10, SG.npoints)

    return LearnedLandAlbedo{SG.NF, typeof(brdf_nn), typeof(brdf_params), typeof(brdf_states), typeof(snow_cover)}(;
        brdf_nn = brdf_nn,
        brdf_params = brdf_params,
        brdf_states = brdf_states,
        snow_cover = snow_cover,
        input_buffer = input_buffer,
        kwargs...
    )
end

Adapt.@adapt_structure SpeedyWeather.LearnedLandAlbedo
function SpeedyWeather.initialize!(brdf::LearnedLandAlbedo, ::PrimitiveEquation)
    params = RingGrids.get_asset(
        brdf.path,
        from_assets = brdf.from_assets,
        name = "brdf", # Adjust this if your asset name differs
        ArrayType = Array,
        FileFormat = Dict,
        version = brdf.version
    )

    brdf.norm_means .= params["x_mean"]
    brdf.norm_stds .= params["x_std"]
    brdf.unnorm_means .= params["y_mean"]
    brdf.unnorm_stds .= params["y_std"]

    load_brdfnet_parameters!(brdf.brdf_params, brdf.brdf_states, params)
    return nothing
end

"""
    Calculates the SZA-independent White-Sky Albedo (WSA).
    Based on the Lucht et al. (2000) analytical integrals.
"""
Base.@propagate_inbounds function calculate_white_sky_albedo(f_iso, f_vol, f_geo)
    NF = typeof(f_iso)
    wsa_vol_const = NF(0.189184)
    wsa_geo_const = NF(-1.377622)
    return f_iso + (wsa_vol_const * f_vol) + (wsa_geo_const * f_geo)
end

"""
Calculate SAL using the Lucht et al. (2000) polynomial.
"""
Base.@propagate_inbounds function calculate_black_sky_albedo(f_iso, f_vol, f_geo, θ)
    NF = typeof(f_iso)
    c_iso = NF(1.0)
    c_vol = evalpoly(θ, (NF(-0.007574), NF(0.0), NF(-0.070987), NF(0.307588)))
    c_geo = evalpoly(θ, (NF(-1.284909), NF(0.0), NF(-0.166314), NF(0.04184)))
    
    return (c_iso * f_iso) + (c_vol * f_vol) + (c_geo * f_geo)
end

""" Calculate BSA from SAL and WSA using fraction of direct radiation at the surface"""
Base.@propagate_inbounds function calculate_bsa_from_fraction(f_iso, f_vol, f_geo, θ, fraction_direct)
    wsa = calculate_white_sky_albedo(f_iso, f_vol, f_geo)
    bsa = calculate_black_sky_albedo(f_iso, f_vol, f_geo, θ)
    return (fraction_direct * bsa) + ((1 - fraction_direct) * wsa)
end

const vis_weight, nir_weight = 0.5308f0, 0.4771f0

Base.@propagate_inbounds function brdf(ij, vars, scheme::LearnedLandAlbedo)
    # Calculate snow cover
    snow_depth = vars.prognostic.land.snow_depth[ij]
    snow_cover = scheme.snow_cover(snow_depth, scheme.snow_depth_scale) * 100

    # Normalise inputs
    vegh = normalise(vars.parameterizations.land.vegetation_high[ij], scheme.norm_means[1], scheme.norm_stds[1])
    vegl = normalise(vars.parameterizations.land.vegetation_low[ij], scheme.norm_means[2], scheme.norm_stds[2])
    soil_moisture1 = normalise(vars.prognostic.land.soil_moisture[ij, begin], scheme.norm_means[3], scheme.norm_stds[3])
    soil_temperature1 = normalise(vars.prognostic.land.soil_temperature[ij, begin], scheme.norm_means[4], scheme.norm_stds[4])
    soil_moisture2 = normalise(vars.prognostic.land.soil_moisture[ij, end], scheme.norm_means[5], scheme.norm_stds[5])
    soil_temperature2 = normalise(vars.prognostic.land.soil_temperature[ij, end], scheme.norm_means[6], scheme.norm_stds[6])
    geopotential = normalise(vars.grid.geopotential[ij, end], scheme.norm_means[7], scheme.norm_stds[7])
    lai_hv = normalise(vars.parameterizations.land.lai_vegetation_high[ij], scheme.norm_means[8], scheme.norm_stds[8])
    lai_lv = normalise(vars.parameterizations.land.lai_vegetation_low[ij], scheme.norm_means[9], scheme.norm_stds[9])
    snow_cover = normalise(snow_cover, scheme.norm_means[10], scheme.norm_stds[10])

    scheme.input_buffer[:] .= (
        vegh, vegl, soil_moisture1,
        soil_temperature1, soil_moisture2, soil_temperature2,
        geopotential, lai_hv, lai_lv, snow_cover,
    )

    input_matrix = reshape(scheme.input_buffer, :, 1)
    (out1, out2, out3), _ = Lux.apply(scheme.brdf_nn, input_matrix, scheme.brdf_params, scheme.brdf_states)

    # Unpack predictions
    vis_iso = out1[1]
    vis_vol = out1[2]
    vis_geo = out2[1]
    nir_iso = out2[2]
    nir_vol = out3[1]
    nir_geo = out3[2]

    # Unnorm outputs
    vis_iso = denormalise(vis_iso, scheme.unnorm_means[1], scheme.unnorm_stds[1])
    vis_vol = denormalise(vis_vol, scheme.unnorm_means[2], scheme.unnorm_stds[2])
    vis_geo = denormalise(vis_geo, scheme.unnorm_means[3], scheme.unnorm_stds[3])
    nir_iso = denormalise(nir_iso, scheme.unnorm_means[4], scheme.unnorm_stds[4])
    nir_vol = denormalise(nir_vol, scheme.unnorm_means[5], scheme.unnorm_stds[5])
    nir_geo = denormalise(nir_geo, scheme.unnorm_means[6], scheme.unnorm_stds[6])

    # Do the above but more efficiently


    sw_pred_iso = (vis_iso * vis_weight) + (nir_iso * nir_weight)
    sw_pred_vol = (vis_vol * vis_weight) + (nir_vol * nir_weight)
    sw_pred_geo = (vis_geo * vis_weight) + (nir_geo * nir_weight)

    μ = vars.parameterizations.cos_zenith[ij]
    θ = acos(μ)
    fraction_direct = vars.parameterizations.direct_radiation_fraction[ij]

    return calculate_bsa_from_fraction(sw_pred_iso, sw_pred_vol, sw_pred_geo, θ, fraction_direct)
end

# Base.@propagate_inbounds function SpeedyWeather.albedo!(ij, vars, scheme::LearnedLandAlbedo, model)
#     land_fraction = model.land_sea_mask.mask[ij]

#     # Don't run for fully ocean cells
#     if land_fraction == 0
#         vars.parameterizations.land.albedo[ij] = zero(land_fraction)
#         return nothing
#     end

#     # Don't run for night-time cells
#     μ = vars.parameterizations.cos_zenith[ij]
#     if μ <= 0
#         return nothing
#     end

#     vars.parameterizations.land.albedo[ij] = clamp(brdf(ij, vars, scheme), 0, 1)
#     return nothing
# end

# The global, ij independent albedo!
Base.@propagate_inbounds function SpeedyWeather.albedo!(vars::Variables, scheme::LearnedLandAlbedo, model)
    NF = model.spectral_grid.NF

    vegh = vars.parameterizations.land.vegetation_high
    vegl = vars.parameterizations.land.vegetation_low
    sm1 = vars.prognostic.land.soil_moisture[:, begin]
    st1 = vars.prognostic.land.soil_temperature[:, begin]
    sm2 = vars.prognostic.land.soil_moisture[:, end]
    st2 = vars.prognostic.land.soil_temperature[:, end]
    geopotential = vars.grid.geopotential[:, end]
    lai_hv = vars.parameterizations.land.lai_vegetation_high
    lai_lv = vars.parameterizations.land.lai_vegetation_low
    snow_depth = vars.prognostic.land.snow_depth

    snow_cover = scheme.snow_cover.(snow_depth, scheme.snow_depth_scale) * 100

    X = scheme.input_buffer

    @views @. X[1, :] = normalise(vegh, scheme.norm_means[1], scheme.norm_stds[1])
    @views @. X[2, :] = normalise(vegl, scheme.norm_means[2], scheme.norm_stds[2])
    @views @. X[3, :] = normalise(sm1, scheme.norm_means[3], scheme.norm_stds[3])
    @views @. X[4, :] = normalise(st1, scheme.norm_means[4], scheme.norm_stds[4])
    @views @. X[5, :] = normalise(sm2, scheme.norm_means[5], scheme.norm_stds[5])
    @views @. X[6, :] = normalise(st2, scheme.norm_means[6], scheme.norm_stds[6])
    @views @. X[7, :] = normalise(geopotential, scheme.norm_means[7], scheme.norm_stds[7])
    @views @. X[8, :] = normalise(lai_hv, scheme.norm_means[8], scheme.norm_stds[8])
    @views @. X[9, :] = normalise(lai_lv, scheme.norm_means[9], scheme.norm_stds[9])
    @views @. X[10, :] = normalise(snow_cover, scheme.norm_means[10], scheme.norm_stds[10])

    (out1, out2, out3), _ = Lux.apply(scheme.brdf_nn, X, scheme.brdf_params, scheme.brdf_states)

    v_iso, v_vol = @view(out1[1, :]), @view(out1[2, :])
    n_iso, n_vol = @view(out2[2, :]), @view(out3[1, :])
    v_geo, n_geo = @view(out2[1, :]), @view(out3[2, :])

    μ = vars.parameterizations.cos_zenith
    fraction_direct = vars.parameterizations.direct_radiation_fraction
    albedo_out = vars.parameterizations.land.albedo

    # Pre-allocate arrays
    sw_pred_iso = similar(albedo_out)
    sw_pred_vol = similar(albedo_out)
    sw_pred_geo = similar(albedo_out)
    θ = similar(μ)
    α = similar(albedo_out)

    @. albedo_out = begin
        v_iso = denormalise(v_iso, scheme.unnorm_means[1], scheme.unnorm_stds[1])
        v_vol = denormalise(v_vol, scheme.unnorm_means[2], scheme.unnorm_stds[2])
        v_geo = denormalise(v_geo, scheme.unnorm_means[3], scheme.unnorm_stds[3])
        n_iso = denormalise(n_iso, scheme.unnorm_means[4], scheme.unnorm_stds[4])
        n_vol = denormalise(n_vol, scheme.unnorm_means[5], scheme.unnorm_stds[5])
        n_geo = denormalise(n_geo, scheme.unnorm_means[6], scheme.unnorm_stds[6])

        sw_pred_iso = (v_iso * vis_weight) + (n_iso * nir_weight)
        sw_pred_vol = (v_vol * vis_weight) + (n_vol * nir_weight)
        sw_pred_geo = (v_geo * vis_weight) + (n_geo * nir_weight)

        θ = acos(μ)

        α = calculate_bsa_from_fraction(sw_pred_iso, sw_pred_vol, sw_pred_geo, θ, fraction_direct)
        clamp(α, zero(NF), one(NF))
    end

    return nothing
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    vars::SpeedyWeather.Variables, 
    albedo::OceanLandAlbedo{Ocean, <:LearnedLandAlbedo}, 
    model::SpeedyWeather.PrimitiveEquation
) where {Ocean}
    
    for ij in eachindex(model.land_sea_mask.mask) 
        SpeedyWeather.albedo!(ij, vars, albedo.ocean, model)
    end

    SpeedyWeather.albedo!(vars, albedo.land, model)

    # Blend ocean and land albedo
    mask = model.land_sea_mask.mask
    @. vars.parameterizations.albedo = mask * vars.parameterizations.land.albedo + 
                                       (1 - mask) * vars.parameterizations.ocean.albedo

    return nothing
end

Base.@propagate_inbounds function SpeedyWeather.albedo!(
    ij::Int, 
    vars::SpeedyWeather.Variables, 
    scheme::LearnedLandAlbedo, 
    model::SpeedyWeather.PrimitiveEquation
)
    return nothing
end

end
