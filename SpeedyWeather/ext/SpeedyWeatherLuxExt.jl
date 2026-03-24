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
    layer_map = [
        "embed_layer" => :layer_1,
        "layer_1" => :layer_2,
        "layer_2" => :layer_4,
        "layer_3" => :layer_6,
        "output_layer" => :layer_7,
    ]
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

Base.@propagate_inbounds function surface_roughness_land(ij, diagn, progn, scheme::LearnedSurfaceRoughness)
    vₕ = diagn.physics.land.vegetation_high[ij]
    vₗ = diagn.physics.land.vegetation_low[ij]
    vᵦ = 1 - vₕ - vₗ  # bare soil
    g = diagn.grid.geopotential[ij, end]
    sd = progn.land.snow_depth[ij]
    soil_moisture = progn.land.soil_moisture[ij, begin]  # currently top layer
    soil_temperature = progn.land.soil_temperature[ij, end]  # currently bottom layer

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

Base.@propagate_inbounds function surface_roughness_ocean(ij, diagn, progn, scheme::LearnedSurfaceRoughness)
    surface = diagn.nlayers
    ℵ = progn.ocean.sea_ice_concentration[ij]
    Uₛ = diagn.grid.u_grid[ij, surface]
    Vₛ = diagn.grid.v_grid[ij, surface]
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

Base.@propagate_inbounds function SpeedyWeather.surface_roughness!(ij, diagn, progn, scheme::LearnedSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.mask[ij]

    # Compute separate ocean and land surface roughness
    diagn.physics.land.surface_roughness[ij] = land_fraction > 0 ? surface_roughness_land(ij, diagn, progn, scheme) : zero(land_fraction)
    diagn.physics.ocean.surface_roughness[ij] = land_fraction < 1 ? surface_roughness_ocean(ij, diagn, progn, scheme) : zero(land_fraction)

    # Blend the two via arithmetic average
    diagn.physics.surface_roughness[ij] = land_fraction * diagn.physics.land.surface_roughness[ij] + (1 - land_fraction) * diagn.physics.ocean.surface_roughness[ij]
    return nothing
end

# --- Albedo ---
function ResidualBlock(dim::Int, dropout_rate::Float64 = 0.1)
    inner_layers = Lux.Chain(
        Lux.Dense(dim => dim, Lux.gelu),
        Lux.LayerNorm((dim,)),
        Lux.Dropout(dropout_rate),
        Lux.Dense(dim => dim, Lux.gelu),
        Lux.LayerNorm((dim,))
    )
    # Lux.SkipConnection automatically handles the `out + identity` operation
    return Lux.SkipConnection(inner_layers, +)
end

function build_head(in_dim::Int, hidden_dim::Int)
    return Lux.Chain(
        Lux.Dense(in_dim => hidden_dim, Lux.gelu),
        Lux.LayerNorm((hidden_dim,)),
        Lux.Dense(hidden_dim => 3)
    )
end

function BRDFNet(; input_dim::Int = 8, shared_dim::Int = 256, head_dim::Int = 64, dropout_rate::Float64 = 0.1)
    # 1. Trunk (Embedding + ResBlocks)
    trunk = Lux.Chain(
        Lux.Dense(input_dim => shared_dim, Lux.gelu),
        Lux.LayerNorm((shared_dim,)),
        ResidualBlock(shared_dim, dropout_rate),
        ResidualBlock(shared_dim, dropout_rate)
    )

    # 2. Parallel Heads (vis, nir, sw)
    # BranchLayer passes the same input (from the trunk) to all three sub-chains
    heads = Lux.BranchLayer(
        build_head(shared_dim, head_dim),
        build_head(shared_dim, head_dim),
        build_head(shared_dim, head_dim)
    )

    # 3. Combiner
    # In Julia/Lux, the feature dimension is dim=1 (column-major: features x batch_size).
    # `vcat` correctly mirrors PyTorch's `torch.cat([...], dim=1)`
    combiner = Lux.WrappedFunction(x -> vcat(x...))

    # Chain them all together
    return Lux.Chain(trunk, heads, combiner)
end

function SpeedyWeather.LearnedBRDF(
        SG::SpectralGrid;
        snow_cover = SaturatingSnowCover(),
        input_dim::Int = 8,
        shared_dim::Int = 256,
        head_dim::Int = 64,
        dropout_rate::Float64 = 0.1,
        kwargs...
    )

    # Set up Lux NN
    brdf_nn = BRDFNet(; input_dim, shared_dim, head_dim, dropout_rate)

    rng = Random.default_rng()
    brdf_params, rand_states = Lux.setup(rng, brdf_nn)

    # Freeze dropout for inference by default
    brdf_states = Lux.testmode(rand_states)

    return LearnedBRDF{SG.NF, typeof(brdf_nn), typeof(brdf_params), typeof(brdf_states), typeof(snow_cover)}(;
        brdf_nn = brdf_nn,
        brdf_params = brdf_params,
        brdf_states = brdf_states,
        snow_cover = snow_cover,
        kwargs...
    )
end

function load_brdf_parameters!(params, weights::Dict{String, Array{Float32}})
    # Define a mapping from PyTorch keys to the corresponding Lux parameter nested structure.
    # Format: "pytorch_layer_prefix" => (Lux_layer_reference, is_layernorm)

    mapping = [
        # 1. Trunk Embed (Maps to params.layer_1)
        "trunk_embed.0" => (params.layer_1.layer_1, false),
        "trunk_embed.2" => (params.layer_1.layer_2, true),

        # 2. ResBlock 1 (Maps to params.layer_1.layer_3)
        "trunk_res_blocks.0.linear1" => (params.layer_1.layer_3.layer_1, false),
        "trunk_res_blocks.0.norm1" => (params.layer_1.layer_3.layer_2, true),
        "trunk_res_blocks.0.linear2" => (params.layer_1.layer_3.layer_4, false),
        "trunk_res_blocks.0.norm2" => (params.layer_1.layer_3.layer_5, true),

        # 3. ResBlock 2 (Maps to params.layer_1.layer_4)
        "trunk_res_blocks.1.linear1" => (params.layer_1.layer_4.layer_1, false),
        "trunk_res_blocks.1.norm1" => (params.layer_1.layer_4.layer_2, true),
        "trunk_res_blocks.1.linear2" => (params.layer_1.layer_4.layer_4, false),
        "trunk_res_blocks.1.norm2" => (params.layer_1.layer_4.layer_5, true),

        # 4. Heads (Maps to params.layer_2. BranchLayer stores them as layer_1, layer_2, layer_3)
        "vis_head.0" => (params.layer_2.layer_1.layer_1, false),
        "vis_head.2" => (params.layer_2.layer_1.layer_2, true),
        "vis_head.3" => (params.layer_2.layer_1.layer_3, false),

        "nir_head.0" => (params.layer_2.layer_2.layer_1, false),
        "nir_head.2" => (params.layer_2.layer_2.layer_2, true),
        "nir_head.3" => (params.layer_2.layer_2.layer_3, false),

        "sw_head.0" => (params.layer_2.layer_3.layer_1, false),
        "sw_head.2" => (params.layer_2.layer_3.layer_2, true),
        "sw_head.3" => (params.layer_2.layer_3.layer_3, false),
    ]

    for (py_prefix, (lux_layer, is_layernorm)) in mapping
        weight_key = py_prefix * ".weight"
        bias_key = py_prefix * ".bias"

        if haskey(weights, weight_key) && haskey(weights, bias_key)
            if is_layernorm
                # PyTorch LayerNorm uses 'weight', Lux uses 'scale'
                lux_layer.scale .= weights[weight_key]
            else
                # PyTorch and Lux Dense layers both expect (out_features, in_features)
                lux_layer.weight .= weights[weight_key]
            end
            lux_layer.bias .= weights[bias_key]
        else
            @warn "Weights or bias for $py_prefix not found in loaded PyTorch dictionary"
        end
    end

    return nothing
end

Adapt.@adapt_structure SpeedyWeather.LearnedBRDF
function SpeedyWeather.initialize!(brdf::LearnedBRDF, ::PrimitiveEquation)
    weights = RingGrids.get_asset(
        brdf.path,
        from_assets = brdf.from_assets,
        name = "brdf", # Adjust this if your asset name differs
        ArrayType = Array,
        FileFormat = Dict,
        version = brdf.version
    )

    load_brdf_parameters!(brdf.brdf_params, weights)

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
Base.@propagate_inbounds function calculate_black_sky_albedo(f_iso, f_vol, f_geo, mu0)
    # Computed directly to avoid heap allocations
    NF = typeof(f_iso)
    c_iso = NF(1.0) - NF(0.007574) * mu0 - NF(1.284909) * mu0^2
    c_vol = NF(-0.070987) * mu0 - NF(0.166314) * mu0^2
    c_geo = NF(0.307588) * mu0 + NF(0.04184) * mu0^2

    return (c_iso * f_iso) + (c_vol * f_vol) + (c_geo * f_geo)
end

""" Calculate BSA from SAL and WSA using fraction of direct radiation at the surface"""
Base.@propagate_inbounds function calculate_bsa_from_fraction(f_iso, f_vol, f_geo, mu0, fraction_direct)
    wsa = calculate_white_sky_albedo(f_iso, f_vol, f_geo)
    bsa = calculate_black_sky_albedo(f_iso, f_vol, f_geo, mu0)
    return (fraction_direct * bsa) + ((1 - fraction_direct) * wsa)
end

Base.@propagate_inbounds function brdf(ij, diagn, progn, scheme::LearnedBRDF)
    # Calculate snow cover
    snow_depth = progn.land.snow_depth[ij]
    snow_cover = scheme.snow_cover(snow_depth, scheme.snow_depth_scale)

    # Normalise inputs
    vegh = normalise(diagn.physics.land.vegetation_high[ij], scheme.input_means[1], scheme.input_stds[1]) 
    vegl = normalise(diagn.physics.land.vegetation_low[ij], scheme.input_means[2], scheme.input_stds[2])
    soil_moisture = normalise(progn.land.soil_moisture[ij, end], scheme.input_means[3], scheme.input_stds[3])
    soil_temperature = normalise(progn.land.soil_temperature[ij, end], scheme.input_means[4], scheme.input_stds[4])
    snow_cover = normalise(snow_cover, scheme.input_means[5], scheme.input_stds[5])
    geopotential = normalise(diagn.grid.geopotential[ij, end], scheme.input_means[6], scheme.input_stds[6])
    lai_hv = normalise(diagn.physics.land.lai_hv[ij], scheme.input_means[7], scheme.input_stds[7])
    lai_lv = normalise(diagn.physics.land.lai_lv[ij], scheme.input_means[8], scheme.input_stds[8])


    scheme.input_buffer[:] .= (
        vegh, vegl, soil_moisture, soil_temperature, snow_cover, geopotential, lai_hv, lai_lv
    )

    # TODO: sort out this reshaping stuff. Probably fixed by not using layer norm in the network.
    input_matrix = reshape(scheme.input_buffer, :, 1)
    prediction, _ = Lux.apply(scheme.brdf_nn, input_matrix, scheme.brdf_params, scheme.brdf_states)
    # vis_pred = prediction[1:3]
    # nir_pred = prediction[4:6]
    sw_pred_iso = prediction[7]
    sw_pred_vol = prediction[8]
    sw_pred_geo = prediction[9]

    # Unnorm outputs
    sw_pred_iso = denormalise(sw_pred_iso, scheme.output_means[7], scheme.output_stds[7])
    sw_pred_vol = denormalise(sw_pred_vol, scheme.output_means[8], scheme.output_stds[8])
    sw_pred_geo = denormalise(sw_pred_geo, scheme.output_means[9], scheme.output_stds[9])

    mu0 = diagn.physics.cos_zenith[ij]
    fraction_direct = diagn.physics.rad_fraction_direct[ij]

    return calculate_bsa_from_fraction(sw_pred_iso, sw_pred_vol, sw_pred_geo, mu0, fraction_direct)
end

Base.@propagate_inbounds function SpeedyWeather.albedo!(ij, diagn, progn, scheme::LearnedBRDF, model)
    land_fraction = model.land_sea_mask.mask[ij]
    
    # Don't run for fully ocean cells
    if land_fraction == 0
        diagn.physics.land.albedo[ij] = zero(land_fraction)
        return nothing
    end
    
    # Don't run for night-time cells
    mu = diagn.physics.cos_zenith[ij]
    if mu <= 0
        diagn.physics.land.albedo[ij] = one(mu) 
        return nothing
    end

    diagn.physics.land.albedo[ij] = clamp(brdf(ij, diagn, progn, scheme), 0, 1)
    return nothing
end

end
