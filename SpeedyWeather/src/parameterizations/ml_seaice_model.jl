"""
Native Julia re-implementation of the trained `RelaxationMoEClassifier` (PyTorch) forward
pass, for use inside a SpeedyWeather.jl sea-ice component without a Julia<->Python bridge.

Weights are exported once from the trained checkpoint (see
`ml_seaice_weights/` + the export script referenced in project memory) as raw
row-major Float32 binaries plus a `manifest.txt` describing each tensor's shape.

Model recap (must stay in lockstep with the Python `RelaxationMoEClassifier`):
    h = backbone(x)                                   # 4x(Linear+GELU), width 64
    diffs = hist[2:end] .- hist[1:end-1] (append sic_raw as last point)
    mom = momentum_mlp(diffs)                          # Linear+GELU, width 16
    combined = [h; mom]                                # 80-d
    class_logits = class_head(combined)                # 2 (inc, dec)
    p = softmax(class_logits)
    sic_eq_inc = sic + (1-sic)*sigmoid(inc_head(combined))
    sic_eq_dec = sic*sigmoid(dec_head(combined))
    tau_inc = tau_floor + softplus(tau_inc_param); tau_dec similarly
    tendency = p_inc*(sic_eq_inc-sic)*(1-exp(-dt/tau_inc)) + p_dec*(sic_eq_dec-sic)*(1-exp(-dt/tau_dec))
"""
module MLSeaIceModel

export MLWeights, load_ml_weights, ml_forward, ml_normalize_features
export MLWeightsChangeGate, load_ml_weights_changegate

# ---------------------------------------------------------------------------
# erf / gelu -- no SpecialFunctions dependency; Abramowitz & Stegun 7.1.26,
# max abs error ~1.5e-7, plenty accurate for this use.
# ---------------------------------------------------------------------------
@inline function erf_approx(x::Float32)
    a1 = 0.254829592f0; a2 = -0.284496736f0; a3 = 1.421413741f0
    a4 = -1.453152027f0; a5 = 1.061405429f0; p = 0.3275911f0
    s = sign(x)
    ax = abs(x)
    t = 1f0 / (1f0 + p * ax)
    y = 1f0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-ax * ax)
    return s * y
end

@inline gelu(x::Float32) = 0.5f0 * x * (1f0 + erf_approx(x / sqrt(2f0)))
@inline gelu(x::AbstractArray{Float32}) = gelu.(x)

@inline sigmoid(x::Float32) = 1f0 / (1f0 + exp(-x))
@inline softplus(x::Float32) = log1p(exp(-abs(x))) + max(x, 0f0)

# ---------------------------------------------------------------------------
# Weight loading
# ---------------------------------------------------------------------------
struct MLWeights
    # backbone: 4 Linear+GELU layers, stored as (in,out) so that y = W'x + b
    W::Vector{Matrix{Float32}}
    b::Vector{Vector{Float32}}
    # momentum mlp
    Wm::Matrix{Float32}
    bm::Vector{Float32}
    # heads
    Wc::Matrix{Float32}; bc::Vector{Float32}   # class_head
    Wi::Matrix{Float32}; bi::Vector{Float32}   # inc_head
    Wd::Matrix{Float32}; bd::Vector{Float32}   # dec_head
    tau_inc_param::Float32
    tau_dec_param::Float32
    tau_floor::Float32
    dt::Float32
    # normalization (6 phys vars: siconc,tos,radsum,doy_sin,sfcwind,tas order depends on feature_set;
    # here: mean/std over [siconc,tos,rsds,rls,sfcwind,tas] raw cols 0..5)
    norm_mean::Vector{Float32}
    norm_std::Vector{Float32}
    radsum_mean::Float32
    radsum_std::Float32
    cross_mean::Vector{Float32}
    cross_std::Vector{Float32}
end

function _read_bin(dir::String, fname::AbstractString, shape)
    path = joinpath(dir, fname)
    data = reinterpret(Float32, read(path))
    if shape == "0"
        return Float32(data[1])
    end
    dims = Tuple(parse(Int, s) for s in split(shape, ","))
    if length(dims) == 1
        return Vector{Float32}(data)
    else
        out_f, in_f = dims  # PyTorch (out,in) row-major
        # row-major (out,in) flat == column-major (in,out) reshape
        return reshape(Vector{Float32}(data), (in_f, out_f))
    end
end

function load_ml_weights(dir::String; tau_floor::Float32=1.0f0, dt::Float32=5.0f0)
    manifest = readlines(joinpath(dir, "manifest.txt"))
    tensors = Dict{String, Any}()
    for line in manifest
        isempty(strip(line)) && continue
        name, shape, fname = split(line, "\t")
        tensors[name] = _read_bin(dir, fname, shape)
    end

    W = [tensors["backbone.0.weight"], tensors["backbone.2.weight"],
         tensors["backbone.4.weight"], tensors["backbone.6.weight"]]
    b = [tensors["backbone.0.bias"], tensors["backbone.2.bias"],
         tensors["backbone.4.bias"], tensors["backbone.6.bias"]]

    return MLWeights(
        W, b,
        tensors["momentum_mlp.0.weight"], tensors["momentum_mlp.0.bias"],
        tensors["class_head.weight"], tensors["class_head.bias"],
        tensors["inc_head.weight"], tensors["inc_head.bias"],
        tensors["dec_head.weight"], tensors["dec_head.bias"],
        tensors["tau_inc_param"][1], tensors["tau_dec_param"][1],
        tau_floor, dt,
        tensors["norm.mean"], tensors["norm.std"],
        tensors["norm.radsum_mean"], tensors["norm.radsum_std"],
        tensors["norm.cross_mean"], tensors["norm.cross_std"],
    )
end

# ---------------------------------------------------------------------------
# Feature normalization -- mirrors build_feature_matrix_from_raw for
# feature_set="sic+sst+radsum+doy_sin+sfcwind+tas" (+ 4 cross-neighbor cols),
# matching the shifted-notos default: tos NOT shifted, other forcing IS.
# raw order in: [siconc, tos, rsds, rls, sfcwind, tas, doy_sin, doy_cos]
# ---------------------------------------------------------------------------
function ml_normalize_features(w::MLWeights, siconc, tos, rsds, rls, sfcwind, tas, doy_sin, cross::NTuple{N,Float32}) where N
    m = w.norm_mean; s = w.norm_std
    f_sic = (siconc - m[1]) / max(s[1], 1f-6)
    f_tos = (tos - m[2]) / max(s[2], 1f-6)
    f_radsum = ((rsds + rls) - w.radsum_mean) / max(w.radsum_std, 1f-6)
    f_doy = doy_sin  # not normalized
    f_wind = (sfcwind - m[5]) / max(s[5], 1f-6)
    f_tas = (tas - m[6]) / max(s[6], 1f-6)
    cross_norm = ntuple(i -> (cross[i] - w.cross_mean[i]) / max(w.cross_std[i], 1f-6), N)
    # feature_set="sic+sst+radsum+doy_sin+sfcwind+tas" -> 6 cols in THIS order, then N cross cols
    # (N=4: N/S/W/E only; N=8: N/S/W/E + NW/NE/SW/SE, matching the offline 8-cross training convention)
    return vcat(Float32[f_sic, f_tos, f_radsum, f_doy, f_wind, f_tas], collect(Float32, cross_norm))
end

# ---------------------------------------------------------------------------
# Forward pass. hist3 = [sic(t-3), sic(t-2), sic(t-1)] raw (not normalized).
# Returns the predicted tendency (raw units, same scale as sic).
# ---------------------------------------------------------------------------
"""
Weights for the dedicated low-SIC (0.0-0.1) "change-gate" model: replaces the
increase/decrease head split with a 2-way [unchanged, changed] softmax gate plus a
SINGLE unified `changed_head` (sigmoid output is itself the trivial convex combination
of the physical bounds {0,1}), and a single shared `tau_changed_param` (not separate
tau_inc/tau_dec) -- see run_lowsic_changegate_monotonic.py's `RelaxationMoEClassifier`
(change-gate variant) for the matching PyTorch definition. Verified bit-exact
(~4e-8 max abs error) against the trained checkpoint
(checkpoints/north_relaxation_moe_lowsic_changegate_mono100000_tau0005/best.pt).
"""
struct MLWeightsChangeGate
    W::Vector{Matrix{Float32}}
    b::Vector{Vector{Float32}}
    Wm::Matrix{Float32}
    bm::Vector{Float32}
    Wc::Matrix{Float32}; bc::Vector{Float32}     # class_head: [unchanged, changed]
    Wch::Matrix{Float32}; bch::Vector{Float32}   # changed_head
    tau_changed_param::Float32
    tau_floor::Float32
    dt::Float32
    norm_mean::Vector{Float32}
    norm_std::Vector{Float32}
    radsum_mean::Float32
    radsum_std::Float32
    cross_mean::Vector{Float32}
    cross_std::Vector{Float32}
end

function load_ml_weights_changegate(dir::String; tau_floor::Float32=1.0f0, dt::Float32=5.0f0)
    manifest = readlines(joinpath(dir, "manifest.txt"))
    tensors = Dict{String, Any}()
    for line in manifest
        isempty(strip(line)) && continue
        name, shape, fname = split(line, "\t")
        tensors[name] = _read_bin(dir, fname, shape)
    end

    W = [tensors["backbone.0.weight"], tensors["backbone.2.weight"],
         tensors["backbone.4.weight"], tensors["backbone.6.weight"]]
    b = [tensors["backbone.0.bias"], tensors["backbone.2.bias"],
         tensors["backbone.4.bias"], tensors["backbone.6.bias"]]

    return MLWeightsChangeGate(
        W, b,
        tensors["momentum_mlp.0.weight"], tensors["momentum_mlp.0.bias"],
        tensors["class_head.weight"], tensors["class_head.bias"],
        tensors["changed_head.weight"], tensors["changed_head.bias"],
        tensors["tau_changed_param"][1],
        tau_floor, dt,
        tensors["norm.mean"], tensors["norm.std"],
        tensors["norm.radsum_mean"], tensors["norm.radsum_std"],
        tensors["norm.cross_mean"], tensors["norm.cross_std"],
    )
end

function ml_normalize_features(w::MLWeightsChangeGate, siconc, tos, rsds, rls, sfcwind, tas, doy_sin, cross::NTuple{N,Float32}) where N
    m = w.norm_mean; s = w.norm_std
    f_sic = (siconc - m[1]) / max(s[1], 1f-6)
    f_tos = (tos - m[2]) / max(s[2], 1f-6)
    f_radsum = ((rsds + rls) - w.radsum_mean) / max(w.radsum_std, 1f-6)
    f_doy = doy_sin
    f_wind = (sfcwind - m[5]) / max(s[5], 1f-6)
    f_tas = (tas - m[6]) / max(s[6], 1f-6)
    cross_norm = ntuple(i -> (cross[i] - w.cross_mean[i]) / max(w.cross_std[i], 1f-6), N)
    return vcat(Float32[f_sic, f_tos, f_radsum, f_doy, f_wind, f_tas], collect(Float32, cross_norm))
end

function ml_forward(w::MLWeightsChangeGate, x10::Vector{Float32}, sic_raw::Float32, hist3::NTuple{3,Float32})
    h = x10
    for l in 1:4
        h = gelu(w.W[l]' * h .+ w.b[l])
    end
    full_seq = Float32[hist3[1], hist3[2], hist3[3], sic_raw]
    diffs = Float32[full_seq[2]-full_seq[1], full_seq[3]-full_seq[2], full_seq[4]-full_seq[3]]
    mom = gelu(w.Wm' * diffs .+ w.bm)

    combined = vcat(h, mom)  # 80-d

    class_logits = w.Wc' * combined .+ w.bc
    m = max(class_logits[1], class_logits[2])
    e1 = exp(class_logits[1] - m); e2 = exp(class_logits[2] - m)
    p_changed = e2 / (e1 + e2)

    changed_logit = (w.Wch' * combined .+ w.bch)[1]
    sic_eq_changed = sigmoid(changed_logit)

    tau_changed = w.tau_floor + softplus(w.tau_changed_param)
    decay_changed = exp(-w.dt / tau_changed)

    return p_changed * (sic_eq_changed - sic_raw) * (1f0 - decay_changed)
end

function ml_forward(w::MLWeights, x10::Vector{Float32}, sic_raw::Float32, hist3::NTuple{3,Float32})
    h = x10
    for l in 1:4
        h = gelu(w.W[l]' * h .+ w.b[l])
    end
    full_seq = Float32[hist3[1], hist3[2], hist3[3], sic_raw]
    diffs = Float32[full_seq[2]-full_seq[1], full_seq[3]-full_seq[2], full_seq[4]-full_seq[3]]
    mom = gelu(w.Wm' * diffs .+ w.bm)

    combined = vcat(h, mom)  # 80-d

    class_logits = w.Wc' * combined .+ w.bc
    m = max(class_logits[1], class_logits[2])
    e1 = exp(class_logits[1] - m); e2 = exp(class_logits[2] - m)
    p_inc = e1 / (e1 + e2)
    p_dec = e2 / (e1 + e2)

    inc_logit = (w.Wi' * combined .+ w.bi)[1]
    dec_logit = (w.Wd' * combined .+ w.bd)[1]
    sic_eq_inc = sic_raw + (1f0 - sic_raw) * sigmoid(inc_logit)
    sic_eq_dec = sic_raw * sigmoid(dec_logit)

    tau_inc = w.tau_floor + softplus(w.tau_inc_param)
    tau_dec = w.tau_floor + softplus(w.tau_dec_param)
    decay_inc = exp(-w.dt / tau_inc)
    decay_dec = exp(-w.dt / tau_dec)

    tendency_inc = (sic_eq_inc - sic_raw) * (1f0 - decay_inc)
    tendency_dec = (sic_eq_dec - sic_raw) * (1f0 - decay_dec)

    return p_inc * tendency_inc + p_dec * tendency_dec
end

end # module
