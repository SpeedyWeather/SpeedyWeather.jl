"""
ML-based sea ice model for the Arctic (lat >= lat_limit), falling back to
SpeedyWeather's own `ThermodynamicSeaIce` scheme everywhere else (in particular
the whole Antarctic). Wraps the trained `RelaxationMoEClassifier` (see
`ml_seaice_model.jl` for the native Julia forward pass, verified against the
PyTorch checkpoint to ~1e-7 abs error).

Update cadence: every `update_every` (default 5 days), matching the model's own
training block length -- SIC is held constant between updates, same convention
used throughout the Python project's rollout/deployment discussion. Block
timestamps are treated as representing the END of their accumulation window
(not the midpoint) -- purely a labeling convention for interoperability, see
project memory; the model itself only ever learned a discrete block-index
correspondence.

Forcing is a TRUE running mean over each 5-day cycle (not an instantaneous
sample): `timestep!` is called every model timestep regardless of the update
schedule, and always accumulates running sums for rsds/net-longwave/wind/tas/
tos; only when the schedule fires are the sums divided by the step count,
used for the ML update, and reset. This mirrors the offline tumbling-block
cache construction exactly, rather than sampling a single instant.

TOS is deliberately NOT resampled at the same cycle as the other forcing
variables (matches the trained "shifted, notos" convention: other forcing's
mean comes from the block that just completed, TOS's mean comes from the
PREVIOUS block) -- implemented via a one-cycle lag cache, `tos_prev_mean`,
populated from the accumulator that just finished, one cycle after it was
accumulated.

Three trained networks are blended at the output level for Arctic points, matching
the offline-validated recipe: the default (mid-band, 0.1-0.9) model everywhere,
plus a dedicated high-SIC (0.9-1.0) model, combined via a sigmoid gate centered
at `highsic_gate_center` (default 0.9) with width `highsic_gate_width` (default
0.025), PLUS a dedicated low-SIC (0.0-0.1) "change-gate" model, combined via a
second sigmoid gate centered at `lowsic_gate_center` (default 0.05 -- NOT 0.1;
the low-SIC model's own systematic increase-magnitude under-prediction bias grows
sharply toward the 0.1 boundary, so the boundary was pulled in to 0.05 to let the
mid-band model, which already extrapolates reasonably well across the whole
original 0.1-0.9 band, cover the worse-biased 0.05-0.10 zone instead -- see project
memory) with width `lowsic_gate_width` (default 0.025):
`tendency = w_low*tendency_low + (1-w_low-w_high)*tendency_mid + w_high*tendency_high`,
all three models run on every Arctic point each cycle with their OWN normalization
stats. This avoids both the mid-band model's known extrapolation bias for SIC>0.9
(see project memory) and its lack of any dedicated model for SIC<0.1, previously the
single largest source of the online "streaky" prediction artifact (66% of the
worst 1% anomalous points had SIC<0.1 vs a 38.5% base rate).

After the per-point update, an optional spectral (spherical-harmonic) Laplacian diffusion pass
is applied to the WHOLE global SIC field once per cycle (forward-transform to spectral via
`model.spectral_transform`, multiply each degree-l coefficient by the exact one-step diffusion
solution `exp(-D*l*(l+1)/radius^2*Δt)`, inverse-transform back), with only the ARCTIC points'
values committed back (non-Arctic/`ThermodynamicSeaIce` points are left untouched) -- this
mitigates the online "streaky" per-pixel artifact (see project memory: confirmed model/grid-
architecture-driven, not forcing-driven, via a real-data spin-up test) without retraining any of
the three networks. Unlike a naive grid-space finite-difference Laplacian (which has a CFL-type
numerical-instability threshold that bites hardest near the pole, where physical grid spacing on
a lat/lon-like grid is smallest), the spectral exponential-decay form is UNCONDITIONALLY STABLE
for any diffusion coefficient/timestep -- there is no analogous instability cliff to avoid here.
`diffusion_coefficient` (default 26.97 m^2/s) was fit offline (data-optimized against a multi-step
free-running rollout objective, not hand-tuned, not baked into the network) on TRAINING years,
then confirmed on the held-out 1989-1990 val rollout to improve CORR, MAE, AND local-anomaly
(streak) severity simultaneously -- not a smoothness-vs-accuracy tradeoff. That offline fit used a
grid-space finite-difference Laplacian on the 2-degree ERA5 grid (real physical spacing, dx
shrinking ~cos(lat) toward the pole, dy constant) as D_dt=11.65 km^2 per 5-day block; converting
to the SI diffusion coefficient used here: D = 11.65e6 m^2 / (5*86400 s) = 26.97 m^2/s.

SST at Arctic (ML-branch) points is floored every timestep at `tos_floor`
(default 271.46 K) -- NOT SpeedyWeather's own `ThermodynamicSeaIce` constant
(271.35 K), but the empirical floor found directly in the ERA5 `tos` training
data (the field's own 0.1st/1st percentiles coincide at ~271.459 K across every
sampled year, evidently OSTIA's own internal freezing floor). This mirrors the
physical reasoning behind SpeedyWeather's own ice-covered-ocean floor (excess
cooling under ice grows ice instead of cooling the mixed layer further) which the
ML branch had NO equivalent of -- SST there was free to run away arbitrarily far
below freezing with nothing to arrest it, fed into `ml_forward` far outside the
training distribution, and was the proximate cause of a real SST-blowup-to-NaN
failure observed in a T63/8-layer test run once a majority of Arctic points sat
in the (previously unblended) SIC>0.9 extrapolation regime during a real winter.
"""

include("ml_seaice_model.jl")
using .MLSeaIceModel: MLWeights, load_ml_weights, ml_forward, ml_normalize_features
using .MLSeaIceModel: MLWeightsChangeGate, load_ml_weights_changegate

export MLSeaIce

@kwdef mutable struct MLSeaIce <: AbstractSeaIce
    "[OPTION] latitude (degrees north) above which the ML model is used; below, falls back to the thermodynamic scheme"
    lat_limit::Float64 = 60.0

    "[OPTION] update interval, matching the training block length"
    update_every::Second = Day(5)

    "[OPTION] path to the exported weight/manifest directory (mid-band, SIC 0.1-0.9, model).
    Defaults to the monotonicity-penalty-trained variant (lambda_mono=1000) -- a direct
    sensitivity audit found the un-penalized model responded to radsum in the physically wrong
    direction (more radiative heating -> less predicted melt) for 88.2% of real val samples;
    this variant fixes that (99%+ correct sign) at a modest regression-accuracy cost (CORR
    0.750->0.718 offline). The original (un-penalized) weights remain at
    ml_seaice_weights/ for comparison/rollback."
    weights_dir::String = joinpath(@__DIR__, "ml_seaice_weights_mono1000")

    "[OPTION] path to the exported weight/manifest directory for the dedicated high-SIC (0.9-1.0)
    model. Defaults to the monotonicity-penalty-trained variant (lambda_mono=10 -- chosen lower
    than the mid-band's 1000 since it was the smallest value that still clearly beat naive
    mid-band-extrapolation onto this band, preserving the original reason for having a dedicated
    high-SIC model at all). The original (un-penalized) weights remain at
    ml_seaice_weights_highsic/ for comparison/rollback."
    weights_dir_highsic::String = joinpath(@__DIR__, "ml_seaice_weights_highsic_mono10")

    "[OPTION] sigmoid gate center (in SIC units) for blending the mid-band and high-SIC models"
    highsic_gate_center::Float32 = 0.9f0

    "[OPTION] sigmoid gate width (in SIC units) for blending the mid-band and high-SIC models"
    highsic_gate_width::Float32 = 0.025f0

    "[OPTION] path to the exported weight/manifest directory for the dedicated low-SIC (0.0-0.1)
    change-gate model. Adopted config: lambda_mono=100000, unchanged_tau=0.005 -- offline CORR=0.735
    vs naive mid-band-extrapolation's 0.630 on this band. The original increase/decrease-head
    architecture and two post-hoc training-data-resampling attempts were both tried and abandoned
    (ambiguous / net-negative results, see project memory); this is the plain uniform-cache-trained
    change-gate checkpoint."
    weights_dir_lowsic::String = joinpath(@__DIR__, "ml_seaice_weights_lowsic_changegate_mono100000_tau0005")

    "[OPTION] sigmoid gate center (in SIC units) for blending the low-SIC and mid-band models.
    Deliberately 0.05, not the low-SIC model's own training-band edge of 0.1 -- an offline soft-blend
    test at boundary=0.1 caused a net regression (mean CORR 0.9439->0.9407) concentrated in freeze-up
    months, traced to the low-SIC model's increase-magnitude bias worsening sharply toward 0.1; at
    boundary=0.05 the regression disappears entirely (mean CORR back to 0.9439) since the low-SIC
    model is only trusted in its more reliable sub-range and the mid-band model covers 0.05-0.10."
    lowsic_gate_center::Float32 = 0.05f0

    "[OPTION] sigmoid gate width (in SIC units) for blending the low-SIC and mid-band models"
    lowsic_gate_width::Float32 = 0.025f0

    "[OPTION] path to a flat Float32 binary with real ERA5-derived initial SIC + 3-block history
    (layout: hist3, hist2, hist1, current, each `npoints` Float32, in the SpeedyWeather grid's own
    point order for this exact configuration); set to \"\" to skip and use the flat/zero seeding"
    initial_sic_path::String = joinpath(@__DIR__, "ml_seaice_weights", "initial_sic.bin")

    "[OPTION] SST floor [K] applied every timestep at Arctic (ML-branch) points -- kept at the
    physical seawater freezing point (273.15-1.8=271.35K, SpeedyWeather's own ThermodynamicSeaIce
    constant) rather than the empirical ERA5 tos training-data floor (~271.459K); the latter is
    slightly higher than the true physical floor, and the model is expected to extrapolate a small
    (~0.11K) amount below its training distribution rather than being artificially prevented from
    ever reaching the physical freezing point"
    tos_floor::Float64 = 273.15 - 1.8

    "[OPTION] fallback thermodynamic scheme parameters (same defaults as ThermodynamicSeaIce), used
    for non-Arctic ocean points only"
    temp_freeze::Float64 = 273.15 - 1.8
    melt_rate::Float64 = 1.0e-6
    freeze_rate::Float64 = 0.1

    "[OPTION] apply an unconditionally-stable spectral (spherical-harmonic) Laplacian diffusion
    pass to the SIC field once per update cycle, at Arctic points only -- mitigates the online
    per-pixel streaky artifact (confirmed model/grid-architecture-driven, not forcing-driven, see
    project memory) without retraining. Set false to recover the original undiffused behavior."
    apply_spectral_diffusion::Bool = true

    "[OPTION] diffusion coefficient [m^2/s] for the spectral Laplacian pass. Default 26.97 m^2/s
    was data-fitted offline (multi-step free-running rollout objective, not hand-tuned) and
    confirmed on held-out validation data to improve accuracy (CORR/MAE) AND reduce streak
    severity simultaneously -- see weights_dir_lowsic-adjacent project memory / this struct's own
    docstring for the exact fitting methodology and the km^2-to-m^2/s unit conversion."
    diffusion_coefficient::Float64 = 26.97

    "[OPTION] minimum locally-diffused ocean weight required to trust the normalized-convolution
    diffusion estimate at a point; below this, the point keeps its pre-diffusion (ML-predicted)
    value for this cycle instead. Points near a lot of land (thin local ocean support) have a
    small denominator here, which can amplify noise if trusted -- see project memory for the
    empirical coastal-vs-interior tuning of this threshold."
    diffusion_min_weight::Float64 = 0.3

    # ======================================================================
    # Southern Hemisphere (Antarctic) options -- an independently-trained set of three band
    # experts (mid/low/high, same architectures as the Arctic's own), NOT the Arctic checkpoints
    # reused: the two hemispheres' sea-ice regimes differ enough (seasonal vs. perennial ice,
    # open-ocean vs. land-locked geometry -- see project memory) that every hyperparameter below
    # was independently re-derived for Antarctica, not copied from the Arctic defaults above,
    # even where the final adopted value happens to coincide (e.g. the high-SIC gate center).
    # ======================================================================

    "[OPTION] enable the Southern Hemisphere ML branch. When false, the Antarctic (and all other
    non-Arctic) ocean falls back entirely to ThermodynamicSeaIce, matching this component's
    original (Arctic-only) behavior."
    enable_south::Bool = true

    "[OPTION] latitude (degrees, positive number) below which (i.e. south of -lat_limit_south) the
    Southern ML branch is used -- re-derived independently for Antarctica (55), NOT a copy of the
    Arctic's 60: at 60S the candidate band is only ~42% ocean (dominated by the Antarctic
    continent itself) and would clip off real winter sea ice, which extends to ~55-60S -- see
    project memory."
    lat_limit_south::Float64 = 55.0

    "[OPTION] path to the exported weight/manifest directory for the Antarctic mid-band (SIC
    0.1-0.9) model. Same architecture as the Arctic mid-band (RelaxationMoEClassifier) and the
    same monotonicity-penalty fix (radsum sign-inversion independently re-diagnosed for
    Antarctica, not assumed to transfer -- found at similar severity, 66.7% backwards pre-fix)."
    weights_dir_south_mid::String = joinpath(@__DIR__, "ml_seaice_weights_south_mid_mono1000")

    "[OPTION] path to the exported weight/manifest directory for the Antarctic dedicated low-SIC
    (0.0-0.1) change-gate model."
    weights_dir_south_lowsic::String = joinpath(@__DIR__, "ml_seaice_weights_south_lowsic_changegate_mono10000")

    "[OPTION] sigmoid gate center (SIC units) for blending the Antarctic low-SIC and mid-band
    models. Independently re-derived via an offline rollout boundary sweep (0.03/0.05/0.08/0.10/
    0.15) -- CORR improved monotonically as this boundary moved down, landing at 0.03, tighter
    than the Arctic's own adopted 0.05 (the Antarctic low-SIC expert's own advantage over the
    mid-band's extrapolation is smaller than the Arctic's equivalent gap, so less territory is
    given to it here -- see project memory)."
    lowsic_gate_center_south::Float32 = 0.03f0

    "[OPTION] sigmoid gate width (SIC units), Southern low-SIC boundary."
    lowsic_gate_width_south::Float32 = 0.025f0

    "[OPTION] apply the same unconditionally-stable spectral Laplacian diffusion pass to Antarctic
    points, with the Southern Hemisphere's own independently-fitted coefficient (NOT the Arctic's
    D). Applied as a SEPARATE spectral transform from the Arctic pass (see timestep!) since the
    two hemispheres use different D and spectral diffusion is a global (whole-field, per-degree-l)
    operation that cannot mix two different decay rates in one transform."
    apply_spectral_diffusion_south::Bool = true

    "[OPTION] diffusion coefficient [m^2/s] for the Southern Hemisphere spectral Laplacian pass.
    Default 1388.9 m^2/s = 600 km^2/block / 432000s, fit offline via the same multi-step
    free-running-rollout sweep methodology as the Arctic's D=26.97, on TRAINING years with the
    chosen D only then validated on genuinely held-out data. Refit 2026-08-06 for the CORRECTED
    band-expert checkpoints and the 2-way (not 3-way) blend below -- supersedes an earlier
    D*dt=900 (D=2083.3) fit against checkpoints later found to have a train-gradient
    normalization-leakage bug (see project memory); the corrected sweep's peak moved down to
    D*dt~400-600, held-out validation selected D*dt=600 (CORR=0.96942, the best of the tested
    values) with a large safety margin below the ~1400-1600 instability onset. Still substantially
    larger than the Arctic's own D=26.97 -- an empirical finding, not yet explained by a specific
    geometric argument."
    diffusion_coefficient_south::Float64 = 1388.9

    "[OPTION] minimum locally-diffused ocean weight for the Southern Hemisphere pass, same role as
    diffusion_min_weight above."
    diffusion_min_weight_south::Float64 = 0.3

    "[OPTION] path to a flat Float32 binary with real ERA5-derived initial Antarctic SIC + 3-block
    history, same layout as initial_sic_path (hist3, hist2, hist1, current, each `npoints`
    Float32, this exact grid configuration's own point order); set to \"\" to skip and use the
    flat/zero seeding."
    initial_sic_path_south::String = joinpath(@__DIR__, "ml_seaice_weights", "initial_sic_south.bin")

    # -- allocated in initialize!(vars, model) once the grid is known --
    weights::Union{Nothing, MLWeights} = nothing
    weights_highsic::Union{Nothing, MLWeights} = nothing
    weights_lowsic::Union{Nothing, MLWeightsChangeGate} = nothing
    weights_south::Union{Nothing, MLWeights} = nothing
    weights_lowsic_south::Union{Nothing, MLWeightsChangeGate} = nothing
    is_arctic::Vector{Bool} = zeros(Bool, 0)          # per-gridpoint mask, lat >= lat_limit & ocean
    is_antarctic::Vector{Bool} = zeros(Bool, 0)       # per-gridpoint mask, lat <= -lat_limit_south & ocean
    east_idx::Vector{Int} = zeros(Int, 0)             # cross-neighbor indices (0 = missing/land -> zero diff)
    west_idx::Vector{Int} = zeros(Int, 0)
    north_idx::Vector{Int} = zeros(Int, 0)
    south_idx::Vector{Int} = zeros(Int, 0)
    # diagonal neighbors (composed from the 4 cardinal arrays), for the 8-cross (3x3) models
    northwest_idx::Vector{Int} = zeros(Int, 0)
    northeast_idx::Vector{Int} = zeros(Int, 0)
    southwest_idx::Vector{Int} = zeros(Int, 0)
    southeast_idx::Vector{Int} = zeros(Int, 0)

    # running-mean accumulators for the CURRENT (not-yet-finalized) cycle
    accum_rsds::Vector{Float64} = zeros(0)
    accum_lw::Vector{Float64} = zeros(0)              # net longwave (down - up)
    accum_wind::Vector{Float64} = zeros(0)
    accum_tas::Vector{Float64} = zeros(0)
    accum_tos::Vector{Float64} = zeros(0)
    accum_count::Int = 0

    # diagnostics only -- most recent finalized cycle's per-point radiative forcing, for offline
    # comparison against the training-data distribution (not used by the model itself)
    diag_rsds_mean::Vector{Float64} = zeros(0)
    diag_lw_mean::Vector{Float64} = zeros(0)

    tos_prev_mean::Vector{Float64} = zeros(0)         # finalized TOS mean from the PREVIOUS cycle (one-cycle lag)
    sic_hist1::Vector{Float64} = zeros(0)             # SIC 3, 2, 1 updates ago
    sic_hist2::Vector{Float64} = zeros(0)
    sic_hist3::Vector{Float64} = zeros(0)
    initialized::Vector{Bool} = zeros(Bool, 0)        # per-point: has the ML model been seeded with real history yet?

    last_update::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
end

MLSeaIce(::SpectralGrid; kwargs...) = MLSeaIce(; kwargs...)

function initialize!(sea_ice_model::MLSeaIce, model::AbstractModel)
    sea_ice_model.weights = load_ml_weights(sea_ice_model.weights_dir)
    sea_ice_model.weights_highsic = load_ml_weights(sea_ice_model.weights_dir_highsic)
    sea_ice_model.weights_lowsic = load_ml_weights_changegate(sea_ice_model.weights_dir_lowsic)
    if sea_ice_model.enable_south
        sea_ice_model.weights_south = load_ml_weights(sea_ice_model.weights_dir_south_mid)
        sea_ice_model.weights_lowsic_south = load_ml_weights_changegate(sea_ice_model.weights_dir_south_lowsic)
    end
    return nothing
end

function initialize!(vars::Variables, sea_ice_model::MLSeaIce, model::PrimitiveEquation)
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    grid = ℵ.grid
    npoints = length(ℵ)
    londs, latds = RingGrids.get_londlatds(grid)

    (; land_fraction) = model.land_sea_mask
    is_arctic = [(latds[ij] >= sea_ice_model.lat_limit) && (land_fraction[ij] < 1) for ij in 1:npoints]
    is_antarctic = sea_ice_model.enable_south ?
        [(latds[ij] <= -sea_ice_model.lat_limit_south) && (land_fraction[ij] < 1) for ij in 1:npoints] :
        zeros(Bool, npoints)

    # cross-neighbor indices: east/west trivial within a ring; north/south via
    # nearest-longitude match in the adjacent ring (exact for full grids, an
    # approximation for reduced grids where ring lengths differ)
    rings = grid.rings
    ring_of = zeros(Int, npoints)
    for (j, r) in enumerate(rings), ij in r
        ring_of[ij] = j
    end

    east_idx = zeros(Int, npoints)
    west_idx = zeros(Int, npoints)
    north_idx = zeros(Int, npoints)
    south_idx = zeros(Int, npoints)

    for (j, r) in enumerate(rings)
        n = length(r)
        for (k, ij) in enumerate(r)
            east_idx[ij] = r[mod1(k + 1, n)]
            west_idx[ij] = r[mod1(k - 1, n)]
        end
    end

    function nearest_in_ring(ij, target_ring)
        (target_ring < 1 || target_ring > length(rings)) && return 0
        r = rings[target_ring]
        isempty(r) && return 0
        lon0 = londs[ij]
        best = r[1]; bestd = 361.0
        for jk in r
            d = abs(mod(londs[jk] - lon0 + 180, 360) - 180)
            if d < bestd
                bestd = d; best = jk
            end
        end
        return best
    end

    for ij in 1:npoints
        j = ring_of[ij]
        north_idx[ij] = nearest_in_ring(ij, j - 1)   # ring j-1 assumed poleward-adjacent (grid-convention dependent)
        south_idx[ij] = nearest_in_ring(ij, j + 1)
    end

    # diagonal neighbors: compose the cardinal arrays (go north/south first, then west/east).
    # west_idx/east_idx are always valid (same-ring wraparound, no missing sentinel needed);
    # only north_idx/south_idx can be 0 (missing, e.g. beyond the pole), so that's the only
    # guard needed here.
    northwest_idx = zeros(Int, npoints)
    northeast_idx = zeros(Int, npoints)
    southwest_idx = zeros(Int, npoints)
    southeast_idx = zeros(Int, npoints)
    for ij in 1:npoints
        n = north_idx[ij]
        s = south_idx[ij]
        northwest_idx[ij] = n == 0 ? 0 : west_idx[n]
        northeast_idx[ij] = n == 0 ? 0 : east_idx[n]
        southwest_idx[ij] = s == 0 ? 0 : west_idx[s]
        southeast_idx[ij] = s == 0 ? 0 : east_idx[s]
    end

    sea_ice_model.is_arctic = is_arctic
    sea_ice_model.is_antarctic = is_antarctic
    sea_ice_model.east_idx = east_idx
    sea_ice_model.west_idx = west_idx
    sea_ice_model.north_idx = north_idx
    sea_ice_model.south_idx = south_idx
    sea_ice_model.northwest_idx = northwest_idx
    sea_ice_model.northeast_idx = northeast_idx
    sea_ice_model.southwest_idx = southwest_idx
    sea_ice_model.southeast_idx = southeast_idx

    sea_ice_model.accum_rsds = zeros(npoints)
    sea_ice_model.accum_lw = zeros(npoints)
    sea_ice_model.accum_wind = zeros(npoints)
    sea_ice_model.accum_tas = zeros(npoints)
    sea_ice_model.accum_tos = zeros(npoints)
    sea_ice_model.accum_count = 0

    sea_ice_model.diag_rsds_mean = zeros(npoints)
    sea_ice_model.diag_lw_mean = zeros(npoints)
    sea_ice_model.tos_prev_mean = fill(NaN, npoints)
    sea_ice_model.sic_hist1 = fill(NaN, npoints)
    sea_ice_model.sic_hist2 = fill(NaN, npoints)
    sea_ice_model.sic_hist3 = fill(NaN, npoints)
    sea_ice_model.initialized = zeros(Bool, npoints)
    sea_ice_model.last_update[] = vars.prognostic.clock.time

    # seed Arctic points with real ERA5-derived SIC + history instead of the flat/zero default;
    # the file stores hist3, hist2, hist1, current back-to-back, each `npoints` Float32, in this
    # exact grid configuration's own point order (see project memory / export script)
    if !isempty(sea_ice_model.initial_sic_path) && isfile(sea_ice_model.initial_sic_path)
        raw = reinterpret(Float32, read(sea_ice_model.initial_sic_path))
        length(raw) == 4 * npoints || error("initial_sic.bin has $(length(raw)) values, expected $(4 * npoints)")
        hist3 = raw[1:npoints]
        hist2 = raw[npoints+1:2npoints]
        hist1 = raw[2npoints+1:3npoints]
        current = raw[3npoints+1:4npoints]
        sst0 = vars.prognostic.ocean.sea_surface_temperature
        for ij in 1:npoints
            if is_arctic[ij]
                ℵ[ij] = clamp(current[ij], 0f0, 1f0)
                sea_ice_model.sic_hist3[ij] = clamp(hist3[ij], 0f0, 1f0)
                sea_ice_model.sic_hist2[ij] = clamp(hist2[ij], 0f0, 1f0)
                sea_ice_model.sic_hist1[ij] = clamp(hist1[ij], 0f0, 1f0)
                # tos_prev_mean would otherwise stay NaN forever: marking `initialized[ij]=true`
                # below skips the flat-start fallback in timestep! that used to populate it from
                # the first cycle's own accumulated SST -- seed it from the model's actual initial
                # SST instead, one real cycle's worth of lag error at most.
                sea_ice_model.tos_prev_mean[ij] = isfinite(sst0[ij]) ? Float64(sst0[ij]) : sea_ice_model.temp_freeze
                sea_ice_model.initialized[ij] = true
            end
        end
    end

    # seed Antarctic points the same way, from an independent real ERA5-derived file (same
    # layout) -- not the Arctic file above, a different hemisphere's real initial condition
    if sea_ice_model.enable_south && !isempty(sea_ice_model.initial_sic_path_south) && isfile(sea_ice_model.initial_sic_path_south)
        raw_s = reinterpret(Float32, read(sea_ice_model.initial_sic_path_south))
        length(raw_s) == 4 * npoints || error("initial_sic_south.bin has $(length(raw_s)) values, expected $(4 * npoints)")
        hist3_s = raw_s[1:npoints]
        hist2_s = raw_s[npoints+1:2npoints]
        hist1_s = raw_s[2npoints+1:3npoints]
        current_s = raw_s[3npoints+1:4npoints]
        sst0 = vars.prognostic.ocean.sea_surface_temperature
        for ij in 1:npoints
            if is_antarctic[ij]
                ℵ[ij] = clamp(current_s[ij], 0f0, 1f0)
                sea_ice_model.sic_hist3[ij] = clamp(hist3_s[ij], 0f0, 1f0)
                sea_ice_model.sic_hist2[ij] = clamp(hist2_s[ij], 0f0, 1f0)
                sea_ice_model.sic_hist1[ij] = clamp(hist1_s[ij], 0f0, 1f0)
                sea_ice_model.tos_prev_mean[ij] = isfinite(sst0[ij]) ? Float64(sst0[ij]) : sea_ice_model.temp_freeze
                sea_ice_model.initialized[ij] = true
            end
        end
    end

    return nothing
end

function timestep!(vars::Variables, sea_ice_model::MLSeaIce, model::PrimitiveEquation)
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    haskey(vars.prognostic.ocean, :sea_surface_temperature) || return nothing
    sst = vars.prognostic.ocean.sea_surface_temperature

    rsds = vars.parameterizations.surface_shortwave_down
    Rld = vars.parameterizations.surface_longwave_down
    Rlu = vars.parameterizations.ocean.surface_longwave_up
    wind10 = vars.parameterizations.surface_wind_speed
    tas = vars.parameterizations.surface_air_temperature

    npoints = length(ℵ)
    (; land_fraction) = model.land_sea_mask
    (; accum_rsds, accum_lw, accum_wind, accum_tas, accum_tos, is_arctic, is_antarctic) = sea_ice_model
    tos_floor = Float32(sea_ice_model.tos_floor)

    # --- ALWAYS accumulate this timestep's instantaneous values (true running mean) ---
    @inbounds for ij in 1:npoints
        if land_fraction[ij] < 1 && isfinite(sst[ij])
            # Arctic AND Antarctic (ML-branch) points: floor SST every step, mirroring the
            # physical reasoning behind SpeedyWeather's own ice-covered-ocean floor (excess
            # cooling under ice grows ice, not further mixed-layer cooling) -- the ML branch has
            # no other mechanism to arrest runaway cooling, unlike the non-ML fallback below which
            # already has one. This is exactly the bug this project found and fixed for the
            # Arctic branch (an SST-blowup-to-NaN failure at T63/8-layer, project memory) --
            # applied to Antarctic points from the start here, not discovered the hard way twice.
            if (is_arctic[ij] || is_antarctic[ij]) && sst[ij] < tos_floor
                sst[ij] = tos_floor
            end
            accum_rsds[ij] += Float64(rsds[ij])
            accum_lw[ij] += Float64(Rld[ij] - Rlu[ij])
            accum_wind[ij] += Float64(wind10[ij])
            accum_tas[ij] += Float64(tas[ij])
            accum_tos[ij] += Float64(sst[ij])
        end
    end
    sea_ice_model.accum_count += 1

    now = vars.prognostic.clock.time
    (now - sea_ice_model.last_update[]) < sea_ice_model.update_every && return nothing
    sea_ice_model.last_update[] = now

    n = Float64(sea_ice_model.accum_count)
    w = sea_ice_model.weights
    w_high = sea_ice_model.weights_highsic
    w_low = sea_ice_model.weights_lowsic
    gate_center = sea_ice_model.highsic_gate_center
    gate_width = sea_ice_model.highsic_gate_width
    low_gate_center = sea_ice_model.lowsic_gate_center
    low_gate_width = sea_ice_model.lowsic_gate_width
    w_south = sea_ice_model.weights_south
    w_low_south = sea_ice_model.weights_lowsic_south
    low_gate_center_south = sea_ice_model.lowsic_gate_center_south
    low_gate_width_south = sea_ice_model.lowsic_gate_width_south
    (; east_idx, west_idx, north_idx, south_idx) = sea_ice_model
    (; northwest_idx, northeast_idx, southwest_idx, southeast_idx) = sea_ice_model
    (; tos_prev_mean, sic_hist1, sic_hist2, sic_hist3, initialized) = sea_ice_model
    doy = Dates.dayofyear(now)
    doy_sin = Float32(sin(2π * doy / 366))
    Δt_days = Float32(Dates.value(Second(sea_ice_model.update_every)) / 86400)

    for ij in 1:npoints
        if is_arctic[ij] && isfinite(sst[ij])
            sic_now = Float32(ℵ[ij])
            rsds_mean = Float32(accum_rsds[ij] / n)
            lw_mean = Float32(accum_lw[ij] / n)
            wind_mean = Float32(accum_wind[ij] / n)
            tas_mean = Float32(accum_tas[ij] / n)
            tos_mean_this_cycle = Float32(accum_tos[ij] / n)
            sea_ice_model.diag_rsds_mean[ij] = rsds_mean
            sea_ice_model.diag_lw_mean[ij] = lw_mean

            if !initialized[ij]
                # seed history with the current (constant) state -- first update after
                # spin-up has no real history; treated as a flat/no-momentum start.
                sic_hist1[ij] = sic_now; sic_hist2[ij] = sic_now; sic_hist3[ij] = sic_now
                tos_prev_mean[ij] = tos_mean_this_cycle
                initialized[ij] = true
            end

            @inline function neighbor_diff(nidx::Int)
                (nidx == 0 || land_fraction[nidx] >= 1) && return 0f0
                return Float32(ℵ[nidx] - sic_now)
            end
            # full 8-direction (3x3) neighborhood; sliced per-model below to match however many
            # cross dims that model's own loaded weights expect (4: cardinal only: N/S/W/E;
            # 8: cardinal + diagonal, matching the offline 8-cross training convention -- order
            # here MUST match build_cache_8cross.py's cross_diff: N,S,W,E,NW,NE,SW,SE)
            cross8 = (
                neighbor_diff(north_idx[ij]),
                neighbor_diff(south_idx[ij]),
                neighbor_diff(west_idx[ij]),
                neighbor_diff(east_idx[ij]),
                neighbor_diff(northwest_idx[ij]),
                neighbor_diff(northeast_idx[ij]),
                neighbor_diff(southwest_idx[ij]),
                neighbor_diff(southeast_idx[ij]),
            )

            hist3 = (Float32(sic_hist3[ij]), Float32(sic_hist2[ij]), Float32(sic_hist1[ij]))

            # mid-band (default) model
            x10 = ml_normalize_features(
                w, sic_now, tos_prev_mean[ij],
                rsds_mean, lw_mean, wind_mean, tas_mean, doy_sin,
                cross8[1:length(w.cross_mean)],
            )
            tendency_mid = ml_forward(w, x10, sic_now, hist3)

            # dedicated high-SIC model, own normalization stats, same raw features
            x10_high = ml_normalize_features(
                w_high, sic_now, tos_prev_mean[ij],
                rsds_mean, lw_mean, wind_mean, tas_mean, doy_sin,
                cross8[1:length(w_high.cross_mean)],
            )
            tendency_high = ml_forward(w_high, x10_high, sic_now, hist3)

            # dedicated low-SIC change-gate model, own normalization stats, same raw features
            x10_low = ml_normalize_features(
                w_low, sic_now, tos_prev_mean[ij],
                rsds_mean, lw_mean, wind_mean, tas_mean, doy_sin,
                cross8[1:length(w_low.cross_mean)],
            )
            tendency_low = ml_forward(w_low, x10_low, sic_now, hist3)

            # sigmoid-gated 3-way soft blend, offline-validated recipe (see project memory):
            # avoids the mid-band model's extrapolation bias for SIC>0.9 AND its lack of a
            # dedicated model for SIC<0.1 -- low-SIC gate centered at 0.05 (not the low-SIC
            # model's own training-band edge of 0.1), see weights_dir_lowsic docstring for why
            gate_high = 1f0 / (1f0 + exp(-(sic_now - gate_center) / gate_width))
            gate_low = 1f0 / (1f0 + exp(-(low_gate_center - sic_now) / low_gate_width))
            gate_mid = 1f0 - gate_low - gate_high
            tendency = gate_low * tendency_low + gate_mid * tendency_mid + gate_high * tendency_high

            sic_next = clamp(sic_now + tendency, 0f0, 1f0)

            sic_hist3[ij] = sic_hist2[ij]; sic_hist2[ij] = sic_hist1[ij]; sic_hist1[ij] = sic_now
            tos_prev_mean[ij] = tos_mean_this_cycle   # cache THIS cycle's TOS mean for next cycle's update
            ℵ[ij] = sic_next
        elseif is_antarctic[ij] && isfinite(sst[ij])
            # Southern Hemisphere branch -- identical structure to the Arctic branch above (same
            # accumulate-then-5-day-update convention, same cross-neighbor features, same 3-way
            # soft blend), but its OWN independently-trained models and OWN independently-derived
            # gate boundaries (see struct docstring: NOT the Arctic's weights or 0.05/0.9 gates).
            sic_now = Float32(ℵ[ij])
            rsds_mean = Float32(accum_rsds[ij] / n)
            lw_mean = Float32(accum_lw[ij] / n)
            wind_mean = Float32(accum_wind[ij] / n)
            tas_mean = Float32(accum_tas[ij] / n)
            tos_mean_this_cycle = Float32(accum_tos[ij] / n)
            sea_ice_model.diag_rsds_mean[ij] = rsds_mean
            sea_ice_model.diag_lw_mean[ij] = lw_mean

            if !initialized[ij]
                sic_hist1[ij] = sic_now; sic_hist2[ij] = sic_now; sic_hist3[ij] = sic_now
                tos_prev_mean[ij] = tos_mean_this_cycle
                initialized[ij] = true
            end

            @inline function neighbor_diff_s(nidx::Int)
                (nidx == 0 || land_fraction[nidx] >= 1) && return 0f0
                return Float32(ℵ[nidx] - sic_now)
            end
            cross8_s = (
                neighbor_diff_s(north_idx[ij]),
                neighbor_diff_s(south_idx[ij]),
                neighbor_diff_s(west_idx[ij]),
                neighbor_diff_s(east_idx[ij]),
                neighbor_diff_s(northwest_idx[ij]),
                neighbor_diff_s(northeast_idx[ij]),
                neighbor_diff_s(southwest_idx[ij]),
                neighbor_diff_s(southeast_idx[ij]),
            )

            hist3 = (Float32(sic_hist3[ij]), Float32(sic_hist2[ij]), Float32(sic_hist1[ij]))

            # 2-way blend (low-SIC change-gate + mid-band): the dedicated Antarctic high-SIC expert
            # was DROPPED 2026-08-06 -- once its own train-gradient normalization-leakage bug was
            # fixed, its honest held-out CORR (0.351) was clearly beaten by simply extrapolating
            # the (also corrected) mid-band model onto the 0.9-1.0 range (CORR=0.473) -- see
            # project memory. The mid-band model now covers the whole range above the low-SIC gate.
            x10_mid_s = ml_normalize_features(
                w_south, sic_now, tos_prev_mean[ij],
                rsds_mean, lw_mean, wind_mean, tas_mean, doy_sin,
                cross8_s[1:length(w_south.cross_mean)],
            )
            tendency_mid_s = ml_forward(w_south, x10_mid_s, sic_now, hist3)

            x10_low_s = ml_normalize_features(
                w_low_south, sic_now, tos_prev_mean[ij],
                rsds_mean, lw_mean, wind_mean, tas_mean, doy_sin,
                cross8_s[1:length(w_low_south.cross_mean)],
            )
            tendency_low_s = ml_forward(w_low_south, x10_low_s, sic_now, hist3)

            gate_low_s = 1f0 / (1f0 + exp(-(low_gate_center_south - sic_now) / low_gate_width_south))
            gate_mid_s = 1f0 - gate_low_s
            tendency_s = gate_low_s * tendency_low_s + gate_mid_s * tendency_mid_s

            sic_next_s = clamp(sic_now + tendency_s, 0f0, 1f0)

            sic_hist3[ij] = sic_hist2[ij]; sic_hist2[ij] = sic_hist1[ij]; sic_hist1[ij] = sic_now
            tos_prev_mean[ij] = tos_mean_this_cycle
            ℵ[ij] = sic_next_s
        elseif land_fraction[ij] < 1 && isfinite(sst[ij])
            # fallback: SpeedyWeather's own thermodynamic scheme (Southern Hemisphere / non-Arctic ocean).
            # Uses this cycle's mean SST (not instantaneous) for consistency with the ML branch, though
            # the original ThermodynamicSeaIce runs every step on the instantaneous value -- a harmless
            # difference since this scheme has no memory/history term to be sensitive to the distinction.
            sst_mean = Float32(accum_tos[ij] / n)
            dT = sst_mean - sea_ice_model.temp_freeze
            F = -sea_ice_model.melt_rate * max(dT, 0) - sea_ice_model.freeze_rate / (Δt_days * 86400) * min(dT, 0)
            sst[ij] = max(sst[ij], sea_ice_model.temp_freeze)
            ℵ[ij] = clamp(ℵ[ij] + Δt_days * 86400 * F, 0, 1)
        end
    end

    # optional spectral Laplacian diffusion pass (Arctic points only), unconditionally stable --
    # see struct docstring for the fitting methodology and why this is preferred over a grid-space
    # finite-difference Laplacian (which has a CFL-type instability threshold near the pole).
    #
    # Normalized-convolution (weight-masked) diffusion, NOT a raw transform of ℵ directly: land
    # points hold a fixed literal 0 SIC (never written by this component, left at SpeedyWeather's
    # own default) -- that is "no data", not "zero ice", so transforming ℵ as-is would treat every
    # coastline as a real, hard SIC discontinuity and pollute the diffusion there with spurious
    # Gibbs ringing, on top of (and far more numerous than) genuine ice-edge fronts. Instead diffuse
    # ℵ*weight and weight separately (weight = ocean fraction) and divide -- the same technique the
    # offline cache-building pipeline already uses for missing/invalid data
    # (spatial_gaussian_smooth_coarse's filled/weight normalization), just in spectral space here.
    if sea_ice_model.apply_spectral_diffusion
        S = model.spectral_transform
        radius = model.planet.radius
        D = Float32(sea_ice_model.diffusion_coefficient)
        Δt_seconds = Float32(Dates.value(Second(sea_ice_model.update_every)))

        weight_field = copy(ℵ)
        filled_field = copy(ℵ)
        for ij in 1:npoints
            wij = 1f0 - Float32(land_fraction[ij])
            weight_field[ij] = wij
            filled_field[ij] = ℵ[ij] * wij
        end

        coeffs_filled = transform(filled_field, S)
        coeffs_weight = transform(weight_field, S)
        (; l_indices) = coeffs_filled.spectrum
        for lm in eachindex(coeffs_filled)
            l = l_indices[lm]
            decay = exp(-D * l * (l + 1) / radius^2 * Δt_seconds)
            coeffs_filled[lm] *= decay
            coeffs_weight[lm] *= decay
        end
        filled_diffused = transform(coeffs_filled, S)
        weight_diffused = transform(coeffs_weight, S)

        min_weight = Float32(sea_ice_model.diffusion_min_weight)
        for ij in 1:npoints
            if is_arctic[ij]
                wd = weight_diffused[ij]
                if wd > min_weight
                    # clamp: the spectral (band-limited) reconstruction can still overshoot/
                    # undershoot [0,1] near sharp real ice-edge fronts (Gibbs ringing at this
                    # truncation) -- an unclamped value here previously fed a corrupted SIC into
                    # the sea-ice-insulation heat flux term, cascading into an SST blowup
                    ℵ[ij] = clamp(filled_diffused[ij] / wd, 0f0, 1f0)
                end
                # else: insufficient local ocean support (mostly land-surrounded) for a reliable
                # diffused estimate -- leave ℵ[ij] at this cycle's pre-diffusion (ML-predicted) value
            end
        end
    end

    # Southern Hemisphere spectral diffusion pass -- a SEPARATE spectral transform from the
    # Arctic one above, with its own D, run afterward on the (by now Arctic-updated, but
    # Antarctic-still-pre-diffusion) ℵ field. Not combinable into one transform: spectral
    # diffusion applies one decay rate per spherical-harmonic degree l to the WHOLE field at
    # once, so two different D's for two different regions genuinely need two transforms, not
    # one. Arctic and Antarctic point sets are disjoint (>=60N vs. <=-55N), so running this
    # second pass after the first, on the partially-updated ℵ, is safe -- each pass only ever
    # writes back its own hemisphere's points.
    if sea_ice_model.enable_south && sea_ice_model.apply_spectral_diffusion_south
        S = model.spectral_transform
        radius = model.planet.radius
        D_s = Float32(sea_ice_model.diffusion_coefficient_south)
        Δt_seconds = Float32(Dates.value(Second(sea_ice_model.update_every)))

        weight_field_s = copy(ℵ)
        filled_field_s = copy(ℵ)
        for ij in 1:npoints
            wij = 1f0 - Float32(land_fraction[ij])
            weight_field_s[ij] = wij
            filled_field_s[ij] = ℵ[ij] * wij
        end

        coeffs_filled_s = transform(filled_field_s, S)
        coeffs_weight_s = transform(weight_field_s, S)
        (; l_indices) = coeffs_filled_s.spectrum
        for lm in eachindex(coeffs_filled_s)
            l = l_indices[lm]
            decay_s = exp(-D_s * l * (l + 1) / radius^2 * Δt_seconds)
            coeffs_filled_s[lm] *= decay_s
            coeffs_weight_s[lm] *= decay_s
        end
        filled_diffused_s = transform(coeffs_filled_s, S)
        weight_diffused_s = transform(coeffs_weight_s, S)

        min_weight_s = Float32(sea_ice_model.diffusion_min_weight_south)
        for ij in 1:npoints
            if is_antarctic[ij]
                wd_s = weight_diffused_s[ij]
                if wd_s > min_weight_s
                    ℵ[ij] = clamp(filled_diffused_s[ij] / wd_s, 0f0, 1f0)
                end
            end
        end
    end

    # reset accumulators for the next cycle
    fill!(accum_rsds, 0.0); fill!(accum_lw, 0.0); fill!(accum_wind, 0.0)
    fill!(accum_tas, 0.0); fill!(accum_tos, 0.0)
    sea_ice_model.accum_count = 0

    return nothing
end
