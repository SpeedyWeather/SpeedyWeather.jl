export SimpleSpectralLongwave

"""
Simple Spectral Model (SSM) for clear-sky longwave radiative transfer.

Introduced by Williams (2026), *Journal of Advances in Modeling Earth Systems*,
doi:10.1029/2025MS005405. Bridges the gap between gray radiation schemes (simple
but inaccurate) and correlated-k schemes (accurate but opaque).

## Physical model

The SSM solves Schwarzschild's two-stream equations (Eqs. 1–2 in Williams 2026):

    dF↑/dτ = F↑ − πB(T)
    dF↓/dτ = πB(T) − F↓

at each of `nwavenumber` discrete wavenumbers, then integrates spectrally to obtain
broadband heating rates. The optical depth is

    τ(p) = D ∫₀ᵖ κ(ν̃, T, p') qGHG dp'/g       (Eq. 3)

where D = 1.5 is the two-stream diffusivity factor (Armstrong 1968).

## Absorption coefficient parameterization

Spectrally resolved mass absorption coefficients κ are represented analytically:

- **H₂O line absorption** – piecewise exponential fit to the rotation band (200–
  1000 cm⁻¹), vibration–rotation band (1000–1700 cm⁻¹) and combination band
  (1700–2500 cm⁻¹); see `h2o_line_kappa_ref` and Eq. 4 of Williams 2026.
- **H₂O continuum** – two gray bands separated at 1700 cm⁻¹ (Eq. 6); scaled by
  water-vapor partial pressure (self-broadening) and temperature (Eq. 8).
- **CO₂** – single Lorentzian centered at 667 cm⁻¹ (the "15 μm bending mode");
  Eq. 5, activated when `do_co2 = true`.

All reference coefficients (Table 1, Williams 2026) are given at
(T_ref, p_ref, RH_ref) = (260 K, 500 hPa, 100 %). Pressure broadening scales
H₂O line and CO₂ absorption as p/p_ref (Eqs. 7, 9). Continuum self-broadening
scales as pᵥ/pᵥ_ref × exp(σ_cont(T_ref − T)) (Eq. 8, Mlawer et al. 1997).

## Discretization

The two-stream integral is solved numerically at each grid column using a sequential
scan (identical to the `OneBandLongwave` scheme in SpeedyWeather, but repeated for
every spectral bin). Upward and downward beam integrations each use a single running
scalar, so no per-column storage proportional to `nwavenumber` is required.

## Usage

```julia
using SpeedyWeather
spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
longwave = SimpleSpectralLongwave(spectral_grid)          
model = PrimitiveWetModel(spectral_grid; longwave_radiation = longwave)
simulation = initialize!(model)
run!(simulation, period = Day(20))
```

Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SimpleSpectralLongwave{NF} <: AbstractLongwave

    # ---- Spectral grid ---------------------------------------------------
    "[OPTION] Number of equally-spaced wavenumber quadrature points"
    nwavenumber::Int = 41

    "[OPTION] Minimum wavenumber of spectral integration range [cm⁻¹]"
    @param wavenumber_min::NF = 10 (bounds = Positive,)

    "[OPTION] Maximum wavenumber of spectral integration range [cm⁻¹]"
    @param wavenumber_max::NF = 2510 (bounds = Positive,)

    # ---- H₂O line absorption (Table 1, Williams 2026) --------------------
    "[OPTION] Peak absorption of the pure-rotation band κ_rot [m² kg⁻¹]"
    @param κ_rot::NF = 37 (bounds = Positive,)

    "[OPTION] e-folding decay length of the rotation band l_rot [cm⁻¹]"
    @param l_rot::NF = 56 (bounds = Positive,)

    "[OPTION] Peak absorption of the vibration–rotation band κ_vr [m² kg⁻¹]"
    @param κ_vr::NF = 5 (bounds = Positive,)

    "[OPTION] e-folding decay length of vibration–rotation band (low-ν side) l_vr1 [cm⁻¹]"
    @param l_vr1::NF = 37 (bounds = Positive,)

    "[OPTION] e-folding decay length of vibration–rotation band (high-ν side) l_vr2 [cm⁻¹]"
    @param l_vr2::NF = 52 (bounds = Positive,)

    # ---- H₂O continuum (Table 1, Williams 2026) --------------------------
    "[OPTION] Continuum absorption below 1700 cm⁻¹ κ_cnt1 [m² kg⁻¹]"
    @param κ_cnt1::NF = 0.004 (bounds = Nonnegative,)

    "[OPTION] Continuum absorption above 1700 cm⁻¹ κ_cnt2 [m² kg⁻¹]"
    @param κ_cnt2::NF = 0.0002 (bounds = Nonnegative,)

    # ---- CO₂ (Table 1, Williams 2026) ------------------------------------
    "[OPTION] Include CO₂ longwave absorption"
    do_co2::Bool = true

    "[OPTION] CO₂ concentration [ppmv]"
    @param co2_ppmv::NF = 280 (bounds = Positive,)

    "[OPTION] Peak absorption coefficient of the CO₂ 15 μm band κ_CO₂ [m² kg⁻¹]"
    @param κ_co2::NF = 110 (bounds = Positive,)

    "[OPTION] e-folding half-width of CO₂ band l_CO₂ [cm⁻¹]"
    @param l_co2::NF = 12 (bounds = Positive,)

    "[OPTION] Centre wavenumber of CO₂ bending mode ν̃_CO₂ [cm⁻¹]"
    @param ν̃_co2::NF = 667 (bounds = Positive,)

    # ---- Physical reference values ---------------------------------------
    "[OPTION] Two-stream diffusivity factor D = 1.5 (Armstrong 1968)"
    @param diffusivity::NF = 1.5 (bounds = Positive,)

    "[OPTION] Reference pressure for pressure broadening [Pa]"
    @param p_ref::NF = 50000 (bounds = Positive,)

    "[OPTION] Reference temperature for absorption coefficient fits [K]"
    @param T_ref::NF = 260 (bounds = Positive,)

    "[OPTION] Reference saturation water-vapor pressure at T_ref [Pa], ≈ eₛ(260 K)"
    @param pv_ref::NF = 224.92 (bounds = Positive,)

    "[OPTION] Temperature-scaling exponent for continuum absorption σ [K⁻¹] (Mlawer et al. 1997)"
    @param σ_cont::NF = 0.02 (bounds = Nonnegative,)

    # ---- Surface emissivity ----------------------------------------------
    "[OPTION] Surface emissivity over ocean (1 = blackbody, per paper assumption)"
    @param emissivity_ocean::NF = 1 (bounds = 0 .. 1,)

    "[OPTION] Surface emissivity over land"
    @param emissivity_land::NF = 1 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure SimpleSpectralLongwave

"""$(TYPEDSIGNATURES)
Generator: construct `SimpleSpectralLongwave` from a `SpectralGrid`."""
SimpleSpectralLongwave(SG::SpectralGrid; kwargs...) = SimpleSpectralLongwave{SG.NF}(; kwargs...)

initialize!(::SimpleSpectralLongwave, ::PrimitiveEquation) = nothing

"""
    h2o_line_kappa_ref(ν̃, rad) -> κ  [m² kg⁻¹]

Reference H₂O line absorption coefficient at (T_ref, p_ref) = (260 K, 500 hPa).
Piecewise exponential fit to the rotation band, vibration–rotation band, and
combination band. Williams (2026), Eq. 4.
"""
@inline function h2o_line_kappa_ref(ν̃::NF, rad::SimpleSpectralLongwave{NF}) where {NF}
    (; κ_rot, l_rot, κ_vr, l_vr1, l_vr2) = rad
    if ν̃ <= 200                          # flat rotation band
        return κ_rot
    elseif ν̃ <= 1000                     # exponential decay of rotation band
        return κ_rot * exp(-(ν̃ - 200) / l_rot)
    elseif ν̃ <= 1450                     # rising vibration–rotation band
        return κ_vr * exp(-(1450 - ν̃) / l_vr1)
    elseif ν̃ <= 1700                     # flat top of vibration–rotation band
        return κ_vr
    elseif ν̃ <= 2500                     # exponential decay of combination band
        return κ_vr * exp(-(ν̃ - 1700) / l_vr2)
    else
        return zero(NF)
    end
end

"""
    co2_kappa_ref(ν̃, rad) -> κ  [m² kg⁻¹]

Reference CO₂ absorption coefficient. Lorentzian (exponential in wavenumber)
centred at the 15 μm bending mode (~667 cm⁻¹). Active only in [500, 850] cm⁻¹.
Williams (2026), Eq. 5.
"""
@inline function co2_kappa_ref(ν̃::NF, rad::SimpleSpectralLongwave{NF}) where {NF}
    (; κ_co2, l_co2, ν̃_co2) = rad
    return (ν̃ > 500 && ν̃ < 850) ? κ_co2 * exp(-abs(ν̃ - ν̃_co2) / l_co2) : zero(NF)
end

"""
    h2o_cont_kappa_ref(ν̃, rad) -> κ  [m² kg⁻¹]

Reference H₂O continuum absorption coefficient. Treated as two gray values
separated at 1700 cm⁻¹ (the main atmospheric window below vs. above).
Williams (2026), Eq. 6.
"""
@inline function h2o_cont_kappa_ref(ν̃::NF, rad::SimpleSpectralLongwave{NF}) where {NF}
    return ν̃ <= 1700 ? rad.κ_cnt1 : rad.κ_cnt2
end

"""
    planck_wavenumber(T, ν̃) -> B  [W m⁻² sr⁻¹ (cm⁻¹)⁻¹]

Spectral Planck function at temperature `T` [K] and wavenumber `ν̃` [cm⁻¹].
The two-stream flux source term is π × B(T, ν̃). Integrating π × B over all
wavenumbers recovers the Stefan–Boltzmann law σ T⁴.
"""
@inline function planck_wavenumber(T::NF, ν̃::NF) where {NF}
    # Physical constants (CODATA 2018)
    h   = NF(6.62607015e-34)   # Planck constant  [J s]
    c   = NF(2.99792458e8)     # speed of light   [m s⁻¹]
    k_B = NF(1.380649e-23)     # Boltzmann        [J K⁻¹]
    ν_m = ν̃ * 100              # cm⁻¹ → m⁻¹
    # B in W m⁻² sr⁻¹ m⁻¹, multiplied by 100 to convert per-m⁻¹ → per-cm⁻¹
    return NF(100) * 2 * h * ν_m^3 * c^2 / (exp(h * c * ν_m / (k_B * T)) - 1)
end

"""
    ssm_delta_tau(ij, k, ν̃, T, q, pₛ, rad, model) -> Δτ

Optical depth increment Δτ through model layer `k` at wavenumber `ν̃` [cm⁻¹].
Combines H₂O line absorption (Eq. 7), H₂O continuum (Eq. 8), and—if
`rad.do_co2 == true`—CO₂ (Eq. 9) from Williams (2026).

The integral uses the two-stream diffusivity factor D:

    Δτ = D · ∫_{pₖ}^{pₖ₊₁} κ(ν̃, T, p) qGHG dp/g

All three species are pressure-broadened (κ ∝ p / p_ref). The H₂O continuum
additionally scales with water-vapor partial pressure pᵥ / pᵥ_ref and an
exponential temperature dependence exp(σ_cont (T_ref − T)).
"""
@inline function ssm_delta_tau(
        ij, k, ν̃::NF,
        T, q, pₛ::NF,
        rad::SimpleSpectralLongwave{NF},
        model,
    ) where {NF}
    σ_full  = model.geometry.σ_levels_full   # sigma at layer midpoints
    σ_half  = model.geometry.σ_levels_half   # sigma at half levels (k and k+1)
    σ_thick = model.geometry.σ_levels_thick  # Δσ per layer
    g       = model.planet.gravity

    p_full_k = σ_full[k] * pₛ                  # [Pa] pressure at layer midpoint
    Δp_k     = σ_thick[k] * pₛ                 # [Pa] layer pressure thickness
    T_k      = T[ij, k]
    q_k      = q[ij, k]

    # ---- H₂O line absorption (Eq. 7): κ ∝ p / p_ref ---------------------
    κ_line = h2o_line_kappa_ref(ν̃, rad)
    Δτ_h2o_line = κ_line * (p_full_k / rad.p_ref) * q_k * Δp_k / g

    # ---- H₂O continuum (Eq. 8): self-broadening + temperature scaling ----
    # Vapor partial pressure: pᵥ = q p / (ε + (1 − ε) q), ε = Mw/Md ≈ 0.622
    p_v_k    = q_k * p_full_k / (NF(0.622) + NF(0.378) * q_k)
    κ_cont   = h2o_cont_kappa_ref(ν̃, rad)
    Δτ_h2o_cont = κ_cont * (p_v_k / rad.pv_ref) *
                  exp(rad.σ_cont * (rad.T_ref - T_k)) * q_k * Δp_k / g

    # ---- CO₂ (Eq. 9): well-mixed, pressure broadening -------------------
    # For a well-mixed gas with κ ∝ p / p_ref, the column integral from TOA
    # to pressure p yields τ = D κ q_CO₂ p² / (2 g p_ref).  The layer
    # increment is the difference of τ at the two bounding half levels:
    #   Δτ = D κ q_CO₂ (p_{k+1}² − pₖ²) / (2 g p_ref)
    Δτ_co2 = if rad.do_co2
        q_co2    = NF(rad.co2_ppmv * 1e-6 * 44 / 29)  # ppmv → kg kg⁻¹
        κ_co2_v  = co2_kappa_ref(ν̃, rad)
        p_half_k   = σ_half[k]   * pₛ
        p_half_kp1 = σ_half[k+1] * pₛ
        κ_co2_v * q_co2 * (p_half_kp1^2 - p_half_k^2) / (2 * g * rad.p_ref)
    else
        zero(NF)
    end

    return rad.diffusivity * (Δτ_h2o_line + Δτ_h2o_cont + Δτ_co2)
end

"""$(TYPEDSIGNATURES)
Column longwave radiative transfer for the Simple Spectral Model.

For each of `nwavenumber` wavenumber bins, solves Schwarzschild's equations
with a single upward and a single downward scan. Temperature tendencies and
surface diagnostics are accumulated spectrally into `vars.tendencies` and
`vars.parameterizations`.

References:
- Williams (2026), doi:10.1029/2025MS005405, Eqs. 1–9 and Table 1.
- Armstrong (1968), doi:10.1016/0022-4073(68)90052-6 (diffusivity factor D).
- Mlawer et al. (1997), doi:10.1029/97JD00237 (continuum temperature scaling).
"""
@propagate_inbounds function parameterization!(
        ij, vars,
        rad::SimpleSpectralLongwave{NF},
        model,
    ) where {NF}

    T    = vars.grid.temperature_prev   # [K]   (npoints × nlayers)
    q    = vars.grid.humidity_prev      # [kg/kg]
    dTdt = vars.tendencies.grid.temperature
    pₛ   = vars.grid.pressure_prev[ij] # [Pa] surface pressure

    σ   = model.atmosphere.stefan_boltzmann  # Stefan–Boltzmann [W m⁻² K⁻⁴]
    cₚ  = model.atmosphere.heat_capacity

    nlayers = size(T, 2)

    # ---- Surface temperature (land–sea-mask weighted average) ------------
    sst          = vars.prognostic.ocean.sea_surface_temperature[ij]
    lst          = vars.prognostic.land.soil_temperature[ij, 1]
    land_frac    = model.land_sea_mask.mask[ij]
    ϵ_ocean      = rad.emissivity_ocean
    ϵ_land       = rad.emissivity_land

    # Broadband surface upward flux [W m⁻²] (Stefan–Boltzmann, for diagnostics)
    U_sfc_ocean  = ifelse(isfinite(sst), ϵ_ocean * σ * sst^4, zero(sst))
    U_sfc_land   = ifelse(isfinite(lst), ϵ_land  * σ * lst^4, zero(lst))
    U_sfc_bb     = (1 - land_frac) * U_sfc_ocean + land_frac * U_sfc_land

    # Effective surface temperature for spectral Planck function
    # Use actual SST/LST; fall back to lowest layer if not finite.
    T_surf_ocean = ifelse(isfinite(sst), sst, T[ij, nlayers])
    T_surf_land  = ifelse(isfinite(lst), lst, T[ij, nlayers])

    # ---- Wavenumber spacing [cm⁻¹] --------------------------------------
    nν   = rad.nwavenumber
    dν̃   = (rad.wavenumber_max - rad.wavenumber_min) / (nν - 1)

    # ---- Broadband diagnostic accumulators [W m⁻²] ----------------------
    olr_sum::NF      = zero(NF)   # outgoing longwave radiation at TOA
    D_surf_sum::NF   = zero(NF)   # downwelling flux at surface

    # ---- Spectral loop ---------------------------------------------------
    for iv in 1:nν
        ν̃ = rad.wavenumber_min + (iv - 1) * dν̃    # [cm⁻¹] wavenumber of bin iv

        # Spectral surface upward flux [W m⁻² (cm⁻¹)⁻¹] × dν̃ = [W m⁻²]
        # π × B(T_sfc) is the hemispherical flux source term (two-stream approx.).
        B_sfc_ocean = ifelse(isfinite(sst), planck_wavenumber(T_surf_ocean, ν̃), zero(NF))
        B_sfc_land  = ifelse(isfinite(lst), planck_wavenumber(T_surf_land,  ν̃), zero(NF))
        # Combined spectral surface flux [W m⁻²]
        U_spec::NF  = dν̃ * NF(π) * (
            (1 - land_frac) * ϵ_ocean * B_sfc_ocean +
             land_frac      * ϵ_land  * B_sfc_land
        )

        # ==================================================================
        # UPWARD BEAM  (scan k = nlayers → 1)
        #
        # U is the upward spectral flux [W m⁻²] at the *bottom* of the
        # current layer (= top of the layer below, or surface for k=nlayers).
        # After each layer update, U is the flux at the *top* of that layer.
        #
        # Tendency accounting (same convention as OneBandLongwaveRadiativeTransfer):
        #   • U enters layer nlayers from below (surface BC).
        #   • After passing through layer k, U_new exits layer k at the top,
        #     removing energy from k and depositing it into k-1.
        # ==================================================================
        U::NF = U_spec
        # Upward surface flux enters the bottom of layer nlayers
        dTdt[ij, nlayers] += surface_flux_to_tendency(U / cₚ, pₛ, model)

        for k in nlayers:-1:1
            Δτ_k   = ssm_delta_tau(ij, k, ν̃, T, q, pₛ, rad, model)
            trans_k = exp(-Δτ_k)
            B_k     = planck_wavenumber(T[ij, k], ν̃)
            U_new::NF = U * trans_k + dν̃ * NF(π) * B_k * (1 - trans_k)

            if k > 1
                # U_new leaves layer k at the top, enters layer k-1 at the bottom
                dTdt[ij, k]     -= flux_to_tendency(U_new / cₚ, pₛ, k,   model)
                dTdt[ij, k - 1] += flux_to_tendency(U_new / cₚ, pₛ, k-1, model)
            else
                # k=1: U_new is the OLR leaving the atmosphere at TOA
                dTdt[ij, 1] -= flux_to_tendency(U_new / cₚ, pₛ, 1, model)
                olr_sum     += U_new
            end
            U = U_new
        end

        # ==================================================================
        # DOWNWARD BEAM  (scan k = 1 → nlayers)
        #
        # TOA boundary condition: no longwave incoming from space (F↓_TOA = 0).
        # D is the downward spectral flux at the *top* of the current layer.
        # After each layer, D is the flux at the *bottom*.
        # ==================================================================
        D::NF = zero(NF)

        for k in 1:(nlayers - 1)
            Δτ_k    = ssm_delta_tau(ij, k, ν̃, T, q, pₛ, rad, model)
            trans_k  = exp(-Δτ_k)
            B_k      = planck_wavenumber(T[ij, k], ν̃)
            D_new::NF = D * trans_k + dν̃ * NF(π) * B_k * (1 - trans_k)

            # D_new leaves layer k at the bottom, enters layer k+1 at the top
            dTdt[ij, k]     -= flux_to_tendency(D_new / cₚ, pₛ, k,   model)
            dTdt[ij, k + 1] += flux_to_tendency(D_new / cₚ, pₛ, k+1, model)
            D = D_new
        end

        # Surface layer: downward flux reaching the surface
        Δτ_nl     = ssm_delta_tau(ij, nlayers, ν̃, T, q, pₛ, rad, model)
        trans_nl   = exp(-Δτ_nl)
        B_nl       = planck_wavenumber(T[ij, nlayers], ν̃)
        D_surf::NF = D * trans_nl + dν̃ * NF(π) * B_nl * (1 - trans_nl)

        # Downward flux leaves the bottom of layer nlayers (absorbed at surface)
        dTdt[ij, nlayers] -= surface_flux_to_tendency(D_surf / cₚ, pₛ, model)
        D_surf_sum        += D_surf
    end # spectral loop

    # ---- Store diagnostics -----------------------------------------------
    vars.parameterizations.outgoing_longwave[ij]  = olr_sum
    vars.parameterizations.surface_longwave_down[ij] = D_surf_sum
    # Stefan–Boltzmann surface emission for consistency with other LW schemes
    vars.parameterizations.ocean.surface_longwave_up[ij] = U_sfc_ocean
    vars.parameterizations.land.surface_longwave_up[ij]  = U_sfc_land
    vars.parameterizations.surface_longwave_up[ij]       = U_sfc_bb

    return nothing
end
