# Energy-conserving thermodynamic sea ice with conduction through the ice

> Status: **in progress**. Implemented and tested; opened as a draft PR on top of #1406 while the maintainers
> decide between this full version and a smaller conservation-only fix.

Date of initial draft: 2026-11-15 (as option 2 of the #1406 plan); split into its own plan on 2026-11-17

Base revision: 0de00d4365636eaf6367bbee3aeef1891a61d8c8 (#1406 after the review changes)

## Originating prompt

> Would you consider any of these changes due to bugs in SpeedyWeather that have been unreported?

> Test it first [...] it is important that we are careful and check this thoroughly.

> Rerun the tests and do the Earth runs for all four

> Please remove option 3 from the draft. I do not think it reasonable to suggest his climate model break
> conservation of energy.

> If Milan asks for option 2, then we will first apply it to our model and re-run the coupler.

> The PR for Speedy Weather, it appears, should have a verification section in the plan. I should compare
> nonlinear radiative blending values using the code's linear area-weighted effective temperature and true
> grid-cell average flux. I should also show that the method used to calculate Bulk Richardson Drag by means
> of liquid seawater temperature underneath the ice instead of the cold ice skin temperature is of little
> consequence due to the insulation factor. I should unit test mathematical invariants for open-water limit,
> cold-air sea ice limit, warm air cap and land isolation.

> Let's lift the 2.5m cap [...]

> Let's start on the option 2 update for the PR

> (Sign-off, 2026-11-16) Signed off. `PrescribedSeaIce`: air temperature, capped. Prescribed SST: nominal 50 m layer.

> Another part of verification can be refactoring of code. In this case, I think a subagent tasked with
> reviewing the code in the new PR once it is drafted would be appropriate to challenge how well it is
> documented for readability, its maintainability, extensibility, simplification, and speed.

> (2026-11-17) Yes, split it. [...] Can we complete the unfinished checks on a draft PR tonight?
> [...] Go ahead, push #1406 and use option A [the Jeevanjee ice surface temperature stays an open question]

## Revision log

- 2026-11-14: #1403 reported. The surface longwave over sea ice is emitted at the sea surface temperature,
  the water under the ice (about −1.8 °C), however cold the air above it.
- 2026-11-14: the fix was audited across all four ocean-emitting schemes (`OneBandLongwave`,
  `JeevanjeeRadiation`, and `AnalyticBandLongwave` and `ClearSkyEcCKDRadiation` from the
  NumericalRadiation extension). There is no separate sea-ice temperature, as asked for in #1403.
  #1406 was opened.
- 2026-11-14: the Earth A/B runs show a side effect with `SlabOcean` + `ThermodynamicSeaIce`. Summer sea ice
  vanishes, because the old SST emission had acted as the heat loss through the ice. #1406 offered two
  options: as is, or with a conductive heat loss through the ice. An option that does not conserve energy
  was considered and rejected.
- 2026-11-15: review of #1406. The conductive flux (option 2) is to land in the same PR, or the PR stays a
  draft. Other review points:
  - dispatch the emission temperature on `model.sea_ice`;
  - rely on the parameterization order instead of the `> 0` check;
  - rename the helper to `ocean_surface_emission_temperature`;
  - name all four schemes in the CHANGELOG;
  - add this plan;
  - comment on the test SST.

  The version is already `0.25.0-DEV` (#1262).
- 2026-11-15: option 2 ties the ice thickness to the concentration, `h(ℵ) = h_thin + (h_pack − h_thin)·ℵ`,
  instead of a fixed effective thickness. Top melt converts energy through the ice volume ℵ·h(ℵ).
- 2026-11-15: `filter!` freezing uses the same volume relation. The ocean's heat deficit becomes ice
  volume, with the exact inverse ℵ(V), and the heat that can't freeze at the cap stays in the ocean.
  `freeze_rate` is removed (main's 0.5 per K is reproduced at first freeze-up). Bottom melt and
  prescribed-SST freezing are listed as open questions.
- 2026-11-15: bottom melt takes its latent heat from the ocean. It is a basal flux γ·max(SST − T_f, 0),
  applied as a volume adjustment in `filter!` next to freezing and limited to the ice present. `melt_rate` is
  replaced by `basal_melt_coefficient` (default matches main as ℵ → 0). Energy now closes in every season.
- 2026-11-15: the cell's radiation is split locally into an ice part and an open-water part. The emission
  helper returns the flux-equivalent temperature T*, so the Jensen gap is gone. The slab ocean subtracts
  only the top-melt flux. The changes outside the sea-ice code are listed under option 2.
- 2026-11-15: the 2.5 m volume cap is lifted. The Europa prototype showed that the cap leaves the water under full
  ice supercooled (down to −6 °C), because conduction never weakens. Past full cover the ice now thickens (h = V),
  so conduction weakens as it grows. Upstream, `ThermodynamicSeaIce` would then store the ice volume V as its
  prognostic variable instead of ℵ, with ℵ = ℵ(V). It is still one variable.
- 2026-11-15: the verification section was extended:
  - the Jensen gap of the linear blend;
  - the Richardson number over ice;
  - mathematical invariants;
  - energy closure;
  - an independent code review before the final push.
- 2026-11-16: main (#1279, Byrne-O'Gorman longwave transmissivity as default) was merged into #1406. It changes
  only the transmissivity; the surface emission still goes through the helper, and the sea ice code on main is
  unchanged, so the design stands.
- 2026-11-16: prototype results (external coupled model, T42 and T213):
  - with the cap (v0.23-based, 1 year A/B against the SST-emission behaviour): more summer ice than the old
    handling in the north (minimum 20.8 against 19.6 × 10⁶ km²), about the same in the south, energy closed to
    0.07 W/m², but water under full cover supercooled to −6 °C (hence the lifted cap);
  - cap lifted (10-year T213 run, 7 years so far): no supercooling, energy residual −0.09 W/m² over the run,
    ice cover steady at 3.8–4.2 % of the ocean, and the ice volume levelling off as predicted (yearly growth
    0.038, 0.023, 0.018, 0.013, 0.012, 0.012, 0.007 m of sea-mean volume).
- 2026-11-16: plan signed off. The two open questions are decided:
  - `PrescribedSeaIce` keeps #1406's emission temperature, `(1 − ℵ)·SST + ℵ·min(SST, T_sa)`. It has no thickness
    or ice budget, so no T_i is solved for it;
  - `ThermodynamicSeaIce` over a prescribed SST grows ice from the deficit below freezing with a nominal 50 m mixed
    layer and leaves the SST unchanged, as on main.

- 2026-11-16: implementation choices that differ from or detail the plan:
  - **`sea_ice` is a core component.** On the GPU the column kernels only see `model.core_components`, so
    dispatching the emission on `model.sea_ice` (review change 2) needs `:sea_ice` there. `ThermodynamicSeaIce` is
    therefore immutable (it was `mutable`; nothing mutated it) and `Adapt`ed.
  - **The ice albedo is a per-point field** `sea_ice_albedo`, written by the albedo parameterization next to the ocean
    albedo (`albedo_ice` of `OceanSeaIceAlbedo` and `JinOceanAlbedo`, the cell's ocean albedo for other schemes),
    because the radiation kernel can't see the albedo scheme. The blended cell albedo is unchanged.
  - **The concentration stays a variable** (`sea_ice_concentration`, now without tendency, so not time stepped),
    diagnosed from the volume in `filter!` and `initialize!`. The albedo, surface fluxes, roughness and output
    read it unchanged.
  - **The ice absorbs the full downward longwave**, as the slab ocean does in these schemes, so the ice balance is
    `ϵσT_i⁴ = (1 − α_i) Rsd + Rld + F_c` (the plan wrote `ϵ Rld`).
  - **No reordering of `OneBandLongwave`:** the surface downward longwave is computed in a pre-pass (only with
    `ThermodynamicSeaIce`), bitwise equal to the downward beam, so the tendencies keep their summation order and
    runs without ice stay bit-identical to main.
  - **Top melt is the volume tendency** `dV/dt = −ℵ·surplus/(ρ_i L_f)`. Top melt beyond the ice present (V < 0 after
    the step) warms the ocean in `filter!`, so energy stays closed.
  - **Parameters:** `freezing_temperature`, `melting_temperature`, `conductivity`, `thin_ice_thickness`,
    `pack_ice_thickness`, `density`, `latent_heat_fusion` (own parameter, as the dry atmosphere has none),
    `basal_melt_coefficient`, `prescribed_ocean_heat_capacity`.
- 2026-11-17: #1414's review found that `ThermodynamicSeaIce` on main is an energy source (the heat lost under full
  ice is thrown away with the clamped ice, and melting takes no heat). The maintainers proposed keeping #1406 as
  option 1 and fixing the sea-ice conservation in a separate PR. This plan is that PR: option 2 rebuilt on top of
  #1406 after its review changes. The shared helper keeps #1406's signature extended by the surface downward
  longwave and the emissivity; `sea_ice = nothing` and `PrescribedSeaIce` keep #1406's blend.
- 2026-11-17: the Jeevanjee ice surface temperature (below) stays an open question for the maintainers, instead of
  adding a sensible-heat term to the ice surface balance now.

## Problem description

On main, `ThermodynamicSeaIce` does not conserve energy (#1414):

1. `filter!` turns heat lost below freezing into ice (`freeze_rate` per K) and resets the SST to the freezing point.
   Once the concentration reaches 1 the clamp throws the extra ice away, but the SST is still reset, so every W/m²
   the water loses under full ice is created again. That is an energy source of about 5–6 W/m² globally (#1414).
2. Melting (`dℵ/dt = −melt_rate (SST − T_f)`) takes no heat from the water, so summer ice melts for free.

#1406 makes the ice emit at a colder temperature than the water. That is right for the atmosphere, but it removes
the ~100 W/m² of excess emission that had acted as the heat loss through the ice, and with (2) the summer ice
vanishes. Real ice keeps growing in winter because the ocean loses heat by conduction through the ice to its cold
surface. SpeedyWeather has no such flux.

## Background

- **Emission temperature (#1406):**
  `T_eff = (1 − ℵ)·SST + ℵ·min(SST, T_sa)`.
  - ℵ is the sea ice concentration.
  - `T_sa` is `vars.parameterizations.surface_air_temperature`, the lowermost layer temperature lowered
    dry-adiabatically to the surface (`surface_condition.jl`).
  - The ice never emits warmer than the water under it.
  - Over open water, and without sea ice, `T_eff = SST` exactly.
- **`ThermodynamicSeaIce`** has concentration only, with no thickness.
  - Melting is a relaxation, `dℵ/dt = −m·max(SST − T_f, 0)`.
  - Freezing is an adjustment in `filter!`: `ℵ += f·max(T_f − SST, 0)`, after which the SST is restored to
    `T_f` and ℵ is clamped to [0, 1].
- **`SlabOcean`:**
  `C dSST/dt = Rsd − Rsu + Rld − Rlu − Lᵥ H − S` over the ocean part of each grid point.
  - The turbulent heat and humidity fluxes over ice are divided by `1 + ℵ/sea_ice_insulation`, where
    `sea_ice_insulation` is 0.01.
  - Radiation is not insulated.
- **Surface Richardson number:** `bulk_richardson_surface` (`drag.jl`, #1363) uses the SST, and the
  saturation humidity at the SST, for the ocean surface, also under ice. The resulting stability factor
  scales the ocean momentum, heat and humidity drag coefficients. Only heat and humidity are insulated.
- **Conductive flux through ice:** in a zero-layer (Semtner 1976) model, `F_c = k (T_f − T_i)/h`, with
  - k ≈ 2 W/m/K, the ice conductivity;
  - h, the ice thickness;
  - T_i, the ice surface temperature.

  T_i follows from the surface energy balance, and the ocean under the ice loses F_c, which grows ice at
  its base.

## Summary of changes

- **Diagnostic ice surface temperature (zero-layer):** solve the ice surface energy balance for T_i at each
  ice-covered point:

  `ε σ T_i⁴ = (1 − α_i) Rsd + Rld + F_c(T_i)` (the ice absorbs the full downward longwave, as the slab ocean does)

  where `F_c = k (T_f − T_i)/h(ℵ)`. The turbulent fluxes, already insulated, go to the open water (see the
  local radiation split below).
  - Use four Newton steps from `min(SST, T_sa)`, a fixed count, so the code stays branch-free
    for Reactant.
  - Cap T_i at 0 °C (melting).

  T_i is diagnostic, not a new prognostic variable. It replaces `min(SST, T_sa)` as the ice emission
  temperature for `ThermodynamicSeaIce`, which is where the dispatch in review change 2 pays off.
- **Ice thickness tied to the concentration (decided 2026-11-15):**

  `h(ℵ) = h_thin + (h_pack − h_thin)·ℵ`, with `h_thin` = 0.2 m and `h_pack` = 2.5 m.

  - A low concentration stands for new, thin ice, which conducts fast and grows fast. Full cover stands for
    thick pack ice, which conducts slowly. This gives the first-order behaviour of a prognostic thickness
    (growth ∝ 1/h) with no new prognostic variable.
  - The ice volume per unit ocean area is `V = ℵ·h(ℵ)`, so
    `dV/dℵ = h_thin + 2 (h_pack − h_thin)·ℵ`.
  - `h ≥ h_thin > 0` and `dV/dℵ ≥ h_thin > 0`, so neither division can be by zero, and both are smooth in ℵ.
  - V is quadratic in ℵ, so it inverts exactly, in a form with no division by `h_pack − h_thin`:
    `ℵ(V) = 2V / (h_thin + √(h_thin² + 4 (h_pack − h_thin) V))`, for 0 ≤ V ≤ h_pack, and ℵ = 1 beyond.
  - **Beyond full cover (no cap, decided 2026-11-15):** for V > h_pack the cover stays 1 and the ice thickens,
    so `h = h(ℵ) + max(V − h_pack, 0)`, which is V at full cover and continuous at V = h_pack. Conduction
    k ΔT / h then weakens as the ice grows (growth ∝ 1/h), so the growth limits itself.
  - Because ℵ = 1 can't tell a 2.5 m pack from a 4 m pack, **V becomes the prognostic variable** of
    `ThermodynamicSeaIce`, replacing ℵ, which is computed from V with `ℵ(V)` wherever it is needed. This stays one
    prognostic variable. Restarts, output, and anything that sets ℵ directly (initial conditions) set V instead,
    with `V(ℵ)` giving the equivalent volume. `PrescribedSeaIce` keeps prescribing ℵ, with h = h(ℵ).
    The denominator is at least 2 h_thin > 0. It also holds when h_pack = h_thin (a fixed thickness).
    Freezing and melting both step in V and convert back with ℵ(V), so a finite step closes energy exactly.
    A first-order step through dV/dℵ would not.
  - With h_thin = 0.2 m, k/h reaches 10 W/m²/K. The Newton steps must converge in this stiffest case too
    (tested below).
- **Ocean under the ice** loses `ℵ·F_c` instead of the ice part of the surface budget. This bottom heat loss
  freezes ice through the existing `filter!` adjustment.
- **Top melt:** when T_i is capped at 0 °C, the surplus `ℵ·(net surface flux − F_c)` melts ice:
  `V ← max(V − surplus·Δt / (ρ_i L_f), 0)`, then `ℵ = ℵ(V)`. If the ice runs out within the step, the
  rest of the surplus goes into the ocean. Energy is closed by construction:
  `atmosphere = ocean + latent heat of the ice`.
- **Freezing in `filter!` uses the same volume relation (decided 2026-11-15).** This replaces main's
  `ℵ += freeze_rate·max(T_f − SST, 0)`.
  - On main, `freeze_rate` (0.5 per K) sets ℵ with no energy accounting. Restoring the SST to T_f then adds
    heat that nothing pays for, and the clamp at ℵ = 1 drops more.
  - New: the heat the ocean is below freezing, `E = C·max(T_f − SST, 0)`, freezes ice.
    - `C = ρ_w c_w H` is `SlabOcean.heat_capacity_mixed_layer`, 2.09×10⁸ J/m²/K at 50 m.
    - `V ← V + E/(ρ_i L_f)`, then `ℵ = ℵ(V)`.
  - No cap: past h_pack the ice thickens at full cover (see above), and the SST returns to T_f. (The first draft
    stopped V at h_pack and left the unfrozen heat as water below T_f. That conserved energy but supercooled the
    water, by up to 4 K in the Europa prototype.)
  - Freezing slows as the pack thickens, which is physical, since dV/dℵ grows with ℵ.
  - Check against main: from open water, 1 K below freezing makes 0.68 m of ice (E/ρ_i L_f), giving
    ℵ ≈ 0.50. That is main's `freeze_rate` of 0.5 per K, so first freeze-up should stay close to main's.
  - `freeze_rate` is removed. ρ_i = 917 kg/m³ and L_f = 3.34×10⁵ J/kg become parameters. This changes the
    public parameters (the version is already `0.25.0-DEV`), and the CHANGELOG says so.
  - The step stays branch-free (`max`, `min`, `ifelse`, `sqrt` of a positive argument), as Reactant needs.
- **Bottom melt takes its heat from the ocean (decided 2026-11-15).** This replaces main's relaxation
  `dℵ/dt = −melt_rate·max(SST − T_f, 0)` in `sea_ice_kernel!`. That relaxation shrinks ℵ but takes no heat
  from the ocean.
  - Basal melt flux: `F_b = γ·max(SST − T_f, 0)` [W/m²], with `basal_melt_coefficient` γ.
  - It runs in `filter!` next to freezing, as an exact adjustment over the ocean step Δt:
    - `ΔV = min(F_b·Δt / (ρ_i L_f), V)`: it can't melt more ice than there is;
    - `V ← V − ΔV`, `ℵ = ℵ(V)`;
    - `SST ← SST − ρ_i L_f ΔV / C`: the ocean pays the latent heat.
  - Freezing (SST < T_f) and bottom melt (SST > T_f) can't both happen at one point, so one `filter!` pass
    handles both, branch-free.
  - Without ice (V = 0), ΔV = 0 and the SST is untouched. Once the ice is gone, the rest of the ocean's heat
    warms the water.
  - Default `γ = ρ_i L_f h_thin · melt_rate` ≈ 6.1 W/m²/K matches main's melt rate as ℵ → 0. Thicker ice
    then melts more slowly in ℵ, which is physical. Observed ice-ocean exchange (McPhee-type
    ρ_w c_w c_h u*) is faster, about 10–250 W/m²/K, so γ is a tuning parameter in the Earth A/B runs.
  - `melt_rate` is removed and `basal_melt_coefficient` added (CHANGELOG). `sea_ice_kernel!` no longer sets
    a melt tendency.
  - With freezing, bottom melt and top melt all stepping in V, the sum `C·SST − ρ_i L_f V` (ocean heat
    minus the latent heat stored in ice) changes only by the surface fluxes. Energy closes per point, in
    every season.
- **Local radiation split (decided 2026-11-15).** Each ocean cell is split into an ice part ℵ and an
  open-water part 1 − ℵ. Each part has its own albedo and emission temperature, and both see the same Rsd
  and Rld.
  - Ice part: albedo α_i (`albedo_ice`, 0.6), emission at T_i from the balance above.
  - Open-water part: the ocean albedo (0.06), emission at the SST. It takes the cell's whole turbulent heat
    and humidity flux. Those fluxes are already insulated by `1 + ℵ/sea_ice_insulation`, so they are small
    under ice, and the ice balance stays free of them.
  - Both parts use the scheme's own ocean emissivity (e.g. 0.65 in `OneBandLongwave`). Ice gets no
    separate emissivity, which keeps the reflected term (1 − ε)Rld independent of ℵ.
  - **Cell emission:** the helper returns the flux-equivalent temperature
    `T* = ((1 − ℵ) SST⁴ + ℵ T_i⁴)^(1/4)`, so εσT*⁴ is exactly the area-weighted emission of the two parts.
    This replaces the linear blend. The Jensen gap is gone, and the atmosphere's upward longwave equals the
    surfaces' emission in every scheme, including the band schemes, whose band sums close with T*.
  - **Why the slab ocean subtracts only the top-melt flux:** the cell total over the ocean is
    `(1 − ℵ)·(open-water budget) + ℵ·(ice surface budget) − turbulent`. The ice surface budget is
    `surplus − F_c`, which is zero minus F_c whenever T_i is below 0 °C. So the cell total already equals
    what the ocean is due, `(1 − ℵ)·(open-water budget) − turbulent − ℵ F_c`, plus ℵ·surplus. The slab
    ocean keeps its line and subtracts ℵ·surplus, which goes into top melt.
- **Changes outside the sea-ice code** (needed by the split):
  - `longwave_radiation.jl`: the helper (`ocean_surface_emission_temperature`) returns T* and takes T_i
    from the sea ice model through the dispatch of review change 2. `OneBandLongwave` and
    `JeevanjeeRadiation` are reordered so that the surface's downward longwave is known before the surface
    upward emission is computed. The ice balance then uses this step's Rld.
  - `ext/` (`AnalyticBandLongwave`, `ClearSkyEcCKDRadiation`): they compute the upward flux first, so they
    use the previous step's surface Rld, stored in `vars.parameterizations`. On the first step it is zero,
    and the Newton start `min(SST, T_sa)` keeps T_i bounded. Both schemes call the same helper.
  - `ocean_albedo.jl`: an accessor for the ice albedo alone, next to today's blend. The blended cell albedo
    is unchanged, so the shortwave totals stay as on main.
  - `ocean.jl` (`SlabOcean` kernel): `dsst` subtracts ℵ·surplus. With no sea ice the surplus is zero, so
    the line is bit-identical to main.
  - New diagnostic variables (`ParameterizationVariable`s): the ice surface temperature T_i, the top-melt
    flux ℵ·surplus, and the stored surface Rld for the band schemes.
  - The run order is unchanged: turbulent fluxes run before radiation, and shortwave before longwave
    (`radiation.jl`), so Rsd and the insulated turbulent fluxes are set when T_i is solved.
- **New `ThermodynamicSeaIce` parameters:**
  - `conductivity` = 2.0 W/m/K;
  - `thin_ice_thickness` = 0.2 m (h at ℵ → 0);
  - `pack_ice_thickness` = 2.5 m (h at ℵ = 1);
  - ice albedo and density as needed.

  The mean thickness at ℵ ≈ 0.57 is about 1.5 m. Tune h_thin and h_pack so that Earth's annual sea ice
  cycle stays close to main's.
- **Decided at sign-off (2026-11-16), formerly open questions:**
  - `PrescribedSeaIce` emits with `(1 − ℵ)·SST + ℵ·min(SST, T_sa)` (#1406's blend), with no T_i solve;
  - freezing under a prescribed SST (`haskey(vars.tendencies.ocean, :sea_surface_temperature)` false): there is
    no slab heat capacity, so the ice grows with a nominal 50 m mixed layer and the SST is left unchanged, as on
    main.
- **Prototype first (done, see the revision log):** the scheme was tried in an external coupled model (the project that found
  #1403), where the same missing flux shows up, before it is ported here.

## Testing and verification

All tests run on the CPU with Julia 1.13 and `--check-bounds=yes`.

### 1. Unit tests (`test/parameterizations/sea_ice.jl`, `longwave_radiation.jl`, `numericalradiation.jl`)

Run on the rebuilt branch (commit 7d7fdbf): `radiation.jl`, `longwave_radiation.jl`, `sea_ice.jl`,
`numericalradiation.jl`, `ocean_sea_ice.jl`, `albedo.jl`, `surface_fluxes.jl`, `surface_condition.jl`,
`surface_roughness.jl`, `ocean_land_time_stepping.jl`, `all_parametrizations.jl`, `variables.jl` and
`output/output_variables.jl`: 1281 tests, all pass. They cover:

- F_c = 0 at T_i = T_f.
- F_c flows upward (it cools the ocean) when T_i < T_f.
- F_c is linear in T_f − T_i.
- The limits work: k → ∞ gives T_i → T_f, and k → 0 makes the ice fully insulating.
- The ocean loses exactly ℵ·F_c. The top-melt surplus changes ℵ by its energy equivalent. The total energy
  closes, including latent heat.
- With fixed forcing, the Newton T_i matches the analytic root of the balance within tolerance, including
  at h = h_thin (the stiffest case).
- Thickness: h(0) = h_thin, h(1) = h_pack, h grows with ℵ; F_c falls as ℵ grows at fixed T_i.
- Volume: `ℵ(V(ℵ)) == ℵ` to round-off over [0, 1], also for h_pack = h_thin. A top melt changes V by
  exactly surplus·Δt/(ρ_i L_f). A full melt from ℵ = 1 uses ρ_i L_f h_pack in total.
- Freezing (`filter!`):
  - the ocean heat removed, C·(T_f − SST), equals ρ_i L_f ΔV;
  - beyond full cover V keeps growing with ℵ = 1, the SST returns to T_f (never below), and the sum
    `C·SST − ρ_i L_f V` is unchanged by `filter!`;
  - with fixed polar-night forcing, growth per unit time falls as the ice thickens beyond h_pack (self-limiting);
  - `h(ℵ(V)) + max(V − h_pack, 0) == V` at full cover, and is continuous at V = h_pack;
  - SST ≥ T_f gives no change;
  - from open water, 1 K below freezing gives ℵ ≈ 0.50 (main's `freeze_rate`), as a regression anchor.
- Bottom melt (`filter!`):
  - the ocean cools by exactly ρ_i L_f ΔV / C;
  - ΔV never exceeds V, so ℵ ≥ 0 without a clamp;
  - no ice gives no change to the SST;
  - SST ≤ T_f gives no melt;
  - as ℵ → 0, dℵ/dt matches main's `melt_rate`·(SST − T_f), as a regression anchor;
  - over a full freeze-melt cycle with no surface flux, `C·SST − ρ_i L_f V` stays constant to round-off.
- Local radiation split:
  - T* bounds: `min(SST, T_i) ≤ T* ≤ max(SST, T_i)`; ℵ = 0 gives `T* == SST` exactly, ℵ = 1 gives
    `T* == T_i`;
  - the ice-only albedo accessor returns `albedo_ice`, and the blended cell albedo is unchanged from main;
  - in `OneBandLongwave` and `JeevanjeeRadiation`, the Rld that the ice balance uses equals the Rld the
    surface receives in the same step (the reorder);
  - in the band schemes, the stored Rld equals the previous step's surface Rld, and the first step after
    `initialize!` and after a restart gives finite T_i;
  - the slab ocean's `dsst` differs from main's by exactly ℵ·surplus/C, and with no sea ice it is
    bit-identical to main.
- Earth A/B (below): the annual cycle of the sea ice extent stays close to main's (tolerance to agree in
  review), while the air over the ice is colder.

### 2. Energy closure

- In every grid point, the upward longwave the atmosphere receives over the ocean equals the `Rlu` the slab
  ocean deducts, for all four schemes.
- The area-weighted global sum of (atmosphere surface flux − ocean budget − land budget) is zero to
  round-off, each step.
- With option 2, the closure includes the conductive flux and the latent heat of the ice.
- Local split: per point, the atmosphere's net surface flux over the ocean equals the slab ocean's budget
  plus ℵ·surplus (the top-melt flux), for all four schemes. With T_i uncapped, the surplus is zero and the
  two are equal.

### 3. Earth A/B runs

- Default `PrimitiveWetModel`, T32 L8, 3 years, each of the four schemes, main against the PR.
- Reported: winter air over sea ice, Arctic and Antarctic sea ice extent (March, August, February,
  September), global surface and top-of-atmosphere energy balance.
- Before option 2, with OneBand (#1406):
  - air over ice: −19 to −23 °C on main, −35 to −43 °C with the PR;
  - Arctic ice in March/August: 25.4/13.2 on main, 23.5/0.0 with the PR (10⁶ km²);
  - Antarctic ice in February/September: 13.4/32.5 on main, 0.0/28.6 with the PR (10⁶ km²).

  Option 2 should restore the summer ice.
- With option 2 (2026-11-17), three arms against the new base 7a76967 / b06af911 (Byrne–O'Gorman default):
  main, #1406 before option 2, option 2. 12 runs, all stable, no NaNs. Last of 3 years; ice in 10⁶ km²,
  air = surface air poleward of 60° in local winter:

  | scheme | arm | Arctic Mar / Sep | Antarctic Feb / Sep | air N / S °C | SST °C (drift/yr) |
  |---|---|---|---|---|---|
  | OneBand | main | 17.2 / 0.0 | 8.2 / 15.7 | −17.3 / −18.2 | 22.0 (+1.2) |
  | OneBand | before opt. 2 | 16.5 / 0.0 | 0.3 / 15.1 | −19.9 / −20.5 | 21.7 (+1.2) |
  | OneBand | option 2 | 18.5 / 10.1 | 12.7 / 18.5 | −29.0 / −22.5 | 21.2 (+0.8) |

  - With OneBand, the default and the only Earth-like climate of the four here, option 2 restores the
    summer ice in both hemispheres (main itself now melts the Arctic out by September) and makes the winter
    air over the ice colder. Summer ice is above Earth's (Arctic ~6 in September): a tuning matter.
  - Rerun of the option 2 arm with OneBand on the rebuilt branch (on top of #1406, with the review fixes, commit
    7d7fdbf): Arctic 18.4 / 10.5, Antarctic 12.9 / 18.1, winter air −27.2 / −23.2 °C, the same as above within the
    year-to-year noise.
  - None of the arms is in equilibrium after 3 years (SST still rising); these are side-by-side comparisons.
  - Jeevanjee, AnalyticBand and ecCKD give climates far from Earth's in all three arms, main included
    (Jeevanjee: air 34 °C over a 4 °C ocean with 45–80 ×10⁶ km² of ice; AnalyticBand: SST 8 °C and cooling;
    clear-sky ecCKD: SST 33 °C). That is this setup on main, not #1406; for these schemes the runs show
    stability, not realism. Full table: `earth_ab2_summary.md` (outside the repository).
  - Open question: with Jeevanjee the option-2 ice surface falls to −60 to −71 °C under +20 °C air (upward
    longwave over ice 83 against 200 W/m² on main). The zero-layer balance has no turbulent (sensible or
    latent) heat exchange with the air, since those fluxes stay insulated, so under warm air the ice surface can
    only cool by radiation. Adding a sensible-heat term to the balance would change the signed-off design.

### 4. GPU

NVIDIA RTX 4050, CUDA, `ThermodynamicSeaIce` on half the points (random SST, surface air temperature and
concentration, volume from the concentration, ice albedo 0.6, 100 W/m² shortwave):

- the ice surface temperature differs between GPU and CPU by at most 3 × 10⁻⁵ K;
- the ocean surface longwave differs by at most 3 × 10⁻⁴ W/m² over ice and 3 × 10⁻⁵ W/m² over open water
  (relative at most 1 × 10⁻⁶, Float32 round-off; open water, where nothing changed, differs too);
- 20 full-model GPU time steps with `OneBandLongwave` and with `JeevanjeeRadiation`: temperature and ice volume
  finite, concentration in [0, 1].

Rerun after the review fixes (commit 7d7fdbf): the same results.

### 5. Benchmark

The repository's regression check (`benchmark/regression/regression.jl check`, debug sweep with `--length 2`) on
CPU (Intel Core Ultra 7 155H) and GPU (RTX 4050) for this PR (code as in commit 7d7fdbf) and the references
main (19a6ba1), the v0.24.1 release and the latest benchmarked revision (a0df41c). Geometric mean of the SYPD
ratio over the seven configurations, one run each:

| this PR against | CPU | GPU |
|---|---|---|
| main | 0.967 | 0.975 |
| v0.24.1 release | 0.832 | 0.967 |
| latest benchmarked | 0.848 | 0.974 |

- Against main, −3 % on CPU and −2.5 % on GPU: within the noise, no regression (cutoff 0.85).
- Against the release and the latest benchmarked revision the CPU ratio is just below 0.85, but main itself is at
  about 0.86 of the release on CPU in the same run, so that gap is main's, not this PR's.
- The extra work of this PR is the downward-beam pre-pass in `OneBandLongwave` (about `nlayers` extra `T⁴` per column)
  and four Newton steps at ice points.
- The script's own trust check failed on Windows: it compares the loaded package paths with `realpath(tree) * "/"`,
  but Windows paths use backslashes. The recorded package paths are those of the right worktrees, so the timings
  are valid; the ratios above are computed from the result files. No confirmation run or bisection was made.

### 6. Independent code review

A sub-agent without the development context reviewed the diff against #1406. It confirmed the energy bookkeeping
of the interactive slab ocean (`C·SST − ρL·V` changes only by the surface fluxes) and the GPU and Reactant
requirements, and found one bug, fixed:

- **Unbounded freezing over a prescribed SST.** The SST stays below freezing, so the same deficit froze again every
  ocean step, with no cap (main's concentration clamp at 1 had bounded it). Now the ice grows towards full cover but
  never thicker than `pack_ice_thickness` (thicker ice is kept). Multi-step test added.

Should-fix items applied:
- `pack_ice_thickness < thin_ice_thickness` gave a NaN in the concentration: pack is taken as at least thin (test);
- the top-melt surplus is never negative (the Newton residual stays with the ocean, which already counts the ice
  part of the cell budget);
- the conductive flux and top melt are zero, not 0·NaN, without ice (e.g. over land with a NaN SST);
- output variables `SeaIceVolume` and `SeaIceSurfaceTemperature`;
- docs: the band schemes' one-step lag (ecCKD also in the shortwave), the restart with zero stored downward
  longwave, the ice albedo without an albedo scheme, and how to set the volume after restarting from an old file;
- the docstring of `ocean_surface_emission_temperature` for `ThermodynamicSeaIce` says it writes the diagnostics;
- an always-true guard removed; this plan corrected (full Rld, four Newton steps, `conductivity`).

Not applied:
- reordering `OneBandLongwave` so the downward beam runs first (avoids the duplicate pre-pass): kept, so that the
  tendencies keep main's summation order and runs without ice stay bit-identical to main;
- reading old restart files without a volume automatically (documented instead);
- fusing the ocean mask into the broadcasts (already on main).

## Documentation changes

- `docs/src/sea_ice.md`: the ice volume, thickness and concentration, the conductive flux, the diagnostic ice
  surface temperature, freezing, bottom and top melt, the new parameters and the energy budget.
- `docs/src/radiation.md`: the flux-equivalent emission temperature with `ThermodynamicSeaIce`.
- `docs/src/ocean.md`: the slab ocean under ice loses the conductive flux and pays for top melt beyond the ice.
- CHANGELOG.

## Known limitations

- The band schemes (`AnalyticBandLongwave`, `ClearSkyEcCKDRadiation`) solve T_i with the previous step's
  Rld (zero on the first step), a one-step lag that `OneBandLongwave` and `JeevanjeeRadiation` don't have.
The ice and the open water share the scheme's ocean emissivity, and the cell's whole insulated turbulent
  flux is given to the open water.
The surface Richardson number over ice uses the SST, which overstates the momentum drag over ice. Heat
  and humidity are covered by the insulation.
`ThermodynamicSeaIce` has a single prognostic variable, the ice volume V. Up to full cover, cover and
  thickness both follow from V, so a cell can't hold thin ice at full cover (a fresh freeze-up) or thick ice at
  low cover (a broken-up pack). Beyond full cover only the thickness grows.
Brine rejection on freezing (the salt flux into the water below) is not modelled.
`T_sa` comes from a dry-adiabatic extrapolation. It is warm-biased in surface inversions.

## Open question for the maintainers

With `JeevanjeeRadiation`, the ice surface falls to −60 to −71 °C under +20 °C air (upward longwave over ice 83
against 200 W/m² on main). The zero-layer balance has no turbulent (sensible or latent) heat exchange with the air,
since those fluxes stay insulated and go to the open water, so under warm air the ice surface can only cool by
radiation. Adding a sensible-heat term to the ice surface balance would fix this but changes the signed-off design.
`JeevanjeeRadiation`'s default climate in this setup is far from Earth's on main too (air 34 °C over a 4 °C ocean).

## Future work

- An ice-skin surface Richardson number (and humidity) over sea ice.
- Separate prognostic cover and thickness (two variables), so a cell can hold thin ice at full cover or thick
  ice at low cover.
- Snow on sea ice (insulation, albedo).
