# Hybrid sigma-pressure coordinates, part 2: the dynamical core

> Status: **planned**. Generalises the sigma-coordinate dynamical core (mass-weighted vertical
> integrals, vertical velocity, vertical advection, adiabatic conversion) to hybrid
> sigma-pressure coordinates, and fixes the nominal-σ accessors that currently make
> `SigmaPressureCoordinates` with a non-trivial transition produce `NaN` coefficients.

Date of initial draft: 2026-08-30

Base revision: 20467269399f40212fee413a545496226de999ca

## Originating prompt

> Can you identify a merged pull request from the last months that implemented the first part
> of hybrid sigma-pressure vertical coordinates?
>
> Create a plan to continue this work. There are open tasks in that pull request, compare
> against recent changes in the dynamical core and generalise the sigma coordinates to hybrid
> sigma-pressure wherever possible. Don't bother testing thoroughly locally but feel free to
> test the backwards compatibility with sigma coordinates on individual functions. Then create
> a new draft pull request, plan with opus execute with sonnet

## Revision log

- 2026-08-30: initial draft.
- 2026-08-30: confirmed the `NaN` bug empirically on `main` for `CubicSigmaPressureCoordinates(SpectralGrid())`
  — `B_half = [0, 0, 0.0049, 0.077, …]`, so the geopotential constants come out as
  `Δp_geopot_full = [NaN, 198.9, 181.2, …]` and `Δp_geopot_half = [Inf, 608.6, …]`.
  Verified fixed after the `get_σ_*` change (all coefficients finite, nominal σ identical
  across `SigmaCoordinates`, `SigmaPressureCoordinates(transition = _ -> 1)` and
  `CubicSigmaPressureCoordinates`; `Σ_k pressure_thickness_ratio = 1` at pₛ = 1e5 and 8e4).
- 2026-08-30: added task 10 — `SigmaPressureCoordinates` does not convert the A/B coefficients
  back to `spectral_grid.NF`, so a `transition` returning `Float64` silently produces a
  `Float64` coordinate inside a `Float32` model.

## Problem description

[#1137](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1137) ("Sigma-pressure
coordinates part 1", merged 2026-06-26) introduced `SigmaPressureCoordinates` and
`CubicSigmaPressureCoordinates`, a coordinate-agnostic interface (`pressure`, `pressure_half`,
`pressure_thickness`, `sigma`), and adapted `Geometry`, the grid geopotential, the
flux-to-tendency conversion, SPPT tapering and several parameterizations. It left five tasks
open:

- [ ] `vertical_integration!`, `vertical_velocity!` — Δσ-weighted integrals are baked into the
      sigma form of the continuity equation
- [ ] adiabatic conversion — `σ_lnp_A`, `σ_lnp_B` are sigma-specific
- [ ] the implicit solver — **R**, **L**, **W** depend on Δσ level spacing
- [ ] remaining parameterizations with TODO comments
- (`LinearDrag` and the surface condition have since been done)

Since then the dynamical core has been restructured substantially (#1139 scaled time step and
physics inside the dycore, #1150 dimensioned fields, #1151 fewer dynamic dispatches, #1177
1-based truncation, #1099/#1101 batched transforms), but none of these touched the
sigma-coordinate assumptions, so the five items are still open exactly as stated.

### The blocking bug found while scoping

`get_σ_full`/`get_σ_half`/`get_σ_thickness` return the **B coefficients**
(`B_full`, `B_half`, `B_thickness`) for `SigmaPressureCoordinates`, not the nominal
σ = A + B:

```julia
get_σ_half(σ::SigmaPressureCoordinates) = σ.B_half        # ← ∂p/∂pₛ, not σ
```

`Geometry` builds `σ_levels_half/full/thick` from `get_σ_half`, so for hybrid coordinates
every consumer of `geometry.σ_levels_*` silently receives ∂p/∂pₛ instead of σ. In the
pressure-dominated layers near the model top B → 0, so

- `initialize!(::Geopotential, ...)` evaluates `log(σ_levels_full/σ_levels_half) = log(0/0)`
  → `NaN` in `Δp_geopot_full`/`Δp_geopot_half`,
- `initialize!(::AdiabaticConversion, ...)` likewise → `NaN` in `σ_lnp_A`/`σ_lnp_B`,
- `HeldSuarez` precomputes `log_σ = log(0)` → `-Inf`.

This was not caught in #1137 because every dynamical test there used
`transition = _ -> 1.0`, i.e. the pure-sigma limit where A ≡ 0 and B ≡ σ.
`CubicSigmaPressureCoordinates` — the constructor the documentation advertises — is
therefore currently unusable in a running model.

## Background

### Notation

Half levels are indexed `k = 1 … nlayers+1` in the code (k = ½ … N+½ in the maths),
full levels `k = 1 … nlayers`.

```
p_{k+½} = A_{k+½} p_ref + B_{k+½} pₛ          (half levels)
p_k     = A_k p_ref     + B_k pₛ              (full levels, A/B averaged from half levels)
Δp_k    = p_{k+½} - p_{k-½} = ΔA_k p_ref + ΔB_k pₛ
σ_k     ≡ A_k + B_k,   Δσ_k ≡ ΔA_k + ΔB_k    (nominal sigma, pₛ-independent)
```

Useful identities (all follow from `A_half = σ_half(1 - transition(σ_half))`,
`B_half = σ_half transition(σ_half)` with `transition(1) = 1`, `σ_half[1] = 0`):

```
Σ_k ΔB_k = B_{N+½} - B_{½} = 1
Σ_k ΔA_k = A_{N+½} - A_{½} = 0
Σ_k Δσ_k = 1
```

Sigma coordinates are the special case A ≡ 0, B ≡ σ, so **every hybrid expression below
reduces to the existing sigma expression term by term**. That is the backwards-compatibility
argument used throughout: for `SigmaCoordinates` the code path is left literally unchanged
(dispatch on the coordinate type), and the new hybrid path is verified to agree with it
numerically when `transition = _ -> 1`.

### Three normalised layer weights

The dynamical core works with quantities divided by pₛ (it carries lnpₛ as the prognostic
variable). Three distinct weights appear, all equal to Δσ_k in the sigma limit — which is
exactly why the current code can use one array for all three:

| weight | meaning | sigma limit | pₛ-dependent? |
|---|---|---|---|
| `δ_k(pₛ) ≡ Δp_k/pₛ = ΔA_k (p_ref/pₛ) + ΔB_k` | mass of layer k per unit pₛ | Δσ_k | **yes** |
| `ΔB_k = ∂Δp_k/∂pₛ` | sensitivity of layer mass to pₛ | Δσ_k | no |
| `B_{k+½} = ∂p_{k+½}/∂pₛ` | sensitivity of interface pressure to pₛ | σ_{k+½} | no |

### The continuity equation (Simmons & Burridge 1981, §3)

Per-layer mass flux divergence, normalised by pₛ:

```
c_k ≡ ∇·(v_k Δp_k)/pₛ = δ_k D_k + ΔB_k (v_k·∇lnpₛ)
```

so that

```
∂lnpₛ/∂t = -Σ_{k=1}^{N} c_k
```

and the vertical mass flux at the interfaces (what the code calls `w`, i.e. M/pₛ):

```
w_k ≡ M_{k+½}/pₛ = B_{k+½} Σ_{r=1}^{N} c_r  -  Σ_{r=1}^{k} c_r ,   k = 1 … N-1 ;   w_N = 0
```

Both reduce to the current code when A ≡ 0.

### Adiabatic conversion (Simmons & Burridge 1981, eq. 3.13)

```
α_k = 1 - (p_{k-½}/Δp_k) ln(p_{k+½}/p_{k-½}),     α_1 = ln 2  (since p_{½} = 0)

(ω/p)_k = v_k·∇ln p_k
        - (1/Δp_k) ln(p_{k+½}/p_{k-½}) Σ_{r<k} ∇·(v_r Δp_r)
        - α_k (1/Δp_k) ∇·(v_k Δp_k)
```

In the stored, pₛ-normalised quantities (`div_sum_above[ij,k] = Σ_{r<k} δ_r D_r`,
`pres_flux_sum_above[ij,k] = Σ_{r<k} ΔB_r (v_r·∇lnpₛ)`, `pres_flux[ij,k] = v_k·∇lnpₛ`):

```
𝒜_k(pₛ) = -(1/δ_k) ln(p_{k+½}/p_{k-½})                      replaces σ_lnp_A[k],  𝒜_1 = 0
ℬ_k(pₛ) = -α_k                                              replaces σ_lnp_B[k]

(ω/p)_k = 𝒜_k (div_sum_above + pres_flux_sum_above)
        + ℬ_k (D_k + (ΔB_k/δ_k) pres_flux_k)
        + (B_k pₛ/p_k) pres_flux_k
```

Sigma limit: δ_k = Δσ_k, ΔB_k/δ_k = 1, B_k pₛ/p_k = 1 → exactly the current three terms.

### Why the geopotential and the implicit solver need no change

**Geopotential integration constants.** `Δp_geopot_full[k] = R ln(σ_{k+½}/σ_k)` is the
sigma form of `R ln(p_{k+½}/p_k)`. In the pure-pressure limit (B = 0, A = σ) the p_ref
factors cancel in the ratio and the two are *identical*; in the pure-sigma limit (A = 0,
B = σ) the pₛ factors cancel and they are again *identical*. Only in the blend zone is there
a discrepancy, of order ΔA·ΔB/σ² · (pₛ/p_ref - 1). So once `σ_levels_*` carry the nominal σ
(the fix above), the precomputed spectral geopotential constants are exact in both limits and
a good approximation in between. The exact (pₛ-dependent, grid-space) kernel already exists
from #1137 and is used for the grid geopotential.

**Implicit solver.** **R**, **L**, **U**, **W** are the linearisation of the system about a
reference state. Linearising the hybrid weights about pₛ = p_ref gives
δ_k(p_ref) = ΔA_k + ΔB_k = Δσ_k, and p_k(p_ref) = σ_k p_ref. So with nominal σ in
`σ_levels_*` the existing operators *are* the correct linearisation of the hybrid system —
no reformulation is needed, only the `get_σ_*` fix and a test that documents this. The
residual between the linearised and the exact hybrid operator is handled explicitly, as for
every other nonlinearity in the semi-implicit scheme.

This is a deliberate scope decision, recorded under "Known limitations".

## Summary of changes

The guiding principle: **`SigmaCoordinates` code paths are left byte-identical** (same
arrays, same kernels, same broadcasts, no extra `exp`/`log` per grid point). Everything new
is reached by dispatch on `::SigmaPressureCoordinates`.

### 1. `SpeedyWeather/src/dynamics/vertical_coordinates.jl`

Fix the nominal-σ accessors and add the coordinate-agnostic quantities the dycore needs.

```julia
# nominal sigma, pₛ-independent — the fix
get_σ_half(σ::SigmaPressureCoordinates)      = σ.A_half .+ σ.B_half
get_σ_full(σ::SigmaPressureCoordinates)      = σ.A_full .+ σ.B_full
get_σ_thickness(σ::SigmaPressureCoordinates) = σ.A_thickness .+ σ.B_thickness
```

New scalar accessors (each with a `SigmaCoordinates` and a `SigmaPressureCoordinates`
method, `@inline`, GPU-safe, no allocation):

| function | returns | sigma method |
|---|---|---|
| `sigma_half(k, coord)` | nominal σ at half level k | `coord.σ_half[k]` |
| `pressure_thickness_ratio(k, pₛ, coord)` | δ_k = Δp_k/pₛ | `coord.σ_thickness[k]` |
| `pressure_sensitivity(k, coord)` | ∂p_k/∂pₛ = B_k | `coord.σ_full[k]` |
| `pressure_sensitivity_half(k, coord)` | ∂p_{k-½}/∂pₛ = B_half[k] | `coord.σ_half[k]` |
| `pressure_thickness_sensitivity(k, coord)` | ∂Δp_k/∂pₛ = ΔB_k | `coord.σ_thickness[k]` |

`pressure_thickness_ratio` is written as `ΔA[k]*(p_ref/pₛ) + ΔB[k]` (not
`pressure_thickness/pₛ`) so it stays one FMA in the sigma-limit-free hybrid kernels.

Export the ones users may want (`sigma_half`); keep the sensitivity accessors internal but
documented.

### 2. `SpeedyWeather/src/dynamics/geometry.jl`

- Fix `geometry.vertical_coordinate` → `geometry.vertical_coordinates` in the two forwarding
  methods at the bottom of the file (currently a latent `ErrorException`, no method is
  reachable today).
- Add forwarding for `pressure_half`, `pressure_above`, `pressure_below`, `sigma`,
  `sigma_half`, `pressure_thickness_ratio`.
- Document in the `σ_levels_*` field docstrings that these are the **nominal** σ = A + B and
  are only the mass weights for `SigmaCoordinates`.

### 3. Grid-space surface pressure at the dynamical-core step

`vars.parameterizations.surface_pressure` is evaluated at Leapfrog step 1 (previous), the
dycore needs step 2 (current) — see `which_prognostic_step(var, ::AbstractLeapfrog, ...)`.
Add

```julia
DynamicsVariable(:surface_pressure, Grid2D(), desc = "Surface pressure at the dynamical core
    time step", units = "Pa"),
```

to `variables(::PrimitiveDry)` in `SpeedyWeather/src/models/primitive_dry.jl` (inherited by
`PrimitiveWet`), and fill it in `dynamics_tendencies!(::Variables, ::PrimitiveEquation)`
right after `pressure_gradient_flux!`:

```julia
surface_pressure_grid!(vars, model)   # no-op for SigmaCoordinates, exp(lnpₛ) for hybrid
```

Dispatch on `model.geometry.vertical_coordinates` so the sigma path pays nothing.

### 4. `vertical_integration!` (`dynamics/tendencies.jl`)

Split the single `Δσₖ` weight into the three weights of the table above.

- `u_mean_grid`, `v_mean_grid` ← weighted by `ΔB_k` (`pressure_thickness_sensitivity`)
- `div_mean_grid`, `div_sum_above` ← weighted by `δ_k` (`pressure_thickness_ratio`)
- `pres_flux_sum_above` ← unchanged, it is `u_mean_grid·∇lnpₛ`, i.e. already ΔB-weighted
- `div_mean` (spectral) ← **unchanged**, stays Δσ-weighted

Keeping the spectral `div_mean` on the constant Δσ weight is deliberate: the weight δ_k is
spatially varying and cannot be applied in spectral space, and the Δσ-weighted part is the
piece the semi-implicit **W** operator corrects. The difference is added in grid space (§5).

New: accumulate the hybrid correction

```
Ĉ(ij) ≡ Σ_k (δ_k - Δσ_k) D_k = (p_ref/pₛ - 1) Σ_k ΔA_k D_k
```

into a new `DynamicsVariable(:div_mean_correction, Grid2D(), ...)`. It is identically zero
for sigma coordinates (ΔA ≡ 0) and the field is simply never touched on that path.

Implementation: keep the existing CPU loop and the two GPU kernels for `::SigmaCoordinates`
exactly as they are (rename nothing), and add `::SigmaPressureCoordinates` methods /
kernels alongside. Dispatch through the existing
`vertical_integration!(arch, vars, geometry, time_stepping)` layer by adding the coordinate
as an argument.

### 5. `surface_pressure_grid_tendency!` (`dynamics/tendencies.jl`)

```julia
@. pres_tend_grid += u_mean_grid * dpres_dx + v_mean_grid * dpres_dy      # unchanged
@. pres_tend_grid += div_mean_correction                                  # hybrid only
```

Sign check: `surface_pressure_spectral_tendency!` computes
`pres_tend = -pres_tend - div_mean`, so a `+Ĉ` in the grid tendency becomes `-Ĉ` in the
final tendency, giving `∂lnpₛ/∂t = -D̄_Δσ - ū·∇lnpₛ - Ĉ = -Σ_k c_k`. ✔

Dispatch the second line on the coordinate type (no-op for sigma).

### 6. `vertical_velocity!` (`dynamics/tendencies.jl`)

For sigma, keep the current fused broadcast verbatim.

For hybrid, add a kernel over `(ij, k)` implementing

```
w[ij,k] = B_{k+½} (div_mean_grid[ij] + ūv̄∇lnp[ij])
        - (div_sum_above[ij,k]       + δ_k(pₛ) D[ij,k])
        - (pres_flux_sum_above[ij,k] + ΔB_k pres_flux[ij,k])
```

for `k = 1 … nlayers-1`, `w[ij,nlayers] = 0`. `B_{k+½}` is
`pressure_sensitivity_half(k+1, coord)`.

### 7. `vertical_advection!` (`dynamics/vertical_advection.jl`)

The kernel needs `1/δ_k` instead of `1/Δσ_k`. Introduce a tiny weight wrapper so the stencil
kernel stays a single method:

```julia
struct SigmaLayerThickness{V};  Δσ::V; end
struct HybridLayerThickness{C, F}; coordinate::C; surface_pressure::F; end

@inline layer_thickness(t::SigmaLayerThickness,  ij, k) = t.Δσ[k]
@inline layer_thickness(t::HybridLayerThickness, ij, k) =
    pressure_thickness_ratio(k, t.surface_pressure[ij], t.coordinate)
```

Both are `isbits` and `Adapt`-able, so the kernel signature is unchanged in shape and the
sigma path compiles to the same code as today. `vertical_advection!` constructs the
appropriate wrapper from `model.geometry.vertical_coordinates` and
`vars.dynamics.surface_pressure`.

### 8. Adiabatic conversion (`dynamics/adiabatic_conversion.jl`, `dynamics/tendencies.jl`)

Keep `AdiabaticConversion` and its precomputed `σ_lnp_A`/`σ_lnp_B` for `SigmaCoordinates`;
`initialize!` now receives correct nominal σ so it is unchanged.

Add a second `_temperature_tendency_kernel!` method (or a separate kernel) selected on the
coordinate type that computes 𝒜_k, ℬ_k and the `(B_k pₛ/p_k)` factor per grid point from
`pressure_half`, `pressure`, `pressure_thickness_ratio`,
`pressure_thickness_sensitivity`, `pressure_sensitivity`, following the formulas in
Background. Special cases:

- k = 1: `p_{½} = 0` → set 𝒜_1 = 0 and α_1 = ln 2 (Simmons & Burridge eq. 3.19), matching
  the sigma initialisation.
- guard `p_{k-½} > 0` rather than testing `k == 1`, so a non-zero model top also works.

Cost: ~2 `log` per layer per grid point on the hybrid path only. This mirrors the trade-off
already accepted for the hybrid grid geopotential in #1137 and should be called out in the
PR for review.

### 9. Parameterizations and forcings still on raw `σ_levels_*`

Now that `σ_levels_*` is nominal σ, most of these are *correct as a nominal-σ formulation*
but should go through the coordinate interface where an actual pressure or a mass weight is
meant:

| file | change |
|---|---|
| `dynamics/forcing.jl` `HeldSuarez` | `p = exp(log_pₛ + log_σ[k])` → `pressure(k, pₛ, coord)`; keep `sigma(k, coord)` for `temp_relax_freq`. Dispatch so the sigma kernel keeps the cheap `log_σ` path. |
| `dynamics/forcing.jl` `JetStreamForcing` | `σ_levels_full` → `sigma(k, coord)` (nominal, no behaviour change) |
| `parameterizations/convection.jl` Lee–Kim | `p = pₛ*σ[k]` → `pressure(k, pₛ, coord)` |
| `parameterizations/convection.jl` Betts–Miller (wet + dry) | `σ`→`pressure(k,pₛ,coord)/pₛ` where used as p/pₛ; `Δσ`→`pressure_thickness_ratio(k,pₛ,coord)` where used as a mass weight; `σ_half`→`pressure_half(k,pₛ,coord)/pₛ`. Audit each use before swapping. |
| `parameterizations/radiation/shortwave_radiation.jl` | `Δσ[k]` (ozone mass weighting) → `pressure_thickness_ratio`; `σ[k]` (ozone distribution) stays nominal `sigma(k, coord)` |
| `parameterizations/vertical_diffusion.jl` | `σ`/`σ_half` → `sigma`/`sigma_half` via the coordinate; keep the nominal-σ ∂σ(K ∂σ) operator and replace the TODO with a note that it is a nominal-σ discretisation |
| `dynamics/initial_conditions.jl` | `pₖ = σ_levels_full[k]*pres_grid[ij]` → `pressure(k, pₛ, coord)`; Jablonowski's η stays nominal σ (it is the test-case definition) |
| `output/writers/netcdf_output.jl`, `dynamics/particle_advection.jl`, `dynamics/horizontal_diffusion.jl`, `variables/set.jl` | nominal σ, no change needed — verify only |

Where a case genuinely needs a reformulation rather than a swap, leave a *specific* TODO
naming what would have to change (not the generic one).

### 10. Number format of the A/B coefficients

`SigmaPressureCoordinates(spectral_grid, σ_half; transition)` converts `σ_half` to
`spectral_grid.NF` but then builds `A_half`/`B_half` as `σ_half .* transition.(σ_half)`, so a
`transition` returning `Float64` (e.g. the default `σ -> σ` applied to a `Float64` literal, or
`cubic_transition` with `Float64` thresholds) yields `Float64` coefficient vectors inside a
`Float32` model. Convert `A_half`, `B_half` (and hence the derived full/thickness vectors) back
to `spectral_grid.NF` after evaluating the transition, and keep `reference_pressure` as it is.

## Testing and verification

The overriding requirement is **backwards compatibility**: `SigmaCoordinates` results must
be bit-identical to `main`. New tests go in
`SpeedyWeather/test/dynamics/vertical_coordinates.jl` (extend) and a new
`SpeedyWeather/test/dynamics/hybrid_coordinates_dycore.jl`.

1. **Accessor identities.** For `SigmaCoordinates` and for
   `SigmaPressureCoordinates(SG; transition = _ -> 1)`:
   `pressure_thickness_ratio(k, pₛ, coord) ≈ Δσ_k`,
   `pressure_sensitivity(k, coord) ≈ σ_k`,
   `pressure_thickness_sensitivity(k, coord) ≈ Δσ_k` for arbitrary pₛ.
   For `transition = _ -> 0` (pure pressure): `pressure_thickness_ratio ≈ Δσ_k p_ref/pₛ`,
   sensitivities ≈ 0.
2. **Nominal σ regression (the bug).** `get_σ_full/half/thickness` of a
   `CubicSigmaPressureCoordinates` equal the `SigmaCoordinates` values for the same `σ_half`;
   `all(isfinite, model.geopotential.Δp_geopot_full)`,
   `all(isfinite, model.adiabatic_conversion.σ_lnp_A)` after `initialize!`.
3. **Weight sums.** `Σ_k pressure_thickness_ratio(k, pₛ, coord) ≈ 1` for any pₛ (mass
   consistency), `Σ_k pressure_thickness_sensitivity ≈ 1`, `Σ_k ΔA_k ≈ 0`.
4. **Per-function backwards compatibility** — run each generalised function on a
   `SigmaCoordinates` model and on a `SigmaPressureCoordinates(…; transition = _ -> 1)`
   model with identical state, and compare the outputs of
   `vertical_integration!`, `vertical_velocity!`, `vertical_advection!`,
   `temperature_grid_tendency!`, `surface_pressure_grid_tendency!` field by field
   (`≈` with tight tolerance; the sigma path must be exactly unchanged vs `main`).
5. **Mass conservation.** For a hybrid model, check that
   `Σ_k c_k` computed from the tendencies equals `-∂lnpₛ/∂t` and that `w[:, nlayers] == 0`
   and the implied `w` at the top interface is zero.
6. **Short integration.** `PrimitiveDryModel` with `CubicSigmaPressureCoordinates`,
   Held-Suarez forcing, `run!(…, period = Day(1))`, assert everything stays finite. Extend
   the existing "Held-Suarez with varying layers and sigma spacing" testset from #1137 with a
   genuine (non-degenerate) transition — this is the case that would have caught the `NaN`.
7. Existing test suites (`dynamics/`, `parameterizations/`) must pass unchanged.

Run with `--check-bounds=yes` per `CLAUDE.md`. Full-suite runs are not required for the draft
PR; per-function checks and the short integration are.

## Documentation changes

- `docs/src/vertical_coordinates.md`: add a section on what the dynamical core does with
  hybrid coordinates — the three weights, which parts are exact and which are linearised
  about pₛ = p_ref (geopotential constants, implicit operators), and the performance
  trade-off.
- `docs/src/primitiveequation.md`: replace the σ-only statements of the continuity /
  vertical velocity / adiabatic conversion equations with the general hybrid form, noting the
  sigma limit.
- Docstrings for the new accessors.
- `CHANGELOG.md` under `## Unreleased`, first bullet:
  `- Hybrid sigma-pressure coordinates in the dynamical core: mass-weighted vertical integrals, vertical velocity, vertical advection and adiabatic conversion; fix nominal σ accessors [#NNNN](…)`
  plus a second bullet for the `NaN` fix.
- Version bump: `SpeedyWeather/Project.toml` — this changes the meaning of
  `get_σ_*(::SigmaPressureCoordinates)` and adds public accessors, so a `-DEV` minor bump.

## Known limitations

- **Spectral geopotential** still uses the precomputed `Δp_geopot_*` constants, i.e. the
  ratio `ln(p_{k+½}/p_k)` evaluated at pₛ = p_ref. Exact in the pure-pressure and pure-sigma
  limits, approximate only in the blend zone. The grid geopotential (used by the
  parameterizations) *is* exact from #1137, so the two differ slightly for hybrid
  coordinates.
- **Implicit solver** operators are the linearisation about pₛ = p_ref. Correct as a
  linearisation (see Background) but the semi-implicit correction is less well conditioned
  for hybrid coordinates when pₛ departs far from p_ref.
- **Betts-Miller convection** and **`BulkRichardsonDiffusion`** remain nominal-σ
  formulations. Reference profiles and the ∂σ(K ∂σ) operator would need to be rederived in
  pressure to be fully coordinate-agnostic.
- Hybrid coordinates are measurably slower (extra `log`/`exp` per grid point in the
  geopotential, adiabatic conversion and vertical velocity, and a grid-space
  `div_mean_correction`). `SigmaCoordinates` performance is unchanged by construction.

## Future work

- Exact spectral geopotential for hybrid coordinates (compute in grid space and transform;
  costs one extra 3D transform per step, and removes the `spectral_geopotential` variable
  TODO in `geopotential.jl`).
- Reference-state-dependent implicit operators, or a two-reference-pressure linearisation.
- Pressure-based Betts-Miller and vertical diffusion.
- GPU benchmarks for the hybrid path (`SpeedyWeather/benchmark/`).
