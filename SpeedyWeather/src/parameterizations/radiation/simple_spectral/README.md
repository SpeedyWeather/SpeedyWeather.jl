# Simple Spectral Model (SSM) for Longwave Radiation

*Williams (2026), JAMES, doi:10.1029/2025MS005405*

---

## Why a radiation scheme?

The atmosphere absorbs and re-emits infrared (longwave) radiation. This is the greenhouse effect: gases like H₂O and CO₂ trap upwelling thermal emission from the surface, warming the lower atmosphere and cooling the upper atmosphere. Without a radiation scheme, a GCM has no way to cool the troposphere, set the tropopause, or respond to changing CO₂.

---

## The landscape of radiation schemes

There is a wide spectrum of schemes, ordered by complexity and accuracy:

| Scheme | What it does | Cost | Accuracy |
|--------|-------------|------|----------|
| **Gray radiation** | One fixed optical depth for the entire spectrum | Trivial | Poor — misses spectral structure |
| **SSM (this paper)** | Analytic fits at ~40 wavenumbers | Low | Good for climate |
| **Correlated-*k*** | RRTMG and friends, used in operational NWP | High | Excellent |

Gray schemes are too simple to respond to CO₂ changes or represent the atmospheric window (the 8–12 μm band where the atmosphere is nearly transparent and directly loses heat to space). Correlated-*k* schemes are accurate but involve large lookup tables and are expensive. The SSM fills the gap with a physically motivated analytic approximation that costs only slightly more than a gray scheme.

---

## The radiative transfer equation (physics)

The atmosphere is in local thermodynamic equilibrium, so each air parcel emits thermal radiation according to the **Planck function** at its temperature *T*:

$$\pi B(T, \tilde{\nu}) = \frac{2\pi h c^2 \tilde{\nu}^3}{\exp(hc\tilde{\nu}/k_B T) - 1}$$

where $\tilde{\nu}$ is wavenumber [cm⁻¹]. Integrating $\pi B$ over all wavenumbers recovers the Stefan–Boltzmann law $\sigma T^4$.

Radiation propagating through the atmosphere is attenuated by absorption and augmented by emission. For the two-stream (upward $F^\uparrow$ and downward $F^\downarrow$) fluxes, the **Schwarzschild equations** describe this:

$$\frac{dF^\uparrow}{d\tau} = F^\uparrow - \pi B(T, \tilde{\nu})$$
$$\frac{dF^\downarrow}{d\tau} = \pi B(T, \tilde{\nu}) - F^\downarrow$$

Here $\tau$ is **optical depth** — a dimensionless measure of how opaque the air column is. $\tau = 0$ means transparent; $\tau \gg 1$ means opaque. The solution layer by layer is exact and analytic: if you know the flux entering a layer and the layer's optical depth $\Delta\tau$, the flux leaving is

$$F^\uparrow_\text{out} = F^\uparrow_\text{in} \cdot e^{-\Delta\tau} + \pi B(T) \cdot (1 - e^{-\Delta\tau})$$

The first term is the attenuated incoming flux; the second is new thermal emission from the layer itself. The transmissivity $e^{-\Delta\tau}$ varies between 1 (transparent) and 0 (fully opaque).

---

## What determines optical depth?

The optical depth increment through a pressure layer $\Delta p$ is

$$\Delta\tau(\tilde{\nu}) = D \cdot \kappa(\tilde{\nu}, T, p) \cdot q \cdot \frac{\Delta p}{g}$$

where $q$ is the mass mixing ratio of the absorbing gas, $g$ is gravity, and $D = 1.5$ is the **diffusivity factor** (Armstrong 1968) that accounts for radiation traveling at all angles rather than straight up.

The mass absorption coefficient $\kappa$ is the hard part. In reality it is a jagged, line-by-line function of billions of molecular transitions. The SSM represents it with smooth analytic fits, calibrated against line-by-line calculations.

---

## The SSM absorption model

Three absorbers are represented, each fitted with simple functions. All reference coefficients are calibrated at $(T_\text{ref}, p_\text{ref}) = (260\,\text{K},\, 500\,\text{hPa})$ — a representative mid-tropospheric state.

### H₂O line absorption (rotation and vibration–rotation bands)

The water vapour spectrum has two main active regions. The SSM fits these as piecewise exponentials (Eq. 4 of the paper):

- **Pure rotation band** (< 1000 cm⁻¹, ~10–100 μm): flat at $\kappa_\text{rot} = 37$ m² kg⁻¹ below 200 cm⁻¹, then exponentially decaying with e-folding length $l_\text{rot} = 56$ cm⁻¹.
- **Vibration–rotation band** (~1000–2500 cm⁻¹, ~4–10 μm): Gaussian-shaped peak centred at 1450–1700 cm⁻¹ with $\kappa_\text{vr} = 5$ m² kg⁻¹.

Pressure broadening means denser (lower) air has stronger absorption: $\kappa \propto p/p_\text{ref}$.

### CO₂ 15 μm bending mode (Eq. 5)

The dominant CO₂ band is a Lorentzian centred at $\tilde{\nu} = 667$ cm⁻¹, active from 500–850 cm⁻¹:

$$\kappa_\text{CO_2}(\tilde{\nu}) = \kappa_\text{CO_2}^\text{peak} \exp\!\left(-\frac{|\tilde{\nu} - 667|}{l_\text{CO_2}}\right)$$

with peak $\kappa_\text{CO_2}^\text{peak} = 110$ m² kg⁻¹ and e-folding half-width $l_\text{CO_2} = 12$ cm⁻¹. Since CO₂ is well-mixed, its optical depth grows as $p^2$ (both abundance and pressure broadening increase with pressure).

### H₂O self-continuum (Eqs. 6–8)

Water vapour also absorbs weakly across the entire infrared through a dimer/continuum mechanism that is crucial in the otherwise nearly-transparent **atmospheric window** (8–12 μm, ~800–1250 cm⁻¹). The SSM treats this as two gray values separated at 1700 cm⁻¹, scaled by water-vapour partial pressure $p_v$ (self-broadening) and an exponential temperature dependence (Mlawer et al. 1997):

$$\Delta\tau_\text{cont} \propto \frac{p_v}{p_{v,\text{ref}}} \cdot e^{\sigma_\text{cont}(T_\text{ref} - T)} \cdot q \cdot \frac{\Delta p}{g}$$

---

## The algorithm

For each atmospheric column independently:

1. **Spectral discretization**: integrate over 41 evenly-spaced wavenumber bins from 10 to 2510 cm⁻¹ (spacing ~62 cm⁻¹).

2. **For each wavenumber bin** $\tilde{\nu}_i$:

   a. **Surface boundary condition**: the surface emits $\pi B(T_\text{sfc}, \tilde{\nu}_i)$ upward (blackbody, weighted by land fraction and emissivity).

   b. **Upward pass** (scan from bottom layer to top):
      - Compute $\Delta\tau_k$ from H₂O line + H₂O continuum + CO₂ opacity
      - Propagate $F^\uparrow$ through each layer using the analytic Schwarzschild formula
      - Flux exiting the top of the atmosphere is the **outgoing longwave radiation (OLR)**

   c. **Downward pass** (scan from top to bottom):
      - Top boundary condition: $F^\downarrow_\text{TOA} = 0$ (no incoming thermal radiation from space)
      - Propagate $F^\downarrow$ downward through each layer
      - Flux arriving at the surface is the **surface downwelling longwave**

   d. **Tendency accumulation**: the divergence of net flux through layer $k$, divided by $c_p$, gives a temperature tendency $\partial T/\partial t$.

3. **Spectrally integrate**: sum tendencies and surface fluxes over all bins to obtain broadband heating rates and OLR.

---

## What makes it novel

The SSM is not the first spectral scheme, but it is the first to represent H₂O *and* CO₂ with analytic functions that are:

- **Physically motivated** — the piecewise fit reflects real spectroscopic band structure
- **Transparent** — 10 tunable parameters with clear physical meaning (Table 1 of the paper)
- **Computationally cheap** — no lookup tables; arithmetic only
- **CO₂-responsive** — because the CO₂ band is resolved explicitly, the scheme produces a realistic radiative forcing (~3–4 W/m² per CO₂ doubling) without any tuning for that purpose

---

## Parameters (Table 1, Williams 2026)

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| $\kappa_\text{rot}$ | 37 | m² kg⁻¹ | Peak H₂O rotation-band absorption |
| $l_\text{rot}$ | 56 | cm⁻¹ | Rotation-band e-folding decay length |
| $\kappa_\text{vr}$ | 5 | m² kg⁻¹ | Peak H₂O vibration–rotation absorption |
| $l_\text{vr1}$ | 37 | cm⁻¹ | Vibration–rotation band low-$\nu$ e-folding |
| $l_\text{vr2}$ | 52 | cm⁻¹ | Vibration–rotation band high-$\nu$ e-folding |
| $\kappa_\text{cnt1}$ | 0.004 | m² kg⁻¹ | H₂O continuum below 1700 cm⁻¹ |
| $\kappa_\text{cnt2}$ | 0.0002 | m² kg⁻¹ | H₂O continuum above 1700 cm⁻¹ |
| $\kappa_\text{CO_2}$ | 110 | m² kg⁻¹ | Peak CO₂ absorption at 667 cm⁻¹ |
| $\tilde{\nu}_\text{CO_2}$ | 667 | cm⁻¹ | CO₂ band centre |
| $l_\text{CO_2}$ | 12 | cm⁻¹ | CO₂ band e-folding half-width |
