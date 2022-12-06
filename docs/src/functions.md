# Function and type index

## Parameters and constants

```@docs
SpeedyWeather.Parameters
SpeedyWeather.Constants
```

## Boundaries and boundary conditions

```@docs
SpeedyWeather.Boundaries
```

## Spherical harmonic transform

```@docs
SpeedyWeather.GeoSpectral
SpeedyWeather.SpectralTransform
SpeedyWeather.spectral
SpeedyWeather.spectral!
SpeedyWeather.gridded
SpeedyWeather.gridded!
SpeedyWeather.triangular_truncation
SpeedyWeather.roundup_fft
SpeedyWeather.spectral_truncation
SpeedyWeather.spectral_truncation!
SpeedyWeather.spectral_interpolation!
SpeedyWeather.get_legendre_polynomials!
SpeedyWeather.∇²!
SpeedyWeather.∇²
SpeedyWeather.∇⁻²!
SpeedyWeather.∇⁻²
SpeedyWeather.gradient_latitude!
SpeedyWeather.gradient_latitude
SpeedyWeather.gradient_longitude!
SpeedyWeather.gradient_longitude
SpeedyWeather.divergence!
SpeedyWeather.curl!
SpeedyWeather._divergence!
SpeedyWeather.curl_div!
SpeedyWeather.UV_from_vordiv!
SpeedyWeather.UV_from_vor!
SpeedyWeather.ϵlm
SpeedyWeather.get_recursion_factors
```

## Dynamics

```@docs
SpeedyWeather.bernoulli_potential!
SpeedyWeather.volume_flux_divergence!
SpeedyWeather.vorticity_fluxes!
SpeedyWeather.vorticity_flux_curl!
SpeedyWeather.vorticity_flux_divergence!
```

## Geometry

```@docs
SpeedyWeather.Geometry
SpeedyWeather.vertical_coordinates
SpeedyWeather.GenLogisticCoefs
SpeedyWeather.generalised_logistic
```

## Time stepping

```@docs
SpeedyWeather.time_stepping!
SpeedyWeather.timestep!
SpeedyWeather.first_timesteps!
SpeedyWeather.leapfrog!
```

## Longwave radiation
```@docs
SpeedyWeather.radset!
SpeedyWeather.radlw_down!
SpeedyWeather.compute_bbe!
SpeedyWeather.radlw_up!
```

## Shortwave radiation
```@docs
SpeedyWeather.shortwave_radiation!
SpeedyWeather.solar!
SpeedyWeather.sol_oz!
SpeedyWeather.cloud!
SpeedyWeather.radsw!
```