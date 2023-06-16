# SpeedyTransforms

SpeedyTransforms is a submodule that has been developed for SpeedyWeather.jl which is technically
independent (SpeedyWeather.jl however imports it) and can also be used without running simulations.
It is just not put into its own respective repository.

### Example transforms

```julia
julia> using SpeedyWeather
julia> alms = zeros(ComplexF64,3,3)    # spectral coefficients
julia> alms[2,2] = 1                   # only l=1,m=1 harmonic
julia> map = gridded(alms)             # convert to grid space
8×4 Matrix{Float64}:
 -0.324541  -0.600363  -0.600363  -0.324541
 -0.134429  -0.248678  -0.248678  -0.134429
  0.134429   0.248678   0.248678   0.134429
  0.324541   0.600363   0.600363   0.324541
  0.324541   0.600363   0.600363   0.324541
  0.134429   0.248678   0.248678   0.134429
 -0.134429  -0.248678  -0.248678  -0.134429
 -0.324541  -0.600363  -0.600363  -0.324541
 
julia> spectral(map)                   # back to spectral space
3×3 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im          0.0+0.0im
 0.0+0.0im  1.0+3.60727e-17im  0.0+0.0im
 0.0+0.0im  0.0+0.0im          0.0+0.0im
```

and we have successfully reobtained the ``l=m=1`` spherical harmonic.

## Functions and type index

```@autodocs
Modules = [SpeedyWeather.SpeedyTransforms]
```