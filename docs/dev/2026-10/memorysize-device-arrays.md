# Report device (GPU) memory in printed sizes of model, simulation and variables

> Status: **completed**. `memorysize` counts device-array memory and is used in all pretty-printing.

Date of initial draft: 2026-10-28

Base revision: e292baf0

## Originating prompt

> Can you have a look at the memory size that's printed for model, simulation and
> simulation.variables, while they seem to work on the CPU and GPU I see large variable arrays
> only taking a few KB of memory -- clearly not right. If you found a way to fix this please open a PR

## Revision log

- Check `device_bytes` before the `isbits` early return so that immutable device-array types are counted too.

## Problem description

`show` for `Simulation`, models and `Variables` uses `Base.summarysize`, which for a `CuArray`
(or `ROCArray`, `MtlArray`, Reactant arrays) only counts the host-side wrapper, not the buffer
in device memory. On GPU the printed sizes are therefore just a few KB.

## Background

The benchmark suite already worked around this for Reactant (`device_bytes` in
`SpeedyWeather/benchmark/reactant_runner.jl`).

## Summary of changes

- `SpeedyWeatherInternals.Architectures`: new exported `memorysize(x)` =
  `Base.summarysize(x)` + bytes of all device arrays reachable from `x` (each object counted
  once), and a `device_bytes(x)` hook (0 by default), extended for `CuArray`, `ROCArray`,
  `MtlArray` and `Reactant.AbstractConcreteArray` in the respective extensions.
- `SpeedyWeather`: `show` of `Simulation`, `AbstractModel` and `Variables` use `memorysize`.

## Testing and verification

- Unit test in `SpeedyWeatherInternals/test/architectures.jl` with a mock device array.
- Manually: printed sizes on CUDA match those on CPU for the same configuration.

## Documentation changes

None.

## Known limitations

Memory of device arrays is only counted for backends with a `device_bytes` method.

## Future work

The Reactant benchmark runner could use `memorysize` instead of its own `device_bytes`.
