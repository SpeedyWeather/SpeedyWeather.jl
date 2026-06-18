# `unmasked_indices` and `copy_unmasked!`

These two functions together let you extract the unmasked subset of a ring-grid field into a plain array, work on it, and write results back.

**Convention:** `mask` is a `Bool` field where `true` = masked (excluded), `false` = unmasked (included).

---

## 1. Build the index vector once

```julia
mask = rand(Bool, grid)                   # or build your own Bool field
indices = unmasked_indices(mask)          # Vector of grid-point indices where mask == false
n = length(indices)                       # number of unmasked points
```

`indices` lives on the same device as `mask` (CPU or GPU) and is sorted.

---

## 2. Copy field → plain array (gather)

```julia
# 2D (single layer)
array = zeros(Float32, n)
copy_unmasked!(array, field, indices)     # array[i] == field[indices[i]]

# 3D (multiple layers)
array = zeros(Float32, n, nlayers)
copy_unmasked!(array, field, indices)     # array[i, k] == field[indices[i], k]
```

---

## 3. Copy plain array → field (scatter)

```julia
copy_unmasked!(field, array, indices)     # field[indices[i], k] = array[i, k]
```

Grid points not referenced by `indices` (the masked ones) are **left unchanged**.

---

## Full round-trip example

```julia
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
grid = spectral_grid.grid
mask = zeros(Bool, grid)                  # all unmasked initially
mask[1:100] .= true                       # mask the first 100 points

indices = unmasked_indices(mask)
n = length(indices)

# extract
temperature = rand(Float32, grid, 8)
subset = zeros(Float32, n, 8)
copy_unmasked!(subset, temperature, indices)

# ... modify subset ...

# write back
copy_unmasked!(temperature, subset, indices)
```
