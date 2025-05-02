# CONVENIENCE/ALLOCATING VERSIONS

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `grids` to a newly allocated `specs::LowerTriangularArray`
using the precomputed spectral transform `S`, for CuArrays specifically."""
function SpeedyTransforms.transform(                         # GRID TO SPECTRAL
    grids::AbstractGridArray{NF, N, <:CuArray},    # input grid
    S::SpectralTransform{NF},               # precomputed spectral transform
) where {NF, N}
    ks = size(grids)[2:end]                 # the non-horizontal dimensions
    specs = CUDA.cu(zeros(LowerTriangularArray{Complex{NF}}, S.lmax+1, S.mmax+1, ks...))
    transform!(specs, grids, S)
    return specs
end


"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `specs` to a newly allocated `grids::AbstractGridArray`
using the precomputed spectral transform `S`, for CuArrays specifically."""
function SpeedyTransforms.transform(                             # SPECTRAL TO GRID
    specs::LowerTriangularArray{Complex{NF}, N, <:CuArray},     # input spectral coefficients
    S::SpectralTransform{NF};                   # precomputed spectral transform
    kwargs...                                   # pass on unscale_coslat=true/false(default)
) where {NF, N}
    ks = size(specs)[2:end]             # the non-horizontal dimensions
    grids = CUDA.cu(zeros(S.Grid{NF}, S.nlat_half, ks...))
    transform!(grids, specs, S; kwargs...)
    return grids
end
