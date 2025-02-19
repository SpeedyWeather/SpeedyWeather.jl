"""
SpeedyTransformsScratchMemory holds scratch memory for the `SpectralTransform` that's used both by the Fourier and Legendre transform. Fields are
$(TYPEDFIELDS)"""
struct SpeedyTransformsScratchMemory{
    NF,
    VectorType,                 # <: ArrayType{NF, 1},
    VectorComplexType,          # <: ArrayType{Complex{NF}, 1},
    ArrayComplexType,           # <: ArrayType{Complex{NF}, 3},
}
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    north::ArrayComplexType
    south::ArrayComplexType
    grid::VectorType                 # scratch memory with 1-stride for FFT output
    spec::VectorComplexType
    column_north::VectorComplexType  # scratch memory for vertically batched Legendre transform
    column_south::VectorComplexType  # scratch memory for vertically batched Legendre transform
end 

function SpeedyTransformsScratchMemory(
    ::Type{NF},                     
    ArrayType::Type{<:AbstractArray}, 
    nfreq_max::Integer, 
    nlayers::Integer, 
    nlat_half::Integer, 
    nlon_max::Integer) where NF

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory_north = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)
    scratch_memory_south = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)

    # SCRATCH MEMORY TO 1-STRIDE DATA FOR FFTs
    scratch_memory_grid  = zeros(NF, nlon_max*nlayers)
    scratch_memory_spec  = zeros(Complex{NF}, nfreq_max*nlayers)

    # SCRATCH MEMORY COLUMNS FOR VERTICALLY BATCHED LEGENDRE TRANSFORM
    scratch_memory_column_north = zeros(Complex{NF}, nlayers)
    scratch_memory_column_south = zeros(Complex{NF}, nlayers)

    return SpeedyTransformsScratchMemory{
        NF,
        ArrayType{NF, 1},
        ArrayType{Complex{NF}, 1}, 
        ArrayType{Complex{NF}, 3}
    }(scratch_memory_north, scratch_memory_south, scratch_memory_grid, scratch_memory_spec, scratch_memory_column_north, scratch_memory_column_south)
end 

"""$(TYPEDSIGNATURES)
Generator function for a `SpeedyTransformsScratchMemory` that holds the scratch memory for SpeedyTransforms.
"""
function SpeedyTransformsScratchMemory(
    ::Type{NF},                     
    ArrayType::Type{<:AbstractArray}, 
    nlat_half::Integer, 
    Grid::Type{<:AbstractGridArray},
    nlayers::Integer) where NF

    Grid = RingGrids.nonparametric_type(Grid)   # always use nonparametric concrete type

    # RESOLUTION PARAMETERS
    nlon_max = get_nlon_max(Grid, nlat_half)    # number of longitudes around the equator
                                            # number of longitudes per latitude ring (one hemisphere only)
    nlons = [RingGrids.get_nlon_per_ring(Grid, nlat_half, j) for j in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    return SpeedyTransformsScratchMemory(NF, ArrayType, nfreq_max, nlayers, nlat_half, nlon_max)
end 

