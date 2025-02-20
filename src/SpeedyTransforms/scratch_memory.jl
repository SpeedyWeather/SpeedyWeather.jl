"""
ScratchMemory holds scratch memory for the `SpectralTransform` that's used both by the Fourier and Legendre transform. Fields are
$(TYPEDFIELDS)"""
struct ScratchMemory{
    NF,
    ArrayComplexType,           # <: ArrayType{Complex{NF}, 3},
}
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    north::ArrayComplexType
    south::ArrayComplexType
end 

function ScratchMemory(
    ::Type{NF},                     
    ArrayType::Type{<:AbstractArray}, 
    nfreq_max::Integer, 
    nlayers::Integer, 
    nlat_half::Integer) where NF

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory_north = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)
    scratch_memory_south = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)

    return ScratchMemory{
        NF,
        ArrayType{Complex{NF}, 3}
    }(scratch_memory_north, scratch_memory_south)
end 

"""$(TYPEDSIGNATURES)
Generator function for a `ScratchMemory` that holds the scratch memory for SpeedyTransforms.
"""
function ScratchMemory(
    ::Type{NF},                     
    ArrayType::Type{<:AbstractArray}, 
    nlat_half::Integer, 
    Grid::Type{<:AbstractGridArray},
    nlayers::Integer) where NF

    Grid = RingGrids.nonparametric_type(Grid)   # always use nonparametric concrete type

    # RESOLUTION PARAMETERS
    nlon_max = get_nlon_max(Grid, nlat_half)    # number of longitudes around the equator
                                            # number of longitudes per latitude ring (one hemisphere only)
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    return ScratchMemory(NF, ArrayType, nfreq_max, nlayers, nlat_half)
end 

