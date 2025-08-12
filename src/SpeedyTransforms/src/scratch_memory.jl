"""
ScratchMemory holds scratch memory for the `SpectralTransform` that's used both by the Fourier and Legendre transform. Fields are
$(TYPEDFIELDS)"""
mutable struct ScratchMemory{ # mutable struct so that referencing the scratch memory in
    NF,                       # SpectralTransform creates a reference and not a copy
    ArrayComplexType,         # <: ArrayType{Complex{NF}, 3},
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
    grid::AbstractGrid,
    nlayers::Integer) where NF

    grid_type = nonparametric_type(grid)   # always use nonparametric concrete type

    # RESOLUTION PARAMETERS
    nlon_max = get_nlon_max(grid_type, grid.nlat_half)    # number of longitudes around the equator
                                            # number of longitudes per latitude ring (one hemisphere only)
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    return ScratchMemory(NF, ArrayType, nfreq_max, nlayers, grid.nlat_half)
end 

