# THE SCRATCH MEMORY INTRODUCED HERE IS TO AVOID TEMPORARAY STORAGE IN THE SPECTRAL TRANSFORM
# THIS IS SO THAT THE SPECTRAL TRANSFORM CAN BE USED WITH CONST ACTIVITY IN ENZYME

# We do this seperate struct here for purely technical reasons to avoid that we have memory 
# alias each other in the call to `_legendre!`.
struct ColumnScratchMemory{VectorComplexType}
    north::VectorComplexType        # scratch memory for vertically batched Legendre transform
    south::VectorComplexType        # scratch memory for vertically batched Legendre transform
end

Adapt.@adapt_structure ColumnScratchMemory

"""
ScratchMemory holds scratch memory for the `SpectralTransform` that's used both by the Fourier and Legendre transform. Fields are
$(TYPEDFIELDS)"""
struct ScratchMemory{ 
    ArrayComplexType,         # <: ArrayType{Complex{NF}, 3},
    VectorType,               # <: ArrayType{NF, 1},
    VectorComplexType,        # <: ArrayType{Complex{NF}, 1},
}
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    north::ArrayComplexType
    south::ArrayComplexType

    column::ColumnScratchMemory{VectorComplexType}
end 

Adapt.@adapt_structure ScratchMemory

function ScratchMemory(
    ::Type{NF},                     
    architecture::AbstractArchitecture, 
    nfreq_max::Integer, 
    nlayers::Integer, 
    nlat_half::Integer,
    nlon_max::Integer) where NF

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory_north = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)
    scratch_memory_south = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)

    # SCRATCH MEMORY COLUMNS FOR VERTICALLY BATCHED LEGENDRE TRANSFORM
    scratch_memory_column_north = zeros(Complex{NF}, nlayers)
    scratch_memory_column_south = zeros(Complex{NF}, nlayers)

    return ScratchMemory{
        array_type(architecture, Complex{NF}, 3), 
        array_type(architecture, NF, 1),
        array_type(architecture, Complex{NF}, 1),
    }(scratch_memory_north, scratch_memory_south, 
    ColumnScratchMemory{array_type(architecture, Complex{NF}, 1)}(
    scratch_memory_column_north, scratch_memory_column_south))
end 

"""$(TYPEDSIGNATURES)
Generator function for a `ScratchMemory` that holds the scratch memory for SpeedyTransforms.
"""
function ScratchMemory(
    ::Type{NF},                     
    architecture::AbstractArchitecture, 
    grid::AbstractGrid,
    nlayers::Integer) where NF

    grid_type = nonparametric_type(grid)   # always use nonparametric concrete type

    # RESOLUTION PARAMETERS
    nlon_max = get_nlon_max(grid_type, grid.nlat_half)    # number of longitudes around the equator
                                            # number of longitudes per latitude ring (one hemisphere only)
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    return ScratchMemory(NF, architecture, nfreq_max, nlayers, grid.nlat_half, nlon_max)
end 

