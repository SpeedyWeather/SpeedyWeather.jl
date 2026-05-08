# THE SCRATCH MEMORY INTRODUCED HERE IS TO AVOID TEMPORARAY STORAGE IN THE SPECTRAL TRANSFORM
# THIS IS SO THAT THE SPECTRAL TRANSFORM CAN BE USED WITH CONST ACTIVITY IN ENZYME

# We do this seperate struct here for purely technical reasons to avoid that we have memory
# alias each other in the call to `_legendre!`.
struct ColumnScratchMemory{VectorComplexType, MatrixComplexType, ArrayComplexType3D}
    north::VectorComplexType        # scratch memory for vertically batched Legendre transform
    south::VectorComplexType        # scratch memory for vertically batched Legendre transform

    # PER-THREAD scratch memory for the multi-threaded CPU Legendre transform.
    # The trailing dimension is the thread index. Shape (nlayers, nthreads).
    # On CPU nthreads = Threads.nthreads(); on GPU nthreads = 1 (unused — the
    # GPU _legendre! ignores ColumnScratchMemory).
    north_threads::MatrixComplexType    # per-thread north column buffers
    south_threads::MatrixComplexType    # per-thread south column buffers

    # PER-THREAD spec accumulators for the (forward) Legendre transform.
    # Shape (nspec, nlayers, nthreads). Each thread accumulates into its own
    # slice; results are summed into the output `specs` after the parallel
    # loop to avoid races. Unused on GPU (nthreads = 1).
    specs_threads::ArrayComplexType3D
end

Adapt.@adapt_structure ColumnScratchMemory

"""
ScratchMemory holds scratch memory for the `SpectralTransform` that's used both by the Fourier and Legendre transform. Fields are
$(TYPEDFIELDS)"""
struct ScratchMemory{
        ArrayComplexType,         # <: ArrayType{Complex{NF}, 3},
        VectorComplexType,        # <: ArrayType{Complex{NF}, 1},
        MatrixComplexType,        # <: ArrayType{Complex{NF}, 2},
        ArrayComplexType3D,       # <: ArrayType{Complex{NF}, 3} (per-thread spec accumulator)
    }
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    north::ArrayComplexType
    south::ArrayComplexType

    column::ColumnScratchMemory{VectorComplexType, MatrixComplexType, ArrayComplexType3D}
end

Adapt.@adapt_structure ScratchMemory

function Architectures.on_architecture(arch::AbstractArchitecture, s::ScratchMemory)
    return ScratchMemory(
        on_architecture(arch, s.north),
        on_architecture(arch, s.south),
        ColumnScratchMemory(
            on_architecture(arch, s.column.north),
            on_architecture(arch, s.column.south),
            on_architecture(arch, s.column.north_threads),
            on_architecture(arch, s.column.south_threads),
            on_architecture(arch, s.column.specs_threads),
        ),
    )
end

# Number of CPU threads to allocate per-thread Legendre scratch buffers for.
# On non-CPU architectures we allocate just one (unused) slice so the struct
# layout stays uniform across architectures.
_legendre_nthreads(::AbstractArchitecture) = 1
_legendre_nthreads(::Architectures.AbstractCPU) = Threads.nthreads()

function ScratchMemory(
        ::Type{NF},
        architecture::AbstractArchitecture,
        nfreq_max::Integer,
        nlayers::Integer,
        nlat_half::Integer,
        nlon_max::Integer,
        nspec::Integer,                                 # size(specs.data, 1) = nonzeros(spectrum)
    ) where {NF}

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory_north = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)
    scratch_memory_south = zeros(Complex{NF}, nfreq_max, nlayers, nlat_half)

    # SCRATCH MEMORY COLUMNS FOR VERTICALLY BATCHED LEGENDRE TRANSFORM
    scratch_memory_column_north = zeros(Complex{NF}, nlayers)
    scratch_memory_column_south = zeros(Complex{NF}, nlayers)

    # PER-THREAD SCRATCH MEMORY for multi-threaded CPU Legendre transform.
    # Trailing dim is the thread index; thread c uses view(buf, ..., c).
    nt = _legendre_nthreads(architecture)
    scratch_memory_column_north_threads = zeros(Complex{NF}, nlayers, nt)
    scratch_memory_column_south_threads = zeros(Complex{NF}, nlayers, nt)
    scratch_memory_specs_threads = zeros(Complex{NF}, nspec, nlayers, nt)

    Arr3T = array_type(architecture, Complex{NF}, 3)
    VecT = array_type(architecture, Complex{NF}, 1)
    MatT = array_type(architecture, Complex{NF}, 2)
    return ScratchMemory{Arr3T, VecT, MatT, Arr3T}(
        on_architecture(architecture, scratch_memory_north),
        on_architecture(architecture, scratch_memory_south),
        ColumnScratchMemory{VecT, MatT, Arr3T}(
            on_architecture(architecture, scratch_memory_column_north),
            on_architecture(architecture, scratch_memory_column_south),
            on_architecture(architecture, scratch_memory_column_north_threads),
            on_architecture(architecture, scratch_memory_column_south_threads),
            on_architecture(architecture, scratch_memory_specs_threads),
        )
    )
end

"""$(TYPEDSIGNATURES)
Generator function for a `ScratchMemory` that holds the scratch memory for SpeedyTransforms.
"""
function ScratchMemory(
        ::Type{NF},
        architecture::AbstractArchitecture,
        grid::AbstractGrid,
        nlayers::Integer,
        nspec::Integer,
    ) where {NF}

    grid_type = nonparametric_type(grid)   # always use nonparametric concrete type

    # RESOLUTION PARAMETERS
    nlon_max = get_nlon_max(grid_type, grid.nlat_half)  # number of longitudes around the equator
    # number of longitudes per latitude ring (one hemisphere only)
    nfreq_max = nlon_max ÷ 2 + 1                      # maximum number of fourier frequencies (real FFTs)

    return ScratchMemory(NF, architecture, nfreq_max, nlayers, grid.nlat_half, nlon_max, nspec)
end
