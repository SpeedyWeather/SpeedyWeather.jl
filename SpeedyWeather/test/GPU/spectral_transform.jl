@testset "GPU spectral transform roundtrip " begin
    spectral_grid = SpectralGrid(truncation = 42, nlayers = 8, architecture = GPU())
    S = SpectralTransform(spectral_grid)

    # first roundtrip
    L = randn(ComplexF32, spectral_grid.spectrum, 8)
    field = zeros(Float32, spectral_grid.grid, 8)
    transform!(field, L, S)
    transform!(L, field, S)

    # 2nd roundtrip
    L2 = deepcopy(L)
    transform!(field, L, S)
    transform!(L, field, S)

    @test L ≈ L2
end

spectral_resolutions = (32,) #, 63, 127)
nlayers_list = (8,) # 8, 32]
# TODO: at the moment only tests for even grids (no ring on equator) pass for some reason
grid_list = [
    FullGaussianGrid,
    # FullClenshawGrid,
    OctahedralGaussianGrid,
    OctahedralClenshawGrid,
    # HEALPixGrid,
    # OctaHEALPixGrid,
]
# CUDA.cu now implicitly converts to Float32 so that's the only relevant type to
# test here
NFs = (Float32,)


# Function to generate random inputs
function get_test_data(; truncation, nlayers, Grid, NF)
    # We use dealiasing=3 to ensure that the transform is exact for both
    # Clenshaw and Gaussian grids
    spectral_grid_cpu = SpectralGrid(; NF, truncation, nlayers, Grid, dealiasing = 3)
    spectral_grid_gpu = SpectralGrid(; NF, truncation, nlayers, Grid, architecture = SpeedyWeather.GPU, dealiasing = 3)

    S_cpu = SpectralTransform(spectral_grid_cpu)
    S_gpu = SpectralTransform(spectral_grid_gpu)

    field_cpu = rand(
        NF, spectral_grid_cpu.grid,
        spectral_grid_cpu.nlayers
    )
    coeffs_cpu = rand(
        LowerTriangularArray{Complex{NF}},
        spectral_grid_cpu.spectrum,
        spectral_grid_cpu.nlayers
    )

    grid_gpu = on_architecture(S_gpu.architecture, field_cpu)
    spec_gpu = on_architecture(S_gpu.architecture, coeffs_cpu)

    return S_cpu, S_gpu, field_cpu, grid_gpu, coeffs_cpu, spec_gpu
end


@testset "Whole transform: test a round trip" begin
    @testset for NF in NFs
        @testset for truncation in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture

                    # Full return journey starting from grid
                    spec_ = SpeedyTransforms.transform(grid_cpu, S_cpu)
                    grid_cpu_roundtrip = SpeedyTransforms.transform(spec_, S_cpu)
                    # @test grid_cpu ≈ grid_test
                    # grid_cpu = grid_test

                    # Full return journey starting from spec
                    grid_ = SpeedyTransforms.transform(spec_cpu, S_cpu)
                    spec_cpu_roundtrip = SpeedyTransforms.transform(grid_, S_cpu)
                    # @test spec_cpu ≈ spec_test
                    # spec_cpu = spec_test

                    # Copy to GPU and repeat the test. We use initally
                    # transformed data so that the imaginary components of the
                    # m=1 specs are zero.
                    grid_gpu = on_architecture(gpu_arch, grid_cpu)
                    spec_gpu = on_architecture(gpu_arch, spec_cpu)
                    grid_gpu_test = on_architecture(gpu_arch, grid_cpu)
                    spec_gpu_test = on_architecture(gpu_arch, spec_cpu)

                    # Full return journey starting from grid on GPU
                    transform!(spec_gpu_test, grid_gpu, S_gpu)
                    transform!(grid_gpu_test, spec_gpu_test, S_gpu)
                    grid_test = on_architecture(cpu_arch, grid_gpu_test)
                    @test grid_cpu_roundtrip ≈ grid_test rtol = sqrt(eps(Float32))

                    # Full return journey starting from spec on GPU
                    transform!(grid_gpu_test, spec_gpu, S_gpu)
                    transform!(spec_gpu_test, grid_gpu_test, S_gpu)
                    spec_test = on_architecture(cpu_arch, spec_gpu_test)
                    @test spec_cpu_roundtrip ≈ spec_test rtol = sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "Whole transform: works in 2D" begin
    S_cpu, S_gpu, field_cpu, field_gpu, spec_cpu, spec_gpu = get_test_data(
        truncation = spectral_resolutions[1], nlayers = nlayers_list[1], Grid = grid_list[1], NF = NFs[1]
    )

    cpu_arch = S_cpu.architecture

    field_2d_cpu = field_cpu[:, 1]
    spec_2d_cpu = spec_cpu[:, 1]

    spec_2d_cpu_res = transform(field_2d_cpu, S_cpu)
    field_2d_cpu_res = transform(spec_2d_cpu, S_cpu)

    field_2d_gpu = field_gpu[:, 1]
    spec_2d_gpu = spec_gpu[:, 1]

    spec_2d_gpu_res = transform(field_2d_gpu, S_gpu)
    field_2d_gpu_res = transform(spec_2d_gpu, S_gpu)

    @test field_2d_cpu_res ≈ on_architecture(cpu_arch, field_2d_gpu_res)
    @test spec_2d_cpu_res ≈ on_architecture(cpu_arch, spec_2d_gpu_res)

    # allocating version
    # GPU
    spec_2d_gpu_res_alloc = transform(field_2d_gpu)
    field_2d_gpu_res_alloc = transform(spec_2d_gpu)

    # CPU
    spec_2d_cpu_res_alloc = transform(field_2d_cpu)
    field_2d_cpu_res_alloc = transform(spec_2d_cpu)

    @test spec_2d_cpu_res_alloc ≈ on_architecture(cpu_arch, spec_2d_gpu_res_alloc)
    @test field_2d_cpu_res_alloc ≈ on_architecture(cpu_arch, field_2d_gpu_res_alloc)
end

@testset "Whole transform: single-layer step view" begin
    # `get_step(vars.prognostic.pressure, 1)` hands `transform` a LowerTriangularMatrix whose
    # `.data` is a 1D view, a slightly different code path, so we test again seperatly here
    NF, truncation, Grid = NFs[1], spectral_resolutions[1], grid_list[1]
    spectral_grid_cpu = SpectralGrid(; NF, truncation, nlayers = 1, Grid, dealiasing = 3)
    spectral_grid_gpu = SpectralGrid(; NF, truncation, nlayers = 1, Grid, architecture = SpeedyWeather.GPU, dealiasing = 3)

    S_cpu = SpectralTransform(spectral_grid_cpu)
    S_gpu = SpectralTransform(spectral_grid_gpu)
    cpu_arch = S_cpu.architecture

    # lm x nsteps, like a prognostic surface pressure with a leapfrog step dimension
    pressure_cpu = rand(LowerTriangularArray{Complex{NF}}, spectral_grid_cpu.spectrum, 2)
    pressure_gpu = on_architecture(S_gpu.architecture, pressure_cpu)

    lnpₛ_cpu = get_step(pressure_cpu, 1)
    lnpₛ_gpu = get_step(pressure_gpu, 1)
    @test ndims(lnpₛ_gpu.data) == 1     # the 1D-data case that broke kernel compilation

    # inverse Legendre transform in isolation, so a failure here localises to the kernel
    # rather than to the Fourier stage that follows
    g_north_cpu, g_south_cpu = zero(S_cpu.scratch_memory.north), zero(S_cpu.scratch_memory.south)
    g_north_gpu, g_south_gpu = zero(S_gpu.scratch_memory.north), zero(S_gpu.scratch_memory.south)
    SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, lnpₛ_cpu, S_cpu.scratch_memory.column, S_cpu)
    SpeedyTransforms._legendre!(g_north_gpu, g_south_gpu, lnpₛ_gpu, S_gpu.scratch_memory.column, S_gpu)
    @test g_north_cpu[:, 1:1, :] ≈ on_architecture(cpu_arch, g_north_gpu)[:, 1:1, :] rtol = sqrt(eps(Float32))
    @test g_south_cpu[:, 1:1, :] ≈ on_architecture(cpu_arch, g_south_gpu)[:, 1:1, :] rtol = sqrt(eps(Float32))

    # spectral -> grid
    field_cpu = transform(lnpₛ_cpu, S_cpu)
    field_gpu = transform(lnpₛ_gpu, S_gpu)
    @test field_cpu ≈ on_architecture(cpu_arch, field_gpu) rtol = sqrt(eps(Float32))

    # and back, grid -> spectral
    spec_cpu = transform(field_cpu, S_cpu)
    spec_gpu = transform(field_gpu, S_gpu)
    @test spec_cpu ≈ on_architecture(cpu_arch, spec_gpu) rtol = sqrt(eps(Float32))
end

@testset "legendre kernels: compare to CPU" begin
    # Both GPU Legendre kernels write every output they are responsible for instead of
    # accumulating into a pre-zeroed array, so the scratch memory and the output coefficients are
    # dirtied first: anything the kernels fail to write shows up as a mismatch rather than being
    # masked by a leftover zero.
    NF = Float32
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
        @testset for nlayers in (1, 4)
            spectral_grid_cpu = SpectralGrid(; NF, trunc = 31, nlayers, Grid)
            spectral_grid_gpu = SpectralGrid(; NF, trunc = 31, nlayers, Grid, architecture = SpeedyWeather.GPU)
            S_cpu = SpectralTransform(spectral_grid_cpu)
            S_gpu = SpectralTransform(spectral_grid_gpu)
            cpu_arch = S_cpu.architecture
            nfreqs = S_cpu.nlons .÷ 2 .+ 1      # frequencies per ring the Fourier transform reads

            specs_cpu = rand(LowerTriangularArray{Complex{NF}}, spectral_grid_cpu.spectrum, nlayers)
            specs_gpu = on_architecture(S_gpu.architecture, specs_cpu)

            # INVERSE, spectral coefficients -> Fourier coefficients per ring
            g_north_cpu, g_south_cpu = zero(S_cpu.scratch_memory.north), zero(S_cpu.scratch_memory.south)
            g_north_gpu, g_south_gpu = S_gpu.scratch_memory.north, S_gpu.scratch_memory.south
            fill!(g_north_gpu, NaN)
            fill!(g_south_gpu, NaN)
            SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, specs_cpu, S_cpu.scratch_memory.column, S_cpu)
            SpeedyTransforms._legendre!(g_north_gpu, g_south_gpu, specs_gpu, S_gpu.scratch_memory.column, S_gpu)
            g_north_test = on_architecture(cpu_arch, g_north_gpu)
            g_south_test = on_architecture(cpu_arch, g_south_gpu)

            # only the frequencies the Fourier transform reads on a ring are defined, which
            # includes those past the Legendre truncation: the kernel has to zero those
            @test all(
                j -> g_north_cpu[1:nfreqs[j], 1:nlayers, j] ≈ g_north_test[1:nfreqs[j], 1:nlayers, j],
                eachindex(nfreqs)
            )
            @test all(
                j -> g_south_cpu[1:nfreqs[j], 1:nlayers, j] ≈ g_south_test[1:nfreqs[j], 1:nlayers, j],
                eachindex(nfreqs)
            )

            # FORWARD, Fourier coefficients per ring -> spectral coefficients
            field_cpu = rand(NF, spectral_grid_cpu.grid, nlayers)
            field_gpu = on_architecture(S_gpu.architecture, field_cpu)
            SpeedyTransforms._fourier!(
                S_cpu.scratch_memory.north, S_cpu.scratch_memory.south, field_cpu, S_cpu
            )
            SpeedyTransforms._fourier!(
                S_gpu.scratch_memory.north, S_gpu.scratch_memory.south, field_gpu, S_gpu
            )
            out_cpu = rand(LowerTriangularArray{Complex{NF}}, spectral_grid_cpu.spectrum, nlayers)
            out_gpu = on_architecture(S_gpu.architecture, out_cpu)      # dirty on purpose
            SpeedyTransforms._legendre!(
                out_cpu, S_cpu.scratch_memory.north, S_cpu.scratch_memory.south,
                S_cpu.scratch_memory.column, S_cpu
            )
            SpeedyTransforms._legendre!(
                out_gpu, S_gpu.scratch_memory.north, S_gpu.scratch_memory.south,
                S_gpu.scratch_memory.column, S_gpu
            )
            @test out_cpu ≈ on_architecture(cpu_arch, out_gpu) rtol = sqrt(eps(NF))
        end
    end

    # The rings contributing to an order do not always reach the equator. The coslat-dependent
    # Legendre shortcuts are not monotonic in latitude, so a ring closer to the equator can retain
    # fewer orders than one closer to the pole, and summing a whole band from the first
    # contributing ring down to the equator would pick up rings the transform has to skip. The
    # forward kernel therefore reads its rings from a table rather than assuming their extent.
    @testset "non-monotonic Legendre shortcut" begin
        NF, nlayers, LegendreShortcut = Float32, 2, SpeedyTransforms.LegendreShortcutLinCubCoslat
        spectral_grid_cpu = SpectralGrid(; NF, trunc = 31, nlayers, Grid = HEALPixGrid)
        spectral_grid_gpu = SpectralGrid(; NF, trunc = 31, nlayers, Grid = HEALPixGrid, architecture = SpeedyWeather.GPU)
        S_cpu = SpectralTransform(spectral_grid_cpu; LegendreShortcut)
        S_gpu = SpectralTransform(spectral_grid_gpu; LegendreShortcut)
        cpu_arch = S_cpu.architecture

        # guard: this really is the non-monotonic case the kernels have to cope with
        @test !issorted(S_cpu.mmax_truncation)

        specs_cpu = rand(LowerTriangularArray{Complex{NF}}, spectral_grid_cpu.spectrum, nlayers)
        specs_gpu = on_architecture(S_gpu.architecture, specs_cpu)

        field_cpu = transform(specs_cpu, S_cpu)
        field_gpu = transform(specs_gpu, S_gpu)
        @test field_cpu ≈ on_architecture(cpu_arch, field_gpu) rtol = sqrt(eps(NF))

        @test transform(field_cpu, S_cpu) ≈
            on_architecture(cpu_arch, transform(field_gpu, S_gpu)) rtol = sqrt(eps(NF))
    end
end

@testset "fourier_batched: compare forward pass to CPU" begin
    @testset for truncation in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    cpu_arch = S_cpu.architecture
                    gpu_arch = S_gpu.architecture

                    # CPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_cpu.scratch_memory.north,
                        S_cpu.scratch_memory.south,
                        grid_cpu,
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_gpu.scratch_memory.north,
                        S_gpu.scratch_memory.south,
                        grid_gpu,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)

                    # NOTE: This is probably unnecessary, we can do a coarser
                    # check for broad comparison of equivalence
                    # for ij in eachindex(grid_cpu, grid_gpu_compare)
                    # @test grid_cpu[ij] ≈ grid_gpu_compare[ij]
                    # end
                    @test grid_cpu ≈ grid_gpu_compare
                end
            end
        end
    end
end

@testset "fourier_batched: CUDA Graphs equivalence and replay" begin
    # the CUDA-Graphs acceleration lives in SpeedyTransformsCUDAExt; only test it when
    # the CUDA backend (and hence the extension) is actually loaded
    ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsCUDAExt)
    if ext !== nothing
        @testset for Grid in grid_list
            spectral_grid = SpectralGrid(; truncation = 32, nlayers = 8, Grid, architecture = GPU(), dealiasing = 3)
            field = rand(Float32, spectral_grid.grid, spectral_grid.nlayers)
            coeffs = rand(ComplexF32, spectral_grid.spectrum, spectral_grid.nlayers)

            # reference: generic (allocating) GPU path, graphs disabled on this transform
            ext.clear_fourier_graph_cache!()
            S_off = SpectralTransform(spectral_grid; gpu_graphs = false)
            spec_off = transform(field, S_off)      # grid -> spectral
            grid_off = transform(coeffs, S_off)      # spectral -> grid

            # GPU-graphs path (default, gpu_graphs = true)
            ext.clear_fourier_graph_cache!()
            S_on = SpectralTransform(spectral_grid)
            spec_on = transform(field, S_on)
            grid_on = transform(coeffs, S_on)

            # graphs path must match the generic path (forward differs only at the level
            # of the non-deterministic atomic Legendre transform, hence the tolerance)
            @test Array(spec_on.data) ≈ Array(spec_off.data) rtol = sqrt(eps(Float32))
            @test Array(grid_on.data) ≈ Array(grid_off.data) rtol = sqrt(eps(Float32))

            # a graph was actually captured and cached
            @test !isempty(ext.GRAPH_CACHES)

            # replaying into the same buffers across repeated calls stays correct
            spec_repeat = similar(spec_on)
            for _ in 1:3
                transform!(spec_repeat, field, S_on)
            end
            @test Array(spec_repeat.data) ≈ Array(spec_off.data) rtol = sqrt(eps(Float32))

            ext.clear_fourier_graph_cache!()
        end
    end
end

@testset "fourier_batched: compare backward pass to CPU" begin
    @testset for truncation in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    cpu_arch = S_cpu.architecture
                    gpu_arch = S_gpu.architecture

                    # Use scratch memory to store mid-transform data, using the
                    # CPU legendre transform to generate the intermediate data
                    # NOTE: assumption of working Legendre transform
                    g_north_cpu = S_cpu.scratch_memory.north
                    g_south_cpu = S_cpu.scratch_memory.south
                    scratch = S_cpu.scratch_memory.column
                    SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, spec_cpu, scratch, S_cpu)
                    # Copy to GPU
                    g_north_gpu = on_architecture(gpu_arch, g_north_cpu)
                    g_south_gpu = on_architecture(gpu_arch, g_south_cpu)

                    # CPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_cpu, g_north_cpu, g_south_cpu, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_gpu, g_north_gpu, g_south_gpu, S_gpu
                    )

                    # Copy back to CPU again for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)

                    # NOTE: again, not necessary to be so fine-grained here
                    # for ij in eachindex(grid_cpu, grid_gpu_compare)
                    #     @test grid_cpu[ij] ≈ grid_gpu_compare[ij]
                    # end
                    @test grid_cpu ≈ grid_gpu_compare
                end
            end
        end
    end
end


@testset "fourier_serial: compare forward pass to CPU" begin
    @testset for truncation in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture

                    # CPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_cpu.scratch_memory.north,
                        S_cpu.scratch_memory.south,
                        grid_cpu,
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_gpu.scratch_memory.north,
                        S_gpu.scratch_memory.south,
                        grid_gpu,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)

                    @test grid_cpu ≈ grid_gpu_compare
                end
            end
        end
    end
end

# NOTE: Currently failing due to problem with Float32 FFTW planning
# @testset "fourier_serial: compare backward pass to CPU" begin
#     @testset for truncation in spectral_resolutions
#         @testset for nlayers in nlayers_list
#             @testset for Grid in grid_list
#                 @testset for NF in (Float32, Float64)
#                     # Generate test data
#                     S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
#                         truncation=truncation, nlayers=nlayers, Grid=Grid, NF=NF
#                     )

#                     # Use scratch memory to store mid-transform data, using the
#                     # CPU legendre transform to generate the intermediate data
#                     # NOTE: assumption of working Legendre transform
#                     g_north_cpu = S_cpu.scratch_memory_north
#                     g_south_cpu = S_cpu.scratch_memory_south
#                     SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, spec_cpu, S_cpu)
#                     # Copy to GPU
#                     g_north_gpu = cu(g_north_cpu)
#                     g_south_gpu = cu(g_south_cpu);

#                     # CPU inverse transform
#                     SpeedyTransforms._fourier_serial!(
#                         grid_cpu, g_north_cpu, g_south_cpu, S_cpu
#                     )
#                     # GPU inverse transform
#                     SpeedyTransforms._fourier_serial!(
#                         grid_gpu, g_north_gpu, g_south_gpu, S_gpu
#                     )

#                     # Copy back to CPU again for comparison
#                     grid_gpu_compare = adapt(Array, grid_gpu)

#                     @test grid_cpu ≈ grid_gpu_compare
#                 end
#             end
#         end
#     end
# end

@testset "legendre: compare inverse transform to CPU" begin
    @testset for NF in NFs
        @testset for truncation in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture

                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_cpu.scratch_memory.north,
                        S_cpu.scratch_memory.south,
                        spec_cpu, S_cpu.scratch_memory.column, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_gpu.scratch_memory.north,
                        S_gpu.scratch_memory.south,
                        spec_gpu, S_gpu.scratch_memory.column, S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored in the
                    # scratch memory. Only the first `nlayers` columns are filled by
                    # `_legendre!`; CPU and GPU scratch may differ in capacity (dim 2),
                    # so slice to the meaningful range.
                    result_gpu = on_architecture(cpu_arch, S_gpu.scratch_memory.north)
                    result_cpu = S_cpu.scratch_memory.north
                    @test result_cpu[:, 1:nlayers, :] ≈ result_gpu[:, 1:nlayers, :] rtol = sqrt(eps(Float32))   # GPU error tolerance always Float32

                    result_gpu = on_architecture(cpu_arch, S_gpu.scratch_memory.south)
                    result_cpu = S_cpu.scratch_memory.south
                    @test result_cpu[:, 1:nlayers, :] ≈ result_gpu[:, 1:nlayers, :] rtol = sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "legendre: compare forward transform to CPU" begin
    @testset for NF in NFs
        @testset for truncation in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        truncation = truncation, nlayers = nlayers, Grid = Grid, NF = NF
                    )

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture

                    # Use scratch memory to store mid-transform data, using the
                    # CPU fourier transform to generate the intermediate data
                    f_north_cpu = S_cpu.scratch_memory.north
                    f_south_cpu = S_cpu.scratch_memory.south
                    SpeedyTransforms._fourier!(f_north_cpu, f_south_cpu, grid_cpu, S_cpu)
                    # Copy to GPU
                    f_north_gpu = on_architecture(gpu_arch, f_north_cpu)
                    f_south_gpu = on_architecture(gpu_arch, f_south_cpu)

                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_cpu,
                        f_north_cpu,
                        f_south_cpu,
                        S_cpu.scratch_memory.column,
                        S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_gpu,
                        f_north_gpu,
                        f_south_gpu,
                        S_gpu.scratch_memory.column,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored spec
                    result_gpu = on_architecture(cpu_arch, spec_gpu)
                    result_cpu = spec_cpu
                    @test result_cpu ≈ result_gpu rtol = sqrt(eps(Float32))   # GPU error tolerance always Float32
                end
            end
        end
    end
end
