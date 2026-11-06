@testset "GPU interpolation" begin
    arch = SpeedyWeather.GPU()
    grid_in = FullClenshawGrid(50, arch)
    grid_out = HEALPixGrid(26, arch)
    interp = RingGrids.interpolator(grid_out, grid_in)

    field_in = on_architecture(arch, rand(grid_in))
    field_out = on_architecture(arch, zeros(grid_out))
    RingGrids.interpolate!(field_out, field_in, interp)

    cpu_arch = SpeedyWeather.CPU()
    grid_in_cpu = on_architecture(cpu_arch, grid_in)
    grid_out_cpu = on_architecture(cpu_arch, grid_out)
    field_in_cpu = on_architecture(cpu_arch, field_in)
    field_out_cpu = on_architecture(cpu_arch, field_out)
    interp_cpu = RingGrids.interpolator(grid_out_cpu, grid_in_cpu)
    RingGrids.interpolate!(field_out_cpu, field_in_cpu, interp_cpu)
    @test on_architecture(cpu_arch, field_out) ≈ field_out_cpu


    field_out_cpu = rand(grid_out_cpu)
    field_out_gpu = on_architecture(arch, field_out_cpu)

    SpeedyWeather.RingGrids.update_locator!(interp_cpu, field_out_cpu)
    SpeedyWeather.RingGrids.update_locator!(interp, field_out_gpu)

    @test interp_cpu.locator.js == on_architecture(cpu_arch, interp.locator.js)
    @test interp_cpu.locator.ij_as == on_architecture(cpu_arch, interp.locator.ij_as)
    @test interp_cpu.locator.ij_bs == on_architecture(cpu_arch, interp.locator.ij_bs)
    @test interp_cpu.locator.ij_cs == on_architecture(cpu_arch, interp.locator.ij_cs)
    @test interp_cpu.locator.ij_ds == on_architecture(cpu_arch, interp.locator.ij_ds)
    # integer indices (above) must match exactly; the float weights are computed from ring
    # latitudes, which are Float64 on CPU but Float32 on Metal, so compare them approximately
    ε = sqrt(eps(Float32))
    @test interp_cpu.locator.Δys ≈ on_architecture(cpu_arch, interp.locator.Δys) rtol = ε atol = ε
    @test interp_cpu.locator.Δabs ≈ on_architecture(cpu_arch, interp.locator.Δabs) rtol = ε atol = ε
    @test interp_cpu.locator.Δcds ≈ on_architecture(cpu_arch, interp.locator.Δcds) rtol = ε atol = ε
end

# helper for the batched tests: interpolate on GPU and on CPU, return both results on CPU
function _batched_gpu_vs_cpu(field_in_cpu, grid_out_cpu)
    arch, cpu_arch = SpeedyWeather.GPU(), SpeedyWeather.CPU()
    NF = eltype(field_in_cpu)
    trailing = size(field_in_cpu)[2:end]

    field_out_cpu = zeros(NF, grid_out_cpu, trailing...)
    interp_cpu = RingGrids.interpolator(grid_out_cpu, field_in_cpu.grid, NF = NF)
    RingGrids.interpolate_2D!(field_out_cpu, field_in_cpu, interp_cpu)

    field_in = on_architecture(arch, field_in_cpu)
    grid_out = on_architecture(arch, grid_out_cpu)
    field_out = on_architecture(arch, zeros(NF, grid_out_cpu, trailing...))
    interp = RingGrids.interpolator(grid_out, field_in.grid, NF = NF)
    RingGrids.interpolate_2D!(field_out, field_in, interp)
    field_out_gpu = on_architecture(cpu_arch, field_out)

    return field_out_gpu, field_out_cpu
end

@testset "GPU batched interpolation of a 3D field" begin
    # all layers in a single kernel launch on GPU, compared against the CPU result
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        grid_in = Grid(8)
        grid_out = FullClenshawGrid(12)
        nlayers = 5
        field_in = randn(Float32, grid_in, nlayers)
        out_gpu, out_cpu = _batched_gpu_vs_cpu(field_in, grid_out)
        @test out_gpu ≈ out_cpu
    end
end

@testset "GPU batched interpolation of a 4D field" begin
    grid_in, grid_out = OctahedralGaussianGrid(8), HEALPixGrid(12)
    field_in = randn(Float32, grid_in, 4, 3)
    out_gpu, out_cpu = _batched_gpu_vs_cpu(field_in, grid_out)
    @test out_gpu ≈ out_cpu
end

@testset "GPU batched interpolation of a constant field" begin
    # constants must survive the per-layer pole averages on GPU too
    grid_in, grid_out = HEALPixGrid(8), FullGaussianGrid(16)
    nlayers = 3
    field_in = zeros(Float32, grid_in, nlayers)
    for k in 1:nlayers
        RingGrids.field_view(field_in, :, k) .= Float32(k)
    end
    out_gpu, _ = _batched_gpu_vs_cpu(field_in, grid_out)
    for k in 1:nlayers
        @test all(RingGrids.field_view(out_gpu, :, k) .≈ Float32(k))
    end
end

@testset "GPU batched interpolation of a contiguous view" begin
    # a field wrapping a contiguous view of a GPU array (e.g. one time slice) is reshaped
    # rather than looped over; the kernel then indexes a reshaped view of a GPU array
    arch, cpu_arch = SpeedyWeather.GPU(), SpeedyWeather.CPU()
    grid_in_cpu, grid_out_cpu = OctahedralGaussianGrid(8), FullGaussianGrid(12)
    nlayers, ntime = 3, 2

    backing_cpu = randn(Float32, RingGrids.get_npoints(grid_in_cpu), nlayers, ntime)
    backing = on_architecture(arch, backing_cpu)
    contiguous = view(backing, :, :, 2)
    @test RingGrids.is_flattenable(contiguous)

    grid_in = on_architecture(arch, grid_in_cpu)
    grid_out = on_architecture(arch, grid_out_cpu)
    field_in = Field(contiguous, grid_in)
    field_out = on_architecture(arch, zeros(Float32, grid_out_cpu, nlayers))
    interp = RingGrids.interpolator(grid_out, grid_in, NF = Float32)
    RingGrids.interpolate_2D!(field_out, field_in, interp)

    # reference: the same slice as a dense CPU field
    dense_cpu = Field(backing_cpu[:, :, 2], grid_in_cpu)
    out_cpu = zeros(Float32, grid_out_cpu, nlayers)
    interp_cpu = RingGrids.interpolator(grid_out_cpu, grid_in_cpu, NF = Float32)
    RingGrids.interpolate_2D!(out_cpu, dense_cpu, interp_cpu)

    @test on_architecture(cpu_arch, field_out) ≈ out_cpu
end

@testset "GPU interpolate! forwards to interpolate_2D!" begin
    arch, cpu_arch = SpeedyWeather.GPU(), SpeedyWeather.CPU()
    grid_in = HEALPixGrid(8, arch)
    grid_out = FullGaussianGrid(12, arch)
    interp = RingGrids.interpolator(grid_out, grid_in, NF = Float32)

    @testset for nlayers in (nothing, 4)
        field_in = on_architecture(arch, isnothing(nlayers) ? randn(Float32, grid_in) : randn(Float32, grid_in, nlayers))
        out1 = on_architecture(arch, isnothing(nlayers) ? zeros(Float32, grid_out) : zeros(Float32, grid_out, nlayers))
        out2 = deepcopy(out1)
        RingGrids.interpolate!(out1, field_in, interp)
        RingGrids.interpolate_2D!(out2, field_in, interp)
        @test on_architecture(cpu_arch, out1) == on_architecture(cpu_arch, out2)
    end
end
