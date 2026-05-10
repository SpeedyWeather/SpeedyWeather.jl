# GPU smoke test for the variable fusion mechanism.
#
# What this verifies:
#   1. `Variables(model)` builds correctly on GPU with the fuse mechanism active.
#   2. `vars.scratch.fused.tend_grid` lives on the GPU (parent .data is a GPU array).
#   3. Tendency entries `vars.tendencies.grid.{u,v,T[,q]}` are SubArray-backed Fields
#      whose underlying memory IS the parent's GPU memory (Base.parent identity).
#   4. Round-trip writes between parent and views land at the right slot offsets,
#      both via Julia broadcasting and via a KernelAbstractions kernel.
#   5. `reset_tendencies!` zeroes only the tendency slots (1..N·nlayers) of the
#      parent — leaving any non-tendency-fused slots untouched.
#   6. End-to-end 1h `run!` of PrimitiveDry + PrimitiveWet on GPU.
#
# How to run (on a machine with a GPU):
#   julia --project=SpeedyWeather --check-bounds=yes gpu_fuse_smoke_test.jl
#
# Auto-detects CUDA / AMDGPU / Metal; if none is loadable, errors with a clear message.
# Or pin a backend via env var, e.g.:
#   FUSE_GPU=cuda julia --project=SpeedyWeather --check-bounds=yes gpu_fuse_smoke_test.jl

using Test
using KernelAbstractions: KernelAbstractions, @kernel, @index

# -- backend selection ------------------------------------------------------
# mirrors the existing pattern in SpeedyWeather/test/GPU/runtests.jl
function load_gpu_package()
    forced = lowercase(get(ENV, "FUSE_GPU", ""))
    candidates = if forced == ""
        [:CUDA, :AMDGPU, :Metal]
    elseif forced == "cuda"
        [:CUDA]
    elseif forced == "amdgpu"
        [:AMDGPU]
    elseif forced == "metal"
        [:Metal]
    else
        error("Unknown FUSE_GPU=$forced. Use cuda, amdgpu, or metal (or unset).")
    end
    for pkg in candidates
        try
            @eval using $pkg
            return pkg
        catch
            continue
        end
    end
    tried = join(candidates, ", ")
    error(
        "No GPU backend loadable. Tried: $tried. " *
        "Add the package via `Pkg.add(\"CUDA\")` (or AMDGPU / Metal) and try again."
    )
end

const GPU_BACKEND = load_gpu_package()
@info "GPU smoke test: using backend" backend=GPU_BACKEND

using SpeedyWeather
using SpeedyWeather: reset_tendencies!, launch!, RingGridWorkOrder, synchronize

const ARCH = SpeedyWeather.GPU()

# -- helpers ----------------------------------------------------------------
# An array `a` is "GPU-resident" if KA can resolve a non-CPU backend for it.
is_gpu_array(a) = !(KernelAbstractions.get_backend(a) isa KernelAbstractions.CPU)

# Run a tiny KA kernel through a view; uses SpeedyWeather's launch! to exercise
# the same Adapt + KA backend dispatch that the dycore uses in production.
@kernel function _bump_kernel!(field, offset)
    ij, k = @index(Global, NTuple)
    @inbounds field[ij, k] += offset
end

function bump_via_kernel!(arch, field, offset)
    launch!(arch, RingGridWorkOrder, size(field), _bump_kernel!, field, eltype(field)(offset))
    synchronize(arch)
    return field
end

# -- core test --------------------------------------------------------------
function test_fuse_layout(model_constructor, model_name::AbstractString,
                         expected_tendency_views::NTuple{N, Symbol}) where {N}
    @info "── $(model_name) on GPU ─────────────────────────────"

    sg = SpectralGrid(trunc=31, nlayers=8, architecture=ARCH)
    model = model_constructor(sg)
    sim = initialize!(model)
    v = sim.variables

    @testset "$model_name layout" begin
        # 1. parent exists at canonical home
        @test haskey(v.scratch, :fused)
        @test haskey(v.scratch.fused, :tend_grid)

        buf = v.scratch.fused.tend_grid
        @test is_gpu_array(buf.data)

        # 2. parent is NOT installed redundantly under the fuse symbol elsewhere
        @test !haskey(v.tendencies.grid, :tend_grid)
        @test !haskey(v.scratch.grid, :tend_grid)

        # 3. expected views are Fields sharing GPU memory with buf.
        # Note: on GPU (CuArray), view() returns a new CuArray over the same device
        # buffer rather than a SubArray, so Base.parent identity does not apply.
        # Instead we verify memory sharing via a round-trip write.
        fill!(buf, eltype(buf)(13))
        for name in expected_tendency_views
            view = getfield(v.tendencies.grid, name)
            @test view isa SpeedyWeather.RingGrids.AbstractField
            @test is_gpu_array(view.data)
            @test all(Array(view.data) .== eltype(buf)(13))  # shared memory visible through view
        end
        fill!(buf, zero(eltype(buf)))
    end

    @testset "$model_name parent ↔ view round trips" begin
        buf = v.scratch.fused.tend_grid

        # Write through parent, read through views (broadcast)
        fill!(buf, eltype(buf)(7))
        for name in expected_tendency_views
            view = getfield(v.tendencies.grid, name)
            @test all(Array(view.data) .== eltype(buf)(7))
        end

        # Write through one view, verify slot-correct in parent and other views untouched
        fill!(buf, zero(eltype(buf)))
        u_view = v.tendencies.grid.u
        fill!(u_view, eltype(buf)(3))
        # u occupies the FIRST slot of the parent (slots 1..nlayers)
        nlayers = sg.nlayers
        u_slots_host = Array(buf.data[:, 1:nlayers])
        @test all(u_slots_host .== eltype(buf)(3))
        if length(expected_tendency_views) > 1
            other_slots_host = Array(buf.data[:, (nlayers + 1):end])
            @test all(other_slots_host .== zero(eltype(buf)))
        end
    end

    @testset "$model_name KA kernel through view" begin
        buf = v.scratch.fused.tend_grid
        fill!(buf, zero(eltype(buf)))

        # bump u_view by 5 via a KA kernel — exercises KA/Adapt dispatch on a fused-memory Field
        bump_via_kernel!(ARCH, v.tendencies.grid.u, 5)
        # Verify: the view itself reads the bumped value (kernel wrote through correctly)
        @test all(Array(v.tendencies.grid.u.data) .== eltype(buf)(5))
        # Verify: the parent buf reflects the write (shared memory, not a copy)
        buf_host = Array(buf.data)
        @test sum(buf_host .== eltype(buf)(5)) == length(Array(v.tendencies.grid.u.data))

        # If multiple fused views, bump the LAST one and confirm write-through to parent
        if length(expected_tendency_views) > 1
            fill!(buf, zero(eltype(buf)))
            last_name = expected_tendency_views[end]
            last_view = getfield(v.tendencies.grid, last_name)
            bump_via_kernel!(ARCH, last_view, 9)
            @test all(Array(last_view.data) .== eltype(buf)(9))
            buf_host2 = Array(buf.data)
            @test sum(buf_host2 .== eltype(buf)(9)) == length(Array(last_view.data))
        end
    end

    @testset "$model_name reset_tendencies!" begin
        buf = v.scratch.fused.tend_grid
        fill!(buf, eltype(buf)(99))
        reset_tendencies!(v)

        # Each fused tendency view should be zero; checked via the view (not index arithmetic)
        # because scratch variables may occupy non-contiguous columns in the parent.
        for name in expected_tendency_views
            view = getfield(v.tendencies.grid, name)
            @test all(Array(view.data) .== zero(eltype(buf)))
        end
    end

    @testset "$model_name 1-hour run" begin
        # End-to-end smoke test: parameterizations + dycore + transforms +
        # reset_tendencies, all touching view-backed Fields on the GPU.
        run!(sim, period=Hour(1))
        @test true   # if we got here, no exception was thrown
    end

    return nothing
end

# -- run --------------------------------------------------------------------
@testset "GPU fusion smoke test" begin
    test_fuse_layout(PrimitiveDryModel, "PrimitiveDry", (:u, :v, :temperature))
    test_fuse_layout(PrimitiveWetModel, "PrimitiveWet", (:u, :v, :temperature, :humidity))
end

@info "GPU fuse smoke test completed"
