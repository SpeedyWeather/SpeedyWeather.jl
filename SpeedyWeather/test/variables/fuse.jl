using SpeedyWeather
using SpeedyWeather: AbstractVariable, GridVariable, TendencyVariable, ScratchVariable,
    DynamicsVariable, PrognosticVariable,
    Grid2D, Grid3D, Grid4D, Spectral2D, Spectral3D, Spectral4D,
    allocate_fused, build_fuse_parents, fuse_family, fused_slots

# Helpers used by all testsets — build a model and unwrap a few useful pieces.
function _testmodel(; trunc = 10, nlayers = 4)
    sg = SpectralGrid(; trunc, nlayers)
    return PrimitiveDryModel(sg)
end

@testset "fuse_family classification" begin
    @test fuse_family(Grid2D()) === :grid
    @test fuse_family(Grid3D()) === :grid
    @test fuse_family(Grid4D()) === :grid
    @test fuse_family(Spectral2D()) === :spectral
    @test fuse_family(Spectral3D()) === :spectral
    @test fuse_family(Spectral4D()) === :spectral
end

@testset "allocate_fused: pure Grid3D group" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        TendencyVariable(:u, Grid3D(), namespace = :grid, fuse = :f),
        TendencyVariable(:v, Grid3D(), namespace = :grid, fuse = :f),
    ]
    parent, views, slots = allocate_fused(vars, model)
    @test size(parent.data, 2) == 2 * 4
    @test slots == [1:4, 5:8]
    # views alias the parent at the declared slots
    for (view, slot) in zip(views, slots)
        @test Base.parent(view.data) === parent.data
        @test parentindices(view.data)[end] == slot
    end
end

@testset "allocate_fused: mixed Grid2D + Grid3D" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        GridVariable(:a3d, Grid3D(), namespace = :g, fuse = :m),
        GridVariable(:b2d, Grid2D(), namespace = :g, fuse = :m),
        GridVariable(:c3d, Grid3D(), namespace = :g, fuse = :m),
        GridVariable(:d2d, Grid2D(), namespace = :g, fuse = :m),
    ]
    parent, views, slots = allocate_fused(vars, model)
    # 3d (4) + 2d (1) + 3d (4) + 2d (1) = 10 columns
    @test size(parent.data, 2) == 10
    @test slots == [1:4, 5:5, 6:9, 10:10]
    # 3D members keep the layer dim (2D Field data), 2D members collapse it (1D Field data)
    @test ndims(views[1].data) == 2   # Grid3D view
    @test ndims(views[2].data) == 1   # Grid2D view
    @test ndims(views[3].data) == 2
    @test ndims(views[4].data) == 1
    # All views are views of the same parent buffer
    for view in views
        @test Base.parent(view.data) === parent.data
    end
end

@testset "allocate_fused: mixed Spectral2D + Spectral3D" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        DynamicsVariable(:a3d, Spectral3D(), namespace = :s, fuse = :ms),
        DynamicsVariable(:b2d, Spectral2D(), namespace = :s, fuse = :ms),
    ]
    parent, views, slots = allocate_fused(vars, model)
    @test size(parent.data, 2) == 5
    @test slots == [1:4, 5:5]
    @test ndims(views[1].data) == 2   # Spectral3D view (LowerTriangularArray{T,2})
    @test ndims(views[2].data) == 1   # Spectral2D view (LowerTriangularMatrix)
    for view in views
        @test Base.parent(view.data) === parent.data
    end
end

@testset "fuse validation: mixing grid + spectral families is rejected" begin
    model = _testmodel(nlayers = 4)
    bad = AbstractVariable[
        TendencyVariable(:u, Grid3D(), namespace = :x, fuse = :bad),
        TendencyVariable(:p, Spectral2D(), namespace = :x, fuse = :bad),
    ]
    @test_throws ErrorException build_fuse_parents(bad, model)
end

@testset "allocate_fused: pure Grid4D group (parent is 4D)" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        PrognosticVariable(:a, Grid4D(n = 2), namespace = :g, fuse = :p4),
        PrognosticVariable(:b, Grid4D(n = 2), namespace = :g, fuse = :p4),
    ]
    parent, views, slots = allocate_fused(vars, model)
    @test ndims(parent.data) == 3
    @test size(parent.data) == (size(parent.data, 1), 2 * 4, 2)   # (npoints, 2·nlayers, n=2)
    @test slots == [1:4, 5:8]
    for (view, slot) in zip(views, slots)
        @test Base.parent(view.data) === parent.data
        @test ndims(view.data) == 3                                # (npoints, nlayers, n)
        @test size(view.data) == (size(parent.data, 1), 4, 2)
        @test parentindices(view.data)[2] == slot
    end
end

@testset "allocate_fused: mixed Grid3D + Grid4D (parent 4D, Grid3D occupies 1 slot)" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        PrognosticVariable(:a4, Grid4D(n = 3), namespace = :g, fuse = :m4),
        TendencyVariable(:b3, Grid3D(),       namespace = :g, fuse = :m4),
        PrognosticVariable(:c4, Grid4D(n = 3), namespace = :g, fuse = :m4),
    ]
    parent, views, slots = allocate_fused(vars, model)
    @test ndims(parent.data) == 3
    # 4 (Grid4D layer slots) + 1 (Grid3D collapsed) + 4 (Grid4D) = 9 along layer axis;
    # trailing dim n = 3 shared.
    @test size(parent.data) == (size(parent.data, 1), 9, 3)
    @test slots == [1:4, 5:5, 6:9]

    # Grid4D members keep the layer dim → 3D view (npoints, nlayers, n)
    @test ndims(views[1].data) == 3
    @test size(views[1].data) == (size(parent.data, 1), 4, 3)
    @test parentindices(views[1].data)[2] == 1:4

    # Grid3D-in-4D member collapses the layer dim → 2D view (npoints, n)
    @test ndims(views[2].data) == 2
    @test size(views[2].data) == (size(parent.data, 1), 3)
    @test parentindices(views[2].data)[2] == 5   # scalar layer index

    @test ndims(views[3].data) == 3
    @test parentindices(views[3].data)[2] == 6:9
end

@testset "allocate_fused: pure Spectral4D + Spectral3D mix" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        PrognosticVariable(:a4, Spectral4D(n = 2), namespace = :s, fuse = :ms4),
        DynamicsVariable(:b3,   Spectral3D(),      namespace = :s, fuse = :ms4),
    ]
    parent, views, slots = allocate_fused(vars, model)
    @test ndims(parent.data) == 3                       # LowerTriangularArray with 3D data
    # layer axis: 4 (Spectral4D) + 1 (Spectral3D collapsed) = 5; trailing n = 2
    @test size(parent.data, 2) == 5
    @test size(parent.data, 3) == 2
    @test slots == [1:4, 5:5]
    @test ndims(views[1].data) == 3                     # Spectral4D view
    @test ndims(views[2].data) == 2                     # Spectral3D-in-4D collapses layer dim
    for view in views
        @test Base.parent(view.data) === parent.data
    end
end

@testset "fuse validation: mixed n in 4D group is rejected" begin
    model = _testmodel(nlayers = 4)
    bad = AbstractVariable[
        PrognosticVariable(:a, Grid4D(n = 2), namespace = :x, fuse = :badn),
        PrognosticVariable(:b, Grid4D(n = 3), namespace = :x, fuse = :badn),
    ]
    @test_throws ErrorException build_fuse_parents(bad, model)
end

@testset "fuse validation: 2D member in 4D group is rejected" begin
    model = _testmodel(nlayers = 4)
    bad = AbstractVariable[
        PrognosticVariable(:a4, Grid4D(n = 2), namespace = :x, fuse = :bad24),
        GridVariable(:b2,       Grid2D(),      namespace = :x, fuse = :bad24),
    ]
    @test_throws ErrorException build_fuse_parents(bad, model)
end

@testset "view ↔ parent aliasing through layer-axis writes (4D parent)" begin
    model = _testmodel(nlayers = 4)
    vars = AbstractVariable[
        PrognosticVariable(:a4, Grid4D(n = 2), namespace = :g, fuse = :alias4),
        TendencyVariable(:b3,   Grid3D(),      namespace = :g, fuse = :alias4),
    ]
    parent, views, _ = allocate_fused(vars, model)
    # Writing through the 3D-in-4D view should populate exactly slot 5 of the parent
    # across all `n` trailing slices (and nothing else).
    fill!(views[2].data, 1)
    @test all(parent.data[:, 5, :] .== 1)
    @test all(parent.data[:, 1:4, :] .== 0)
    # Per-step slice of the parent is contiguous in memory
    @test parent.data[:, :, 1] isa AbstractMatrix
end

@testset "fuse validation: cross-type same name is rejected" begin
    model = _testmodel(nlayers = 4)
    # Same (namespace, name) but different variable types (TendencyVariable vs ScratchVariable)
    bad = AbstractVariable[
        TendencyVariable(:u, Grid3D(), namespace = :z, fuse = :ct),
        ScratchVariable(:u, Grid3D(), namespace = :z, fuse = :ct),
    ]
    @test_throws ErrorException build_fuse_parents(bad, model)
end

@testset "fuse validation: duplicate name within group is rejected" begin
    # If two variables of the same type share (namespace, name), the (namespace, name)
    # dedup in build_fuse_parents silently treats them as the same variable. To trigger
    # the duplicate-name check we need to construct vars where dedup wouldn't fire —
    # which can't happen via build_fuse_parents. So we exercise the check via the inner
    # uniqueness assertion by calling build_fuse_parents with a manually pre-grouped Vector
    # — done implicitly through normal model construction in the broader tests.
    # Here we just verify the cross-type duplicate-name error message contains the variable name.
    model = _testmodel(nlayers = 4)
    bad = AbstractVariable[
        TendencyVariable(:u, Grid3D(), namespace = :w, fuse = :dn),
        ScratchVariable(:u, Grid3D(), namespace = :w, fuse = :dn),
    ]
    err = try
        build_fuse_parents(bad, model)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("u", err.msg)
end

# Activate these test once we actually use fused variables in the model
#=
@testset "view ↔ parent aliasing round trips" begin
    model = _testmodel(nlayers = 4)
    sim = initialize!(model)
    v = sim.variables
    buf = v.fused.tend_grid
    nlayers = model.spectral_grid.nlayers

    # Write through parent, read through views
    fill!(buf, eltype(buf)(7))
    @test all(v.tendencies.grid.u.data .== 7)
    @test all(v.tendencies.grid.v.data .== 7)
    @test all(v.tendencies.grid.temperature.data .== 7)

    # Write through one view, verify the right slot is touched in parent and other views untouched
    fill!(buf, zero(eltype(buf)))
    fill!(v.tendencies.grid.u, eltype(buf)(3))
    @test all(buf.data[:, 1:nlayers] .== 3)
    @test all(buf.data[:, (nlayers + 1):end] .== 0)
    @test all(v.tendencies.grid.v.data .== 0)
end

@testset "reset_tendencies! zeros fused tendency slots through views" begin
    model = _testmodel(nlayers = 4)
    sim = initialize!(model)
    v = sim.variables
    buf = v.fused.tend_grid
    nlayers = model.spectral_grid.nlayers

    fill!(buf, eltype(buf)(99))
    SpeedyWeather.reset_tendencies!(v)

    # Tendency views must be zero
    @test all(v.tendencies.grid.u.data .== 0)
    @test all(v.tendencies.grid.v.data .== 0)
    @test all(v.tendencies.grid.temperature.data .== 0)

    # The first 3·nlayers slots of the parent (the tendency slots: u, v, T) are zero
    @test all(buf.data[:, 1:(3 * nlayers)] .== 0)
    # Non-tendency fused slots (scratch a, b) were NOT touched by reset_tendencies!
    @test all(buf.data[:, (3 * nlayers + 1):end] .== 99)
end

@testset "fused parent lives at vars.fused.<sym>; standalone variables unaffected" begin
    model = _testmodel(nlayers = 4)
    sim = initialize!(model)
    v = sim.variables

    @test haskey(v.fused, :tend_grid)
    # Parent is not duplicated under other groups
    @test !haskey(v.tendencies.grid, :tend_grid)
    @test !haskey(v.scratch.grid, :tend_grid)
    @test !haskey(v.scratch, :fused)

    # Standalone (non-fused) tendency variables still work normally
    @test haskey(v.tendencies, :vorticity)   # spectral tendency, not fused
    @test !(v.tendencies.vorticity.data isa SubArray)
end
=#