import .SpeedyWeather: get_step, get_steps, nsteps, get_nsteps, DEFAULT_NSTEPS,
    which_step, which_prognostic_step, which_tendency_step,
    get_prognostic_step, get_tendency_step,
    prognostic_steps, prognostic_grid_steps, prognostic_spectral_steps,
    tendency_steps, tendency_grid_steps, tendency_spectral_steps,
    DummyTimeStepper, ArrayDimensions

@testset "get_step: fields" begin
    spectral_grid = SpectralGrid(truncation = 11, nlayers = 3)
    grid = spectral_grid.grid
    npoints = spectral_grid.npoints

    # 2D field (horizontal only), no step dimension, step must be 1
    field_2d = rand(spectral_grid.GridVariable2D, grid)
    @test get_step(field_2d, 1) == field_2d
    @test size(get_step(field_2d, 1)) == (npoints,)
    # a variable without a time dimension has a single step: the full variable
    @test get_step(field_2d) == field_2d
    @test nsteps(field_2d) == 1

    # horizontal + time, the 2nd dimension IS the step dimension
    field_xyt = rand(spectral_grid.GridVariableXYT, grid, 2)
    @test ArrayDimensions.hastime(field_xyt)
    @test nsteps(field_xyt) == 2
    for step in 1:2
        @test size(get_step(field_xyt, step)) == (npoints,)
        @test get_step(field_xyt, step) == field_xyt[:, step]
    end
    @test get_step(field_xyt) == get_step(field_xyt, 2)      # default is the last step

    # horizontal + vertical, the 2nd dimension is NOT the step dimension
    # so the full array is returned and `step` ignored
    field_xyz = rand(spectral_grid.GridVariableXYZ, grid, 3)
    @test ArrayDimensions.hastime(field_xyz) == false
    @test nsteps(field_xyz) == 1
    for step in 1:3
        @test get_step(field_xyz, step) == field_xyz
        @test size(get_step(field_xyz, step)) == (npoints, 3)
    end

    # horizontal + vertical + time, steps in the 3rd dimension
    field_xyzt = rand(spectral_grid.GridVariableXYZT, grid, 3, 2)
    for step in 1:2
        @test size(get_step(field_xyzt, step)) == (npoints, 3)
        @test get_step(field_xyzt, step) == field_xyzt[:, :, step]
    end
    @test get_step(field_xyzt) == get_step(field_xyzt, 2)

    # views alias the parent, no copies
    get_step(field_xyt, 1)[1] = 42
    @test field_xyt[1, 1] == 42
end

@testset "get_step: LowerTriangularArrays" begin
    spectral_grid = SpectralGrid(truncation = 11, nlayers = 3)
    spectrum = spectral_grid.spectrum
    n = size(rand(spectral_grid.SpectralVariable2D, spectrum), 1)

    # 2D spectral (horizontal only), no step dimension, step must be 1
    spec_2d = rand(spectral_grid.SpectralVariable2D, spectrum)
    @test get_step(spec_2d, 1) == spec_2d
    @test get_step(spec_2d) == spec_2d
    @test nsteps(spec_2d) == 1

    # horizontal + time, the 2nd dimension IS the step dimension
    spec_lmt = rand(spectral_grid.SpectralVariableXYT, spectrum, 2)
    @test ArrayDimensions.hastime(spec_lmt)
    @test nsteps(spec_lmt) == 2
    for step in 1:2
        @test size(get_step(spec_lmt, step)) == (n,)
        @test get_step(spec_lmt, step) == spec_lmt[:, step]
    end
    @test get_step(spec_lmt) == get_step(spec_lmt, 2)

    # horizontal + vertical, the 2nd dimension is NOT the step dimension
    spec_lmz = rand(spectral_grid.SpectralVariableXYZ, spectrum, 3)
    @test ArrayDimensions.hastime(spec_lmz) == false
    @test nsteps(spec_lmz) == 1
    for step in 1:3
        @test get_step(spec_lmz, step) == spec_lmz
        @test size(get_step(spec_lmz, step)) == (n, 3)
    end

    # horizontal + vertical + time, steps in the 3rd dimension
    spec_lmzt = rand(spectral_grid.SpectralVariableXYZT, spectrum, 3, 2)
    for step in 1:2
        @test size(get_step(spec_lmzt, step)) == (n, 3)
        @test get_step(spec_lmzt, step) == spec_lmzt[:, :, step]
    end
    @test get_step(spec_lmzt) == get_step(spec_lmzt, 2)

    # views alias the parent, no copies
    get_step(spec_lmt, 1)[1] = 42
    @test spec_lmt[1, 1] == 42
end

@testset "get_step: plain arrays (as unwrapped inside GPU kernels)" begin
    # Adapt unwraps Fields/LowerTriangularArrays to their bare data arrays,
    # get_step has to keep working on those with the same semantics
    vector = rand(Float32, 5)
    matrix = rand(Float32, 5, 2)
    array3 = rand(Float32, 5, 3, 2)

    @test get_step(vector, 1) == vector
    @test get_step(matrix, 2) == matrix[:, 2]
    @test get_step(array3, 2) == array3[:, :, 2]
    @test get_step(matrix) == matrix[:, 2]      # default: last index
end

@testset "get_steps" begin
    spectral_grid = SpectralGrid(truncation = 11, nlayers = 3)
    grid = spectral_grid.grid

    # variables WITHOUT a time dimension: a 1-tuple holding the full variable
    field_2d = rand(spectral_grid.GridVariable2D, grid)
    @test length(get_steps(field_2d)) == 1
    @test get_steps(field_2d)[1] == field_2d

    field_xyz = rand(spectral_grid.GridVariableXYZ, grid, 3)
    @test length(get_steps(field_xyz)) == 1      # the 3 layers are NOT steps
    @test get_steps(field_xyz)[1] == field_xyz

    spec_lmz = rand(spectral_grid.SpectralVariableXYZ, spectral_grid.spectrum, 3)
    @test length(get_steps(spec_lmz)) == 1
    @test get_steps(spec_lmz)[1] == spec_lmz

    # variables WITH a time dimension: one view per step

    field_xyt = rand(spectral_grid.GridVariableXYT, grid, 2)
    steps = get_steps(field_xyt)
    @test length(steps) == 2
    @test steps == (get_step(field_xyt, 1), get_step(field_xyt, 2))

    field_xyzt = rand(spectral_grid.GridVariableXYZT, grid, 3, 2)
    steps = get_steps(field_xyzt)
    @test length(steps) == 2
    @test steps == (get_step(field_xyzt, 1), get_step(field_xyzt, 2))

    # compile-time length via Val(N), a subset of the steps
    steps = get_steps(field_xyzt, Val(2))
    @test steps isa NTuple{2}
    @test steps == (get_step(field_xyzt, 1), get_step(field_xyzt, 2))
    @test get_steps(field_xyt, Val(1)) == (get_step(field_xyt, 1),)

    # spectral equivalent
    spectrum = spectral_grid.spectrum
    spec_lmt = rand(spectral_grid.SpectralVariableXYT, spectrum, 2)
    @test get_steps(spec_lmt) == (get_step(spec_lmt, 1), get_step(spec_lmt, 2))
    @test get_steps(spec_lmt, Val(2)) isa NTuple{2}
end

@testset "number of steps per time stepper" begin
    # fallbacks: one step for everything
    dummy = DummyTimeStepper()
    @test prognostic_steps(dummy) == 1
    @test prognostic_grid_steps(dummy) == 1
    @test prognostic_spectral_steps(dummy) == 1
    @test tendency_steps(dummy) == 1
    @test tendency_grid_steps(dummy) == 1
    @test tendency_spectral_steps(dummy) == 1
    @test DEFAULT_NSTEPS == (
        prognostic_grid = 1, prognostic_spectral = 1,
        tendency_grid = 1, tendency_spectral = 1,
    )

    spectral_grid = SpectralGrid(truncation = 11, nlayers = 3)

    # Leapfrog: 2 spectral steps for prognostic variables, 1 for tendencies
    for Model in (ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid, time_stepping = Leapfrog(spectral_grid))
        leapfrog = model.time_stepping
        @test leapfrog isa SpeedyWeather.AbstractLeapfrog

        nsteps = get_nsteps(leapfrog, model)
        @test nsteps.prognostic_spectral == 2
        @test nsteps.tendency_grid == 1
        @test nsteps.tendency_spectral == 1
        # 2D models only need one grid step, primitive equations two (parameterizations
        # are evaluated at the previous step)
        @test nsteps.prognostic_grid == (model isa PrimitiveEquation ? 2 : 1)

        # model-less fallbacks agree with the model-dispatched methods
        @test prognostic_spectral_steps(leapfrog, model) == prognostic_spectral_steps(leapfrog)
        @test tendency_steps(leapfrog, model) == tendency_steps(leapfrog)
        @test tendency_grid_steps(leapfrog) == tendency_steps(leapfrog)
        @test tendency_spectral_steps(leapfrog) == tendency_steps(leapfrog)
    end

    # NCycleLorenz: one prognostic step but 2 spectral tendency steps (F and G terms)
    model = BarotropicModel(spectral_grid, time_stepping = NCycleLorenz(spectral_grid))
    ncycle = model.time_stepping
    nsteps = get_nsteps(ncycle, model)
    @test nsteps.prognostic_grid == 1
    @test nsteps.prognostic_spectral == 1
    @test nsteps.tendency_grid == 1
    @test nsteps.tendency_spectral == 2
    @test prognostic_steps(ncycle) == 1
end

@testset "which_step dispatch" begin
    spectral_grid = SpectralGrid(truncation = 11, nlayers = 3)
    model = PrimitiveDryModel(spectral_grid)
    leapfrog = model.time_stepping
    var = rand(spectral_grid.SpectralVariableXYT, spectral_grid.spectrum, 2)

    # fallback for a component without a specialised method is step 1
    component = model.orography
    @test which_step(var, leapfrog, component) == 1
    @test which_prognostic_step(var, leapfrog, component) == 1
    @test which_tendency_step(var, leapfrog, component) == 1

    # dispatch with model falls back to dispatch without model
    @test which_step(var, leapfrog, component, model) == which_step(var, leapfrog, component)
    @test which_prognostic_step(var, leapfrog, component, model) == which_prognostic_step(var, leapfrog, component)
    @test which_tendency_step(var, leapfrog, component, model) == which_tendency_step(var, leapfrog, component)

    # Leapfrog reads the current (2nd) step for the transforms, but only 1 for parameterizations
    transform = model.spectral_transform
    @test which_prognostic_step(var, leapfrog, transform) == 2
    @test get_prognostic_step(var, leapfrog, transform) == get_step(var, 2)
    @test get_prognostic_step(var, leapfrog, transform, model) == get_step(var, 2)

    parameterization = model.convection
    @test which_prognostic_step(var, leapfrog, parameterization) == 1
    @test get_prognostic_step(var, leapfrog, parameterization) == get_step(var, 1)

    # tendencies have only one step in Leapfrog
    tendency = rand(spectral_grid.SpectralVariableXYT, spectral_grid.spectrum, 1)
    @test which_tendency_step(tendency, leapfrog, transform) == 1
    @test get_tendency_step(tendency, leapfrog, transform) == get_step(tendency, 1)

    # get_step with a component resolves via which_step
    @test get_step(var, leapfrog, component) == get_step(var, which_step(var, leapfrog, component))
    @test get_step(var, leapfrog, component, model) == get_step(var, 1)
end
