# `copy!` and `deepcopy` of the variables on the GPU. A contiguous view of a GPU array (e.g. a member `(:, k)`
# or `(:, slots)` of a fused parent without step dimension) is not a `SubArray` but another GPU array sharing
# the parent's memory: the default `deepcopy` copies it into independent memory (the copy's member would no
# longer be a view of the copy's fused parent, the batched transforms then read stale data), and `copy!`
# between CPU and GPU would `copyto!` between a `SubArray` and a GPU array (scalar indexing).

# The arrays of all variables (except the scratch, which `copy!` skips) on the CPU as "group.name" => Array
function variables_on_host(vars)
    arrays = Dict{String, Array}()
    for group in SpeedyWeather.ALL_VARIABLE_GROUPS
        group in (:scratch, :fused) || _variables_on_host!(arrays, getfield(vars, group), string(group))
    end
    return arrays
end

function _variables_on_host!(arrays, variables::NamedTuple, prefix)
    for (name, variable) in pairs(variables)
        path = string(prefix, ".", name)
        if variable isa NamedTuple
            _variables_on_host!(arrays, variable, path)
        elseif variable isa Union{SpeedyWeather.RingGrids.AbstractField, LowerTriangularArray}
            arrays[path] = Array(on_architecture(SpeedyWeather.CPU(), variable).data)
        end
    end
    return arrays
end

# UTILITIES to check copying and deepcopy of Variables, even with array types where SubArray don't exists (CuArray, JLArray, ..)

# Test write: Write `offset + i` into fused parent `i`: a variable that then holds this value everywhere is a view of parent `i`
function fill_fused_parents!(vars, offset)
    for (i, fuse_name) in enumerate(keys(vars.fused))
        data = parent(vars.fused[fuse_name]).data
        fill!(data, eltype(data)(offset + i))
    end
    return vars
end

# Test read: The variables holding `offset + i` everywhere, i.e. the views of fused parent `i`, as "group.name" => fuse symbol
function fused_members(vars, offset)
    members = Dict{String, Symbol}()
    fuse_names = collect(keys(vars.fused))
    for group in SpeedyWeather.ALL_VARIABLE_GROUPS
        group === :fused || _fused_members!(members, getfield(vars, group), string(group), fuse_names, offset)
    end
    return members
end

function _fused_members!(members, variables::NamedTuple, prefix, fuse_names, offset)
    for (name, variable) in pairs(variables)
        path = string(prefix, ".", name)
        if variable isa NamedTuple
            _fused_members!(members, variable, path, fuse_names, offset)
        elseif variable isa Union{SpeedyWeather.RingGrids.AbstractField, LowerTriangularArray}
            data = on_architecture(SpeedyWeather.CPU(), variable).data
            isempty(data) && continue
            for (i, fuse_name) in enumerate(fuse_names)
                all(==(eltype(data)(offset + i)), data) && (members[path] = fuse_name)
            end
        end
    end
    return members
end

@testset "GPU: deepcopy keeps every fused variable a view of its copied parent" for Model in (BarotropicModel, PrimitiveWetModel)
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.GPU(), nlayers = Model <: Barotropic ? 1 : 8)
    vars = Variables(Model(spectral_grid))
    # mark the original's parents with 1000 + i before the deepcopy, so the copy starts with these values
    members = fused_members(fill_fused_parents!(vars, 1000), 1000)
    @test !isempty(members)

    # mark the copy's parents with 2000 + i: a view of a copied parent then holds 2000 + i, a detached
    # copy of a fused variable still the 1000 + i it was copied with
    copied = deepcopy(vars)
    fill_fused_parents!(copied, 2000)
    @test fused_members(copied, 2000) == members    # every fused variable of the copy is a view of its parent
    @test fused_members(vars, 1000) == members      # the original still holds 1000 + i: no shared memory
end

@testset "GPU: a deepcopy of the variables steps like the original" begin
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.GPU())
    model = PrimitiveWetModel(spectral_grid; feedback = nothing)
    simulation = initialize!(model)
    initialize!(simulation, steps = 10)
    copied = Simulation(deepcopy(simulation.variables), model)

    foreach(_ -> SpeedyWeather.time_step!(simulation), 1:10)
    foreach(_ -> SpeedyWeather.time_step!(copied), 1:10)

    for name in (:u, :v, :temperature, :humidity)
        original_field = on_architecture(SpeedyWeather.CPU(), getproperty(simulation.variables.grid, name))
        copied_field = on_architecture(SpeedyWeather.CPU(), getproperty(copied.variables.grid, name))
        @test copied_field ≈ original_field
    end
end

@testset "GPU: copy! between CPU and GPU variables" begin
    model = PrimitiveWetModel(SpectralGrid(); feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 3)                     # non-trivial data in all groups
    vars = simulation.variables
    reference = variables_on_host(vars)

    gpu_model = PrimitiveWetModel(SpectralGrid(architecture = SpeedyWeather.GPU()); feedback = nothing)
    vars_gpu = copy!(Variables(gpu_model), vars)    # CPU → GPU
    @test variables_on_host(vars_gpu) == reference
    vars_gpu2 = copy!(Variables(gpu_model), vars_gpu)   # GPU → GPU
    @test variables_on_host(vars_gpu2) == reference
    vars_back = copy!(Variables(model), vars_gpu2)  # GPU → CPU
    @test variables_on_host(vars_back) == reference
end
