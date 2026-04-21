using JLArrays
using SpeedyWeatherInternals.Architectures
import SpeedyWeatherInternals.Architectures: AbstractCPU, CPU, GPU

const JLGPU = GPU{JLArrays.JLBackend}

@testset "CPU architecture" begin
    arch = CPU()
    arch_static = CPUStatic()

    @test arch isa AbstractCPU
    @test arch isa AbstractArchitecture
    @test arch_static isa AbstractCPU

    @test array_type(arch) == Array
    @test array_type(CPU) == Array
    @test array_type(arch, Float32, 2) == Array{Float32, 2}

    @test compatible_array_types(arch) == (Array,)
    @test compatible_array_types(CPU) == (Array,)

    @test repr(arch) == "CPU"
end

@testset "JLGPU architecture" begin
    arch = GPU(JLArrays.JLBackend())

    @test arch isa GPU
    @test arch isa AbstractArchitecture

    @test array_type(arch) == JLArray
    @test array_type(arch, Float32, 2) == JLArray{Float32, 2}

    @test JLArray in compatible_array_types(arch)
    @test JLArrays.JLDeviceArray in compatible_array_types(arch)
end

@testset "architecture() inference" begin
    # scalars / nothing
    @test architecture(1.0) === nothing
    @test architecture(1f0) === nothing
    @test architecture() === nothing

    # CPU arrays
    a = rand(Float32, 4)
    @test architecture(a) isa CPU
    @test architecture(Array) isa CPU
    @test architecture(Array{Float32, 2}) isa CPU

    # SubArray of a CPU array: must be CPU
    v = view(a, :)
    @test architecture(v) isa CPU

    # JLArray
    ja = JLArray(a)
    @test architecture(ja) isa JLGPU
    @test architecture(JLArray) isa JLGPU

    # Architecture of an architecture
    cpu = CPU()
    @test architecture(cpu) === cpu
end

@testset "ismatching" begin
    cpu = CPU()
    jlgpu = GPU(JLArrays.JLBackend())

    a_cpu = rand(Float32, 4)
    a_jl = JLArray(a_cpu)

    # plain arrays
    @test ismatching(cpu, a_cpu)
    @test ismatching(cpu, Array{Float32, 1})
    @test ismatching(cpu, Array)
    @test !ismatching(cpu, a_jl)
    @test !ismatching(cpu, JLArray{Float32, 1})

    # SubArray of CPU array must match CPU
    @test ismatching(cpu, view(a_cpu, :))
    @test ismatching(cpu, view(a_cpu, 1:2))

    # JLArray matches JLGPU
    @test ismatching(jlgpu, a_jl)
    @test ismatching(jlgpu, JLArray{Float32, 1})
    @test !ismatching(jlgpu, a_cpu)

    # type-based dispatch
    @test ismatching(CPU, Array)
    @test !ismatching(CPU, JLArray)

    # arch vs arch
    @test ismatching(cpu, CPU())
    @test ismatching(jlgpu, GPU(JLArrays.JLBackend()))
    @test !ismatching(cpu, jlgpu)

    # with a view / subarray 
    @test ismatching(cpu, view(a_cpu, 1:2))
end

@testset "nonparametric_type" begin
    @test nonparametric_type(Array{Float32, 1}) == Array
    @test nonparametric_type(Array{Float64, 3}) == Array
    @test nonparametric_type(JLArray{Float32, 2}) == JLArray

    # SubArray strips down to the parent nonparametric type
    a = rand(Float32, 4)
    v = view(a, :)
    @test nonparametric_type(typeof(v)) == Array

    ja = JLArray(a)
    jv = view(ja, :)
    @test nonparametric_type(typeof(jv)) == JLArray
end

@testset "on_architecture CPU" begin
    cpu = CPU()
    jlgpu = GPU(JLArrays.JLBackend())

    a = rand(Float32, 4, 2)
    ja = JLArray(a)

    # CPU → CPU: identity
    @test on_architecture(cpu, a) === a
    @test on_architecture(cpu, a) isa Array

    # JLGPU → CPU: materialises to Array
    @test on_architecture(cpu, ja) isa Array
    @test Array(ja) ≈ on_architecture(cpu, ja)

    # CPU → JLGPU
    ja2 = on_architecture(jlgpu, a)
    @test ja2 isa JLArray
    @test Array(ja2) ≈ a

    # JLGPU → JLGPU: identity
    @test on_architecture(jlgpu, ja) === ja

    # SubArray of CPU array → JLGPU
    v = view(a, :, 1)
    jv = on_architecture(jlgpu, v)
    @test jv isa JLArray
    @test Array(jv) ≈ v

    # SubArray of JLArray → CPU
    jvsub = view(ja, :, 1)
    av = on_architecture(cpu, jvsub)
    @test av isa Array
    @test av ≈ Array(view(ja, :, 1))

    # pass-throughs
    @test on_architecture(cpu, 3.14) === 3.14
    sr = 1.0:0.1:2.0
    @test on_architecture(cpu, sr) === sr
    @test on_architecture(jlgpu, sr) === sr
end

@testset "on_architecture Tuple / NamedTuple" begin
    cpu = CPU()
    jlgpu = GPU(JLArrays.JLBackend())

    a1 = rand(Float32, 3)
    a2 = rand(Float64, 2)
    t = (a1, a2)

    t_gpu = on_architecture(jlgpu, t)
    @test t_gpu[1] isa JLArray
    @test t_gpu[2] isa JLArray

    t_back = on_architecture(cpu, t_gpu)
    @test t_back[1] isa Array
    @test t_back[2] isa Array

    nt = (x = a1, y = a2)
    nt_gpu = on_architecture(jlgpu, nt)
    @test nt_gpu.x isa JLArray
    @test keys(nt_gpu) == (:x, :y)
end

