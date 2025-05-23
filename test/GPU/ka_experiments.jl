# these tests are not used in any actual test scripts in the CI 
# they are just for testing and experimenting with the KA kernels


using KernelAbstractions
import SpeedyWeather: on_architecture, CPU, CPUStatic

# To-Do write tests for each type of dims_type in the kernel launching util, 
# the tests currently below will be removed when the KA becomes the only one
using SpeedyWeather, BenchmarkTools, Test, CUDA
  
if CUDA.functional()
    arch = SpeedyWeather.CUDAGPU()
else 
    arch = SpeedyWeather.CPU()
end 

arch_static = SpeedyWeather.CPUStatic()

L = 129
M = 128
N = 8
    # ∇²!

L = 513
M = 512
N = 16

NF= Float32
alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))
#alms = rand(LowerTriangularArray{Complex{NF}},33, 32)
alms_cpu = on_architecture(CPU(), alms)

alms2 = deepcopy(alms_cpu)
alms3 = deepcopy(alms)

S = SpectralTransform(alms)
S_cpu = SpectralTransform(alms_cpu)
S_static = SpectralTransform(alms, architecture=arch_static)

# so far: KA 5x-10x faster on CPU
SpeedyWeather.SpeedyTransforms.∇²!(alms2, alms_cpu, S_cpu);
SpeedyWeather.SpeedyTransforms.∇²_KA!(alms3, alms, S);

@test on_architecture(CPU(), alms3) ≈ alms2

# Divergence
alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))
alms2 = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))

alms_cpu = on_architecture(CPU(), alms)
alms2_cpu = on_architecture(CPU(), alms2)

alms3 = copy(alms_cpu)
alms4 = copy(alms)

# so far KA 4x slower on CPU
SpeedyWeather.SpeedyTransforms.divergence!(alms3, alms_cpu, alms2_cpu, S_cpu)
SpeedyWeather.SpeedyTransforms.divergence_KA!(alms4, alms, alms2, S)
SpeedyWeather.SpeedyTransforms.divergence_new!(alms4, alms, alms2, S)
SpeedyWeather.SpeedyTransforms.divergence_KA_veryold_typed!(alms4, alms, alms2, S)

@test on_architecture(CPU(), alms4) ≈ alms3

# ∇! 
alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))
alms_cpu = on_architecture(CPU(), alms)
alms1 = copy(alms_cpu)
alms2 = copy(alms_cpu)
alms3 = copy(alms)
alms4 = copy(alms)

# so far KA 3x slower on CPU
SpeedyWeather.SpeedyTransforms.∇!(alms1, alms2, alms_cpu, S_cpu)

SpeedyWeather.SpeedyTransforms.∇_new!(alms1, alms2, alms_cpu, S_cpu)

SpeedyWeather.SpeedyTransforms.∇_3KA!(alms3, alms4, alms, S)
SpeedyWeather.SpeedyTransforms.∇_KA!(alms3, alms4, alms, S)

@test alms1 ≈ on_architecture(CPU(), alms3)
@test alms2 ≈ on_architecture(CPU(), alms4)
 


# UV_from_vordiv    

alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))
alms2 = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}},L, M, N))

div = divergence(alms, alms2, S)
vor = curl(alms, alms2, S)

div2 = copy(div)
vor2 = copy(vor)

U = similar(div)
V = similar(div)

U2 = similar(div)
V2 = similar(div)

U3 = similar(div)
V3 = similar(div)

# so far KA 4x slower on CPU
SpeedyWeather.SpeedyTransforms.UV_from_vordiv!(U, V, vor, div, S)
SpeedyWeather.SpeedyTransforms.UV_from_vordiv_KA!(U2, V2, vor, div, S)
SpeedyWeather.SpeedyTransforms.UV_from_vordiv_KA_split!(U3, V3, vor, div, S)

@test U ≈ U2
@test V ≈ V2
@test U ≈ U3
@test V ≈ V3
