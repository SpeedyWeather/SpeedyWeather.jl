module SpeedyWeatherCUDAExt

using SpeedyWeather
import CUDA: CUDA, CUDAKernels, CuArray, CUFFT
import AbstractFFTs
using DocStringExtensions

# for RingGrids and LowerTriangularMatrices:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:CuArray}) = CuArray
LowerTriangularMatrices.nonparametric_type(::Type{<:CuArray}) = CuArray

SpeedyWeather.default_array_type(::Type{GPU}) = CuArray

# Override FFT package deciding function
which_FFT_package(::Type{<:CuArray{<:AbstractFloat}}) = CUFFT 

# Wrap function to convert zeros to a CuArray
cu_zeros(args...; kwargs...) = CUDA.cu(zeros(args...; kwargs...)) 
SpeedyWeather.SpeedyTransforms.get_zeros_func(::Type{<:CuArray}) = cu_zeros


"""$(TYPEDSIGNATURES)
Util function to generate FFT plans based on the array type of the fake Grid 
data provided. Directly indexes the fake data grid, which is more allocate-y but 
necessary when using CuArrays due to current limitations with the CUFFT planning 
API."""
function SpeedyWeather.SpeedyTransforms.plan_FFTs!(
    rfft_plans::Vector{AbstractFFTs.Plan},
    brfft_plans::Vector{AbstractFFTs.Plan},
    rfft_plans_1D::Vector{AbstractFFTs.Plan},
    brfft_plans_1D::Vector{AbstractFFTs.Plan},
    fake_grid_data::AbstractGridArray{NF, N, <:CuArray{NF}},
    scratch_memory_north::CuArray{Complex{NF}},
    rings::AbstractArray,
    nlons::Vector{<:Int},
) where {NF<:AbstractFloat, N}
    # Determine the FFT package, which in this particular case is CUFFT
    FFT_package = CUFFT

    # For each ring generate an FFT plan (for all layers and for a single layer)
    for (j, nlon) in enumerate(nlons)
        @inbounds real_matrix_input = fake_grid_data.data[rings[j], :]
        @inbounds complex_matrix_input = scratch_memory_north[1:nlon÷2 + 1, :, j]
        @inbounds real_vector_input = fake_grid_data.data[rings[j], 1]
        @inbounds complex_vector_input = scratch_memory_north[1:nlon÷2 + 1, 1, j]

        rfft_plans[j] = FFT_package.plan_rfft(real_matrix_input, 1)
        brfft_plans[j] = FFT_package.plan_brfft(complex_matrix_input, nlon, 1)
        rfft_plans_1D[j] = FFT_package.plan_rfft(real_vector_input, 1)
        brfft_plans_1D[j] = FFT_package.plan_brfft(complex_vector_input, nlon, 1)
    end

    return rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D
end


# DEVICE SETUP FOR CUDA

"""$(TYPEDSIGNATURES)
Return default used device for internal purposes, either `CPU` or `GPU` if a GPU is available."""
Device() = CUDA.functional() ? GPU() : CPU()
SpeedyWeather.DeviceSetup() = DeviceSetup(Device(), Device_KernelAbstractions(Device()), workgroup_size(Device()))

"""$(TYPEDSIGNATURES)
Return default used device for KernelAbstractions, either `CPU` or `CUDADevice` if a GPU is available."""
SpeedyWeather.Device_KernelAbstractions() = CUDA.functional() ? KernelAbstractions.CUDADevice : KernelAbstractions.CPU
SpeedyWeather.Device_KernelAbstractions(::GPU) = KernelAbstractions.CUDADevice

SpeedyWeather.DeviceArray(::GPU, x) = Adapt.adapt(CuArray, x)

"""$(TYPEDSIGNATURES)
Returns a `CuArray` when `device<:GPU` is used. Doesn't uses `adapt`, therefore always returns CuArray."""
SpeedyWeather.DeviceArrayNotAdapt(::GPU, x) = CuArray(x)

end # module