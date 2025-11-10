module SpeedyTransformsAMDGPUExt
    
    import AMDGPU: AMDGPU, rocFFT, RocArray
    import AbstractFFTs
    import LinearAlgebra
    using DocStringExtensions

    using SpeedyTransforms
    using SpeedyTransforms.RingGrids
    using SpeedyTransforms.LowerTriangularArrays
    
    # Override FFT package deciding function
    SpeedyTransforms.which_FFT_package(::Type{<:RocArray{<:AbstractFloat}}) = rocFFT 

    """$(TYPEDSIGNATURES)
    Util function to generate FFT plans based on the array type of the fake Grid 
    data provided. Uses indexing as we seemingly can't use views for FFT planning 
    with CUFFT."""
    function SpeedyTransforms.plan_FFTs!(
        rfft_plans::Vector{AbstractFFTs.Plan},
        brfft_plans::Vector{AbstractFFTs.Plan},
        rfft_plans_1D::Vector{AbstractFFTs.Plan},
        brfft_plans_1D::Vector{AbstractFFTs.Plan},
        fake_field_data::AbstractField{NF, N, <:RocArray{NF}},
        scratch_memory_north::RocArray{Complex{NF}},
        rings,
        nlons::Vector{<:Int}
    ) where {NF<:AbstractFloat, N}
        # Determine which FFT package to use
        FFT_package = SpeedyTransforms.which_FFT_package(RocArray{NF})

        # For each ring generate an FFT plan (for all layers and for a single layer)
        for (j, nlon) in enumerate(nlons)
            real_matrix_input = fake_field_data.data[rings[j], :]
            complex_matrix_input = scratch_memory_north[1:nlon÷2 + 1, :, j]
            real_vector_input = fake_field_data.data[rings[j], 1]
            complex_vector_input = scratch_memory_north[1:nlon÷2 + 1, 1, j]

            rfft_plans[j] = FFT_package.plan_rfft(real_matrix_input, 1)
            brfft_plans[j] = FFT_package.plan_brfft(complex_matrix_input, nlon, 1)
            rfft_plans_1D[j] = FFT_package.plan_rfft(real_vector_input, 1)
            brfft_plans_1D[j] = FFT_package.plan_brfft(complex_vector_input, nlon, 1)
        end

        return rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D
    end

    """$(TYPEDSIGNATURES)
    (Forward) FFT, applied in zonal direction of `field` provided. This is the 
    GPU/CUDA equivalent of the `apply_batched_fft!` function in the CPU version. 
    Uses indexing as we seemingly can't use views with the FFT planning with CUFFT.
    """
    function SpeedyTransforms._apply_batched_fft!(
        f_out::RocArray{<:Complex, 3},
        field::AbstractField{NF, N, <:RocArray},
        S::SpectralTransform, 
        j::Int,
        nfreq::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true
    ) where {NF<:AbstractFloat, N}
        rfft_plan = S.rfft_plans[j]     # FFT planned wrt nlon on ring
        nlayers = size(field, 2)        # number of vertical layers

        if not_equator
            view(f_out, 1:nfreq, 1:nlayers, j) .= rfft_plan * field.data[ilons, :]
        else
            fill!(f_out[1:nfreq, 1:nlayers, j], 0)
        end
    end

    """$(TYPEDSIGNATURES)
    (Inverse) FFT, applied in zonal direction of `field` provided. This is the
    GPU/CUDA equivalent of the `apply_batched_fft!` function in the CPU version.
    Uses indexing as we seemingly can't use views with the FFT planning with CUFFT.
    """
    function SpeedyTransforms._apply_batched_fft!(
        field::AbstractField{NF, N, <:RocArray},
        g_in::RocArray{<:Complex, 3},
        S::SpectralTransform,
        j::Int,
        nlon::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true
    ) where {NF<:AbstractFloat, N}
        brfft_plan = S.brfft_plans[j]   # FFT planned wrt nlon on ring
        nlayers = size(field, 2)        # number of vertical layers
        nfreq = nlon÷2 + 1              # linear max Fourier frequency wrt to nlon

        if not_equator
            view(field.data, ilons, :) .= brfft_plan * g_in[1:nfreq, 1:nlayers, j]
        end
    end
    
end 
