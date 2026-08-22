import SpeedyWeather: on_architecture, GPU

@testset "Field broadcasting against bare GPU arrays" begin

    # regression test for Field .* Adjoint(CuVector)/CuVector broadcasting on GPU (#TODO PR number):
    # combining a Field's FieldGPUStyle with a bare GPU array's AbstractGPUArrayStyle used to
    # fall back to Broadcast.Unknown(), silently triggering scalar indexing on the GPU array
    arch = SpeedyWeather.GPU()
    grid = FullGaussianGrid(6, arch)
    nlayers = 4

    field = on_architecture(arch, rand(Float32, grid, nlayers))
    tendency = on_architecture(arch, zeros(Float32, grid, nlayers))
    coefs = on_architecture(arch, rand(Float32, nlayers))

    @testset "Field .* Adjoint(vector)" begin
        @test_nowarn @. tendency -= coefs' * field
    end

    @testset "Field .* vector (reshaped)" begin
        coefs_r = reshape(coefs, 1, :)
        @test_nowarn @. tendency -= coefs_r * field
    end

    @testset "Field (3D) .* Field (2D, different N, same grid)" begin
        field_2D = on_architecture(arch, rand(Float32, grid))
        @test_nowarn @. tendency -= field_2D * field
    end
end
