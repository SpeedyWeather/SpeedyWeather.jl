@testset "Horizontal diffusion of random" begin
    @testset for HD in (HyperDiffusion, SpectralFilter)
        @testset for NF in (Float32, Float64)

            spectral_grid = SpectralGrid(; NF)
            horizontal_diffusion = HD(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; horizontal_diffusion)
            simulation = initialize!(model)
            initialize!(simulation)

            (; variables) = simulation
            (; expl, impl) = model.horizontal_diffusion
            vor = get_step(variables.prognostic.vorticity, 1)
            vor_tend = get_step(variables.tendencies.vorticity, 1)

            # apply diffusion
            SpeedyWeather.horizontal_diffusion!(vor_tend, vor, expl, impl)

            # diffusion tendency has opposite sign (real/imag respectively)
            # than prognostic variable to act as a dissipation
            (; spectrum) = model.spectral_transform
            for k in eachmatrix(vor, vor_tend)
                for m in 1:spectrum.mmax
                    for l in max(2, m):(spectrum.lmax - 1)
                        @test -sign(real(vor[l, m, k])) == sign(real(vor_tend[l, m, k]))
                        @test -sign(imag(vor[l, m, k])) == sign(imag(vor_tend[l, m, k]))
                    end
                end
            end

            vor2 = get_step(variables.prognostic.vorticity, 2)
            SpeedyWeather.update_prognostic!(variables, model)

            @test any(vor .!= vor2)    # check that at least some coefficients are different
            @test any(vor .== vor2)    # check that at least some coefficients are identical

            # damping should not increase real or imaginary part of variable
            @test all(abs.(real.(vor2)) .<= abs.(real.(vor)))
            @test all(abs.(imag.(vor2)) .<= abs.(imag.(vor)))
        end
    end
end

@testset "Implicit diffusion uses the prognostic time step" begin
    # The implicit part has to be 1/(1 - Δt*expl) with Δt the step the prognostic variables
    # are actually advanced with (2Δt for leapfrog, = default_time_step), so that
    # (tend + expl*var)*impl followed by the leapfrog update reproduces the exact implicit
    # solution var_new = (var + Δt*tend)/(1 - Δt*expl). Any other factor over- or
    # under-damps the full tendency, not just the diffusion.
    @testset for HD in (HyperDiffusion, SpectralFilter)
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 4)
        horizontal_diffusion = HD(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; horizontal_diffusion)
        simulation = initialize!(model)

        (; expl, impl, expl_div, impl_div) = model.horizontal_diffusion

        # arrays are precomputed with the radius-scaled time step, see initialize!
        Δt = SpeedyWeather.default_time_step(model.time_stepping) / model.planet.radius

        # last degree (trunc+2) is zeroed for both, exclude it here
        @test impl[1:(end - 1), :] ≈ 1 ./ (1 .- Δt .* expl[1:(end - 1), :])
        @test impl_div[1:(end - 1), :] ≈ 1 ./ (1 .- Δt .* expl_div[1:(end - 1), :])

        # a variable with zero tendency must decay by exactly 1/(1 - Δt*expl) per step
        var = rand(ComplexF32, spectral_grid.spectrum, spectral_grid.nlayers)
        tend = zeros(ComplexF32, spectral_grid.spectrum, spectral_grid.nlayers)
        SpeedyWeather.horizontal_diffusion!(tend, var, expl, impl)
        l_indices = Array(var.spectrum.l_indices)
        decay = [1 / (1 - Δt * Array(expl)[l_indices[lm], k]) for lm in eachindex(l_indices), k in 1:spectral_grid.nlayers]
        @test parent(var .+ Δt .* tend) ≈ parent(var) .* decay
    end
end
