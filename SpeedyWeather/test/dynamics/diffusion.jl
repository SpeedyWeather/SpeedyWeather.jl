@testset "Horizontal diffusion of random" begin
    @testset for HD in (HyperDiffusion, SpectralFilter)
        @testset for NF in (Float32, Float64)

            spectral_grid = SpectralGrid(; NF)
            horizontal_diffusion = HD(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; horizontal_diffusion)
            simulation = initialize!(model)

            (; variables) = simulation
            (; expl, impl) = model.horizontal_diffusion
            vor = get_step(variables.prognostic.vor, 1)
            vor_tend = variables.tendencies.vor

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

            vor2 = get_step(variables.prognostic.vor, 2)
            SpeedyWeather.leapfrog!(vor, vor2, vor_tend, model.time_stepping.Δt, 1, model.time_stepping)

            @test any(vor .!= vor2)    # check that at least some coefficients are different
            @test any(vor .== vor2)    # check that at least some coefficients are identical

            # damping should not increase real or imaginary part of variable
            for lmk in eachindex(vor, vor2)
                @test abs(real(vor2[lmk])) <= abs(real(vor[lmk]))
                @test abs(imag(vor2[lmk])) <= abs(imag(vor[lmk]))
            end
        end
    end
end
