@testset "Horizontal diffusion of random" begin
    for NF in (Float32, Float64)

        spectral_grid = SpectralGrid(NF)
        m = PrimitiveDryModel(; spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        (; vor) = p.layers[1].timesteps[1]
        (; vor_tend) = d.layers[1].tendencies
        (; ∇²ⁿ, ∇²ⁿ_implicit) = m.horizontal_diffusion

        vor = randn(typeof(vor), size(vor)...)

        # ∇²ⁿ, ∇²ⁿ_implicit are nlev-vectors, one array per layer, pick surface
        SpeedyWeather.horizontal_diffusion!(vor_tend, vor, ∇²ⁿ[end], ∇²ⁿ_implicit[end])

        # diffusion tendency has opposite sign (real/imag respectively)
        # than prognostic variable to act as a dissipation 
        (; lmax, mmax) = m.spectral_transform
        for m in 1:mmax+1
            for l in max(2, m):lmax
                @test -sign(real(vor[l, m])) == sign(real(vor_tend[l, m]))
                @test -sign(imag(vor[l, m])) == sign(imag(vor_tend[l, m]))
            end
        end

        vor0 = copy(vor)
        vor1 = copy(vor)
        SpeedyWeather.leapfrog!(vor0, vor1, vor_tend, m.time_stepping.Δt, 1, m.time_stepping)

        @test any(vor0 .!= vor1)    # check that at least some coefficients are different
        @test any(vor0 .== vor1)    # check that at least some coefficients are identical

        # damping should not increase real or imaginary part of variable
        for lm in SpeedyWeather.eachharmonic(vor0, vor)
            @test abs(real(vor1[lm])) <= abs(real(vor[lm]))
            @test abs(imag(vor1[lm])) <= abs(imag(vor[lm]))
        end
    end
end