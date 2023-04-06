@testset "Horizontal diffusion of random" begin
    for T in (Float32,Float64)

        p,d,m = initialize_speedy(T)

        (;vor) = p.layers[1].timesteps[1]
        (;vor_tend) = d.layers[1].tendencies
        (;∇²ⁿ,∇²ⁿ_implicit) = m.horizontal_diffusion

        vor = randn(typeof(vor),size(vor)...)

        SpeedyWeather.horizontal_diffusion!(vor_tend,vor,∇²ⁿ,∇²ⁿ_implicit)

        # diffusion tendency has opposite sign (real/imag respectively)
        # than prognostic variable to act as a dissipation 
        (;lmax,mmax) = m.spectral_transform
        for m in 1:mmax+1
            for l in max(2,m):lmax+1
                @test -sign(real(vor[l,m])) == sign(real(vor_tend[l,m]))
                @test -sign(imag(vor[l,m])) == sign(imag(vor_tend[l,m]))
            end
        end

        vor0 = copy(vor)
        vor1 = copy(vor)
        SpeedyWeather.leapfrog!(vor0,vor1,vor_tend,m.constants.Δt,1,m.constants)

        @test any(vor0 .!= vor1)    # check that at least some coefficients are different
        @test any(vor0 .== vor1)    # check that at least some coefficients are identical

        # damping should not increase real or imaginary part of variable
        for lm in SpeedyWeather.eachharmonic(vor0,vor)
            @test abs(real(vor1[lm])) <= abs(real(vor[lm]))
            @test abs(imag(vor1[lm])) <= abs(imag(vor[lm]))
        end
    end
end