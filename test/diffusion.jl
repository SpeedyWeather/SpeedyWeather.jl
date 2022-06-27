@testset "Horizontal diffusion of random" begin
    for T in (Float32,Float64)

        p,d,m = initialize_speedy(T)

        (;vor) = p
        (;vor_tend) = d.tendencies
        (;damping,damping_impl) = m.horizontal_diffusion

        vor       .= randn(Complex{T},size(vor)...)
        vor_tend  .= zeros(Complex{T},size(vor_tend)...)

        vor_lf1 = view(vor,:,:,1,:)

        SpeedyWeather.horizontal_diffusion!(vor_tend,vor_lf1,damping,damping_impl)

        (;nlev) = m.geospectral.geometry
        (;lmax,mmax) = m.geospectral.spectral_transform

        # diffusion tendency has opposite sign (real/imag respectively)
        # than prognostic variable to act as a dissipation 
        for k in [1]
            for m in 1:mmax+1
                for l in max(2,m):lmax+1
                    @test -sign(real(vor_lf1[l,m,k])) == sign(real(vor_tend[l,m,k]))
                    @test -sign(imag(vor_lf1[l,m,k])) == sign(imag(vor_tend[l,m,k]))
                end
            end
        end

        vor0 = copy(vor_lf1)
        SpeedyWeather.leapfrog!(vor,vor_tend,m.constants.Δt,m.constants,1)
        vor1 = view(vor,:,:,2,:)

        @test any(vor0 .!= vor1)    # check that at least some coefficients are different
        @test any(vor0 .== vor1)    # check that at least some coefficients are identical

        # damping should not increase real or imaginary part of variable
        for lm in eachindex(vor0,vor_lf1)
            @test real(vor0[lm]) >= real(vor_lf1[lm])
            @test imag(vor0[lm]) >= imag(vor_lf1[lm])
        end
    end
end