@testset "Horizontal diffusion of random" begin
    for NF in (Float32, Float64)

        spectral_grid = SpectralGrid(; NF)
        model = PrimitiveDryModel(; spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables

        (; ∇²ⁿ, ∇²ⁿ_implicit) = model.horizontal_diffusion
        vor = progn.vor[1]
        (; vor_tend) = diagn.tendencies

        vor .= randn(LowerTriangularArray{eltype(vor)}, size(vor, as=Matrix)...)

        # ∇²ⁿ, ∇²ⁿ_implicit are nlev-vectors, one array per layer, pick surface
        SpeedyWeather.horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)

        # diffusion tendency has opposite sign (real/imag respectively)
        # than prognostic variable to act as a dissipation 
        (; lmax, mmax) = model.spectral_transform
        for k in eachmatrix(vor, vor_tend)
            for m in 1:mmax+1
                for l in max(2, m):lmax
                    @test -sign(real(vor[l, m, k])) == sign(real(vor_tend[l, m, k]))
                    @test -sign(imag(vor[l, m, k])) == sign(imag(vor_tend[l, m, k]))
                end
            end
        end

        vor2 = progn.vor[2]
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