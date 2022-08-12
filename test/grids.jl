@testset "Grid indexing" begin
    for NF in (Float32,Float64)
        L = FullLatLonGrid(randn(NF,96*48),24)          # L24 grid
        F = FullGaussianGrid(randn(NF,96*48),24)        # F24 grid
        O = OctahedralGaussianGrid(randn(NF,3168),24)   # O24 grid
        H = HEALPixGrid(randn(NF,3072),16)              # H16 grid

        # getindex
        for ij in eachindex(L) L[ij] end
        for ij in eachindex(F) F[ij] end
        for ij in eachindex(O) O[ij] end
        for ij in eachindex(H) H[ij] end

        # set index
        for ij in eachindex(L) L[ij] = 0 end
        for ij in eachindex(F) F[ij] = 0 end
        for ij in eachindex(O) O[ij] = 0 end
        for ij in eachindex(H) H[ij] = 0 end

        @test all(L .== 0)
        @test all(F .== 0)
        @test all(O .== 0)
        @test all(H .== 0)
    end
end

@testset "Grid generators" begin
    for NF in (Float32,Float64)
        for G in (  FullLatLonGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    HEALPixGrid,
                    )

            n = 32      # resolution parameter nlat_half/nside
            G1 = zeros(G,n)
            G2 = zero(G1)
            G3 = G(G2)

            @test G1 == G2
            @test G1 == G3
            @test all(G1 .== G2)
            @test all(G1 .== G3)

            # check that G2, G3 are deep copies
            G1[1] = 1
            @test G1[1] == 1
            @test G2[1] == 0
            @test G3[1] == 0
        end
    end
end