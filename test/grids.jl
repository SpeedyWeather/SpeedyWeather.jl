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