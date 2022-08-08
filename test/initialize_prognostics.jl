@testset "Zero generators" begin
    @testset for NF in (Float32,Float64)
        P = Parameters(;NF)
        G = Geometry(P)
        S = SpectralTransform(P)
        
        P = zeros(PrognosticVariables{NF},5,5,3)
        P = zeros(DiagnosticVariables,G,S)
    end
end