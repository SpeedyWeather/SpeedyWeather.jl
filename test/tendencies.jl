using SpeedyWeather
using Test

# @testset "Diagnostic variables (tendencies) initialisation." begin
#     # loop over different precisions
#     for NF in (Float16,Float32,Float64)

#         P = Parameters(NF=NF)
#         G = GeoSpectral{NF}(P)
#         Diag = DiagnosticVariables{NF}(G)

#         @test all(Diag.tendencies.vor_tend .== 0) #Check zero valued.
#         @test all(Diag.grid_variables.vor_grid.== 0)

#     end
# end

@testset "Calculate Tendencies" begin
    NF = Float64 

    P,D,M = initialize_model(NF)

    l2 = 2

    get_tendencies!(P,D,l2,M) 

    # P = Parameters(NF=NF)
    # C = Constants{NF}(P)
    # G = GeoSpectral{NF}(P)
    # B = Boundaries{NF}(P,G)
    # Diag = DiagnosticVariables{NF}(G)

    # #@test all(Diag.tendencies.vor_tend .== 0) #Check zero valued
    
   
    # Prog = initial_conditions(P,B,G)
       
    # 


    
end