@testset "Initialize from rest" begin

    # original speedy T30 with 96x48 grid
    P = Params(NF=Float64,trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    Prog = initialise_from_rest(P,B,G)

    @test all(Prog.vor .= 0)
    @test all(Prog.div .= 0)

    temp_grid = gridded(Prog.temp,G)
    pres_surf_grid = gridded(Prog.pres_surf,G)
    humid_grid = gridded(Prog.humid,G)

    # temperature between 100K and 350K everywhere
    @test all(temp_grid .> 100)
    @test all(temp_grid .< 350)

    # pressure between log(800hPa) and log(1200hPa) everywhere
    @test all(pres_surf_grid .> log(800))
    @test all(pres_surf_grod .< log(1200))

    # humidity non-negative everywhere
    @test all(humid_grid .>= 0) 
end