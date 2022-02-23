@testset "Initialize from rest" begin

    # original speedy T30 with 96x48 grid
    P = Parameters(NF=Float64)
    G = GeoSpectral(P)
    B = Boundaries(P,G)

    Prog = SpeedyWeather.initialize_from_rest(P,B,G)

    @test all(Prog.vor .== 0)
    @test all(Prog.div .== 0)

    # test surface layer only at the moment
    temp_grid = gridded(Prog.temp[:,:,end],G)
    pres_surf_grid = gridded(Prog.pres_surf,G)
    humid_grid = gridded(Prog.humid[:,:,end],G)

    # temperature between 200K and 350K everywhere
    @test all(temp_grid .> 200)
    @test all(temp_grid .< 350)

    # surface pressure between log(300hPa) and log(1200hPa) everywhere
    @test all(pres_surf_grid .> log(300))
    @test all(pres_surf_grid .< log(1200))

    # humidity non-negative everywhere
    # humidity has currently values of O(1e9)...
    @test_skip all(humid_grid .>= 0) 
end
