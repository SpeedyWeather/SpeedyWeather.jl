@testset "Initialize from rest" begin

    # original speedy T30 with 96x48 grid
    P = Parameters(NF=Float64)
    G = GeoSpectral(P)
    B = Boundaries(P)
    S = G.spectral

    Prog = SpeedyWeather.initialize_from_rest(P,B,G)

    @test all(Prog.vor .== 0)
    @test all(Prog.div .== 0)

    # test surface layer only at the moment
    temp_grid = gridded(Prog.temp[:,:,end],S)
    pres_surf_grid = gridded(Prog.pres_surf,S)
    humid_grid = gridded(Prog.humid[:,:,end],S)

    # temperature between 200K and 350K everywhere
    println((sum(temp_grid)/length(temp_grid),minimum(temp_grid),maximum(temp_grid)))
    @test all(temp_grid .> 200)
    @test all(temp_grid .< 350)

    # surface pressure between log(300hPa) and log(2000hPa) everywhere
    println((sum(pres_surf_grid)/length(pres_surf_grid),minimum(pres_surf_grid),maximum(pres_surf_grid)))
    @test all(pres_surf_grid .> log(300))
    @test all(pres_surf_grid .< log(2000))

    # humidity non-negative everywhere
    # humidity has currently values of O(1e9)...
    println((sum(humid_grid)/length(humid_grid),minimum(humid_grid),maximum(humid_grid)))
    @test_skip all(humid_grid .>= 0) 
end
