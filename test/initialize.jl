@testset "Zero generators" begin
    @testset for NF in (Float32,Float64)
        P = Parameters(;NF)
        G = Geometry(P)
        S = SpectralTransform(P)
        
        P = zeros(PrognosticVariables{NF},5,5,3)
        P = zeros(DiagnosticVariables,G,S)
    end
end

@testset "Initialize from rest" begin

    # BAROTROPIC MODEL
    progn, diagn, model = initialize_speedy(Barotropic,initial_conditions=:rest)
    for layers in progn.layers
        for leapfrog in layers.leapfrog
            @test all(leapfrog.vor .== 0)
        end
    end

    # SHALLOW WATER MODEL
    progn, diagn, model = initialize_speedy(ShallowWater,initial_conditions=:rest)
    for layers in progn.layers
        for leapfrog in layers.leapfrog
            @test all(leapfrog.vor .== 0)
            @test all(leapfrog.div .== 0)
        end
    end
    @test all(progn.pres.leapfrog[1] .== 0)
    @test all(progn.pres.leapfrog[2] .== 0)

    # # PRIMITIVE EQUATION MODEL
    # progn, diagn, model = initialize_speedy(initial_conditions=:rest,model=PrimitiveEquation)
    # @test all(progn.vor .== 0)

    """
    S = model.geospectral.spectral_transform
    k = model.parameters.nlev       # test surface layer only at the moment
    lf = 1                          # first leapfrog index
    temp_grid = gridded(progn.temp[:,:,lf,k],S)
    pres_surf_grid = gridded(progn.pres_surf[:,:,lf],S)
    humid_grid = gridded(progn.humid[:,:,lf,k],S)

    # temperature between 200K and 350K everywhere
    # println((sum(temp_grid)/length(temp_grid),minimum(temp_grid),maximum(temp_grid)))
    @test all(temp_grid .> 200)
    @test all(temp_grid .< 350)

    # surface pressure between log(300hPa) and log(2000hPa) everywhere
    # println((sum(pres_surf_grid)/length(pres_surf_grid),minimum(pres_surf_grid),maximum(pres_surf_grid)))
    @test all(pres_surf_grid .> log(300))
    @test all(pres_surf_grid .< log(2000))

    # humidity non-negative everywhere
    # humidity has currently values of O(1e9)...
    # println((sum(humid_grid)/length(humid_grid),minimum(humid_grid),maximum(humid_grid)))
    @test_skip all(humid_grid .>= 0)
    """
end