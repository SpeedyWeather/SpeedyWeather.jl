@testset "Horizontal diffusion of random" begin
    
    mx,nx,nlev = 31,32,8

    for T in (Float16,Float32,Float64)

        A           = rand(Complex{T},mx,nx,nlev)
        tendency    = rand(Complex{T},mx,nx,nlev)
        D1          = rand(T,mx,nx)
        D2          = rand(T,mx,nx)

        A1 = A[:,:,1]
        tendency1 = copy(tendency[:,:,1])

        SpeedyWeather.horizontal_diffusion!(A,tendency,D1,D2)
        SpeedyWeather.horizontal_diffusion!(A1,tendency1,D1,D2)

        @test tendency1 == tendency[:,:,1]
    end
end