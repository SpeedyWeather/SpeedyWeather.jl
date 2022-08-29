using KernelAbstractions    

@testset "KernelAbstractions tests" begin 

    # only on CPU (currently)

    # Tests for the basic setup 
    device_setup = SpeedyWeather.DeviceSetup(SpeedyWeather.CPUDevice())

    @kernel function mul_test!(A, @Const(B), @Const(C))
        i, j = @index(Global, NTuple)
        A[i,j] = B[i,j] * C[i,j]
    end
    
    B = SpeedyWeather.DeviceArray(device_setup, rand(10,10))
    C = SpeedyWeather.DeviceArray(device_setup, rand(10,10))
    A = zero(B)

    ev = SpeedyWeather.launch_kernel!(device_setup, mul_test!, size(B), A, B, C)
    wait(ev)

    @test A ≈ B .* C 

    # ----------------------------------
    # Tests for the individual functions and functions of KernelAbstractions

    # first some preperation 
    progn_vars, diagn_vars, model_setup = initialize_speedy();
    SpeedyWeather.time_stepping!(progn_vars,diagn_vars,model_setup) # evolve the state of speedy, so that we also have nonzero entries for everything
    
    M = model_setup 
    diagn = diagn_vars.layers[1]
    progn = progn_vars.layers[1] 
    lf = 1
    NF = Float32
    device = model_setup.device_setup

    # ---------------------------------
    # transferring all arrays to the device
    diagn_vars = SpeedyWeather.DeviceArray(device, diagn_vars)
    progn_vars = SpeedyWeather.DeviceArray(device, progn_vars)
    model_setup = SpeedyWeather.DeviceArray(device, model_setup)

    # ---------------------------------
    # now the individual functions 

    # horizontal_diffusion! 
    SpeedyWeather.get_tendencies!(diagn, M)   

    (;vor) = progn.leapfrog[lf]
    (;vor_tend) = diagn.tendencies
    (;damping, damping_impl) = M.horizontal_diffusion

    vor_tend_old = deepcopy(vor_tend)
    vor_tend_new = deepcopy(vor_tend)

    SpeedyWeather.horizontal_diffusion!(vor_tend_old, vor, damping, damping_impl)
    SpeedyWeather.horizontal_diffusion!(vor_tend_new, vor, damping, damping_impl, model_setup.device_setup)
    @test vor_tend_old ≈ vor_tend_new
    # ----------------------------------

end 