using KernelAbstractions    
import Parameters: @with_kw, @unpack


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
    # Tests for the individual functions 

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
    # now the individual functions 

    # horizontal_diffusion! 
    SpeedyWeather.get_tendencies!(diagn, M)   

    @unpack vor = progn.leapfrog[lf]
    @unpack vor_tend = diagn.tendencies
    @unpack damping, damping_impl = M.horizontal_diffusion
                                      
    # we have to convert everything to CuArrays, incase we are on GPU, later we will not have to do this anymore when we have fitting `adapt` defined for all structs 
    vor_tend = SpeedyWeather.DeviceArray(device, vor_tend)
    vor = SpeedyWeather.DeviceArray(device, vor)
    damping = SpeedyWeather.DeviceArray(device, damping)
    damping_impl = SpeedyWeather.DeviceArray(device, damping_impl)

    vor_tend_old = deepcopy(vor_tend)
    vor_tend_new = deepcopy(vor_tend)

    SpeedyWeather.horizontal_diffusion!(vor_tend_old, vor, damping, damping_impl)
    SpeedyWeather.horizontal_diffusion!(vor_tend_new, vor, damping, damping_impl, model_setup.device_setup)
    @test vor_tend_old ≈ vor_tend_new
    # ----------------------------------

end 