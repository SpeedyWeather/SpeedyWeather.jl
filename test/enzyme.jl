using Enzyme, KernelGradients
import Parameters: @with_kw, @unpack
 

@testset "AD/Enzyme tests" begin 
    # Initizialize what we need from SpeedyWeather 
    progn_vars, diagn_vars, model_setup = initialize_speedy();
    SpeedyWeather.time_stepping!(progn_vars,diagn_vars,model_setup) # evolve the state of speedy, so that we also have nonzero entries for everything

    M = model_setup 
    diagn = diagn_vars.layers[1]
    progn = progn_vars.layers[1] 
    NF = Float32
    device = model_setup.device_setup

    # ----------------------------------
    # Tests for the individual functions 

    # horizontal_diffusion! 
    SpeedyWeather.get_tendencies!(diagn, M)   
    lf = 1
    @unpack vor = progn.leapfrog[lf]
    @unpack vor_tend = diagn.tendencies
    @unpack damping, damping_impl = M.horizontal_diffusion
                                      
    # we have to convert everything to CuArrays, incase we are on GPU, later we will not have to do this anymore when we have fitting `adapt` defined for all structs 
    vor_tend = SpeedyWeather.DeviceArray(device, vor_tend)
    vor = SpeedyWeather.DeviceArray(device, vor)
    damping = SpeedyWeather.DeviceArray(device, damping)
    damping_impl = SpeedyWeather.DeviceArray(device, damping_impl)

    ∂vor_tend = one(vor_tend)
    ∂vor = zero(vor)
    ∂damping = zero(damping)
    ∂damping_impl = zero(damping_impl);

    ∇! = autodiff(SpeedyWeather.horizontal_diffusion_kernel!(model_setup.device_setup.device_KA(), model_setup.device_setup.n))
    ev = ∇!(Duplicated(vor_tend, ∂vor_tend), Duplicated(vor, ∂vor), Duplicated(damping, ∂damping), 
        Duplicated(damping_impl, ∂damping_impl); ndrange=length(vor_tend))
    wait(ev)

    @test true # right now everything just running without error has to be sufficient
    # ----------------------------------
end