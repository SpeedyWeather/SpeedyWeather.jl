using MLDataDevices
using SpeedyWeatherInternals.Architectures: CPU, CPUStatic

@testset "MLDataDevices.get_device extension" begin
    # CPU architectures map to a CPUDevice. GPU/Reactant are covered by the same
    # method definitions but require their backend packages to be loaded, so they
    # are not exercised here to keep the unit tests fast and dependency-free.
    @test MLDataDevices.get_device(CPU()) == MLDataDevices.cpu_device()
    @test MLDataDevices.get_device(CPU()) isa MLDataDevices.CPUDevice
    @test MLDataDevices.get_device(CPUStatic()) isa MLDataDevices.CPUDevice
end
