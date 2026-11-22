ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsAMDGPUExt)
test_gpu_graphs(ext, "HIP")
