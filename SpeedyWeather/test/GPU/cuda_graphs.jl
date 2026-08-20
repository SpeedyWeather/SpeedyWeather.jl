ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsCUDAExt)
test_gpu_graphs(ext, "CUDA")
