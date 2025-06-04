import SpeedyWeather.SpectralTransforms: ismatching, _divergence!

# _divergence! function barrier for GPU arrays 
SpeedyWeather.SpectralTransforms._divergence!(
    kernel,
    div::LowerTriangularArray{T,N,<:CuArray},
    u::LowerTriangularArray{T,N,<:CuArray},
    v::LowerTriangularArray{T,N,<:CuArray},
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
) = _divergence_KA!(kernel, div, u, v, S; radius)
    
