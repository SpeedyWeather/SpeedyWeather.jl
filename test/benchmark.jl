using Test, BenchmarkTools
include("../src/SpeedyWeather.jl")
using .SpeedyWeather

# original speedy T30 with 96x48 grid
P = Parameters(NF=Float64;trunc=30,nlat=48,nlon=96)
G = GeoSpectral(P)
B = Boundaries(P,G)

geopot_surf_spectral = B.geopot_surf
geopot_surf_grid = gridded(geopot_surf_spectral,G)
geopot_surf_fourier = fourier(geopot_surf_grid,G)

@btime fourier_inverse(fourier($geopot_surf_grid,$G),$G);
@btime legendre_inverse(legendre($geopot_surf_fourier,$G),$G);
@btime gridded(spectral($geopot_surf_grid,$G),$G);