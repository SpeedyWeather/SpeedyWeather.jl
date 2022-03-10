using Test, BenchmarkTools
include("../src/SpeedyWeather.jl")
using .SpeedyWeather

# original speedy T30 with 96x48 grid
trunc = 31
nlon = 96
nlat = 48

P = Parameters(NF=Float64;trunc,nlat,nlon)
G = GeoSpectral(P)
B = Boundaries(P,G)
S = G.spectral

geopot_surf_spectral = B.geopot_surf
geopot_surf_grid = gridded(geopot_surf_spectral)
geopot_surf_spectral2 = spectral(geopot_surf_grid)
geopot_surf_grid2 = gridded(geopot_surf_spectral2)

spectral_truncation!(geopot_surf_spectral)
spectral_truncation!(geopot_surf_spectral2)

@btime spectral!($geopot_surf_spectral,$geopot_surf_grid,$S);
@btime gridded!($geopot_surf_grid,$geopot_surf_spectral,$S);