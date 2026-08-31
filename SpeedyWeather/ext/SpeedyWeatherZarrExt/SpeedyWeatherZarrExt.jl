module SpeedyWeatherZarrExt

using SpeedyWeather
using Zarr
using DocStringExtensions

import SpeedyWeather: ZarrOutput, HEALPixOutput, AbstractZarrOutput,
    AbstractOutput, AbstractOutputVariable,
    AbstractSimulation, AbstractModel, Barotropic, OutputWriterCore,
    OUTPUT_VARIABLES_DICT, OutputVariablesDict, DEFAULT_NLAYERS_SOIL,
    DEFAULT_OUTPUT_NF, DEFAULT_OUTPUT_INTERVAL, DEFAULT_MISSING_VALUE,
    DEFAULT_COMPRESSION_LEVEL, DEFAULT_KEEPBITS,
    Variables, Simulation, SpectralGrid, Field,
    initialize!, finalize!, output!, write_array!, define_variable!, set!, add!, add_default!,
    define_dimension!, vertical_dimension, get_nlayers, get_dimension_length, define_coordinate!,
    is3D, is_land, hastime, get_indices, get_flat_indices, scale!, get_soil_layers,
    get_lond, get_latd, get_npoints, get_nlat, on_architecture, CPU,
    AbstractFullGrid, path_or_nothing, run_folder_name

import SpeedyWeather.RingGrids
import SpeedyWeather.RingGrids: HEALPixGrid, get_londlatds, whichring, grids_match
import SpeedyWeather: round!
import SpeedyWeather.Printf
import SpeedyWeather.Dates: Dates, DateTime, Period, Second, Millisecond

include("shared.jl")            # utilities shared by all Zarr-backed writers
include("zarr_output.jl")       # ZarrOutput: rectangular lon/lat layout
include("healpix_output.jl")    # HEALPixOutput: flat HEALPix (RING) layout

end # module
