module RingGridsNCDatasetsExt

using RingGrids
using NCDatasets

function RingGrids.get_nc_variable_name(ncfile::NCDataset, name::String)

    if !haskey(ncfile, name) && name != ""
        # Helper to show user available var names
        available = join(keys(ncfile), ", ")
        error("Variable '$name' not found in file. Available variables: [$available]")

    elseif name == ""
        candidates = filter(k -> k ∉ keys(ncfile.dim), keys(ncfile))

        if isempty(candidates)
            error("No suitable variable found in asset (all variables matched dimension names)")
        elseif length(candidates) > 1
            msg = join(candidates, ", ")
            error("Ambiguous asset: found $(length(candidates)) variables ($msg). Please specify a `varname` kwarg")
        end

        target_name = candidates[1]
        @warn "No asset variable name provided; using the only available candidate: $target_name"
        return target_name
    else
        return name
    end
end

# lazy load from NCDataset into netCDF variable (actually CommonDataModel.CFVariable)
function RingGrids._get_asset(path::String, name::String, ArrayType::Type{<:NCDataset}, FileFormat::Type{<:NCDataset})
    ds = NCDataset(path)
    # TODO also read lat, lon from file and flip array in case it's not as expected
    target_name = RingGrids.get_nc_variable_name(ds, name)
    return ds[target_name], ds
end

# load from NetCDF into Array
function RingGrids._get_asset(path::String, name::String, ArrayType::Type{<:Array}, FileFormat::Type{<:NCDataset})
    v, ds = RingGrids._get_asset(path, name, NCDataset, FileFormat)
    data = RingGrids.load_shape_preserving(v, Val(ndims(v)))
    if eltype(data) <: AbstractFloat    # exclude case of loading integer data, e.g. land-sea mask
        fill_value = get(v.attrib, "_FillValue", NaN)
        # use === to include NaNs but also expect fill value to have same type as data
        data[data .=== fill_value] .= NaN
    end
    close(ds)
    return data
end

end # module
