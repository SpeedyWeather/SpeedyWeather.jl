module RingGridsNCDatasetsExt

using RingGrids
using NCDatasets

function get_nc_variable_name(ncfile::NCDataset, name::String)

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
function RingGrids._get_asset(path::String, name::String, ArrayType::Type{<:NCDataset}, FileFormat::Type{<:NCDataset}, fill_value)
    ds = NCDataset(path)
    # TODO also read lat, lon from file and flip array in case it's not as expected
    target_name = get_nc_variable_name(ds, name)
    return ds[target_name], ds
end

# load from NetCDF into Array
function RingGrids._get_asset(path::String, name::String, ArrayType::Type{<:Array}, FileFormat::Type{<:NCDataset}, fill_value)
    v, ds = RingGrids._get_asset(path, name, NCDataset, FileFormat, fill_value)
    data = RingGrids.load_shape_preserving(v, Val(ndims(v)))
    if eltype(data) <: AbstractFloat    # exclude case of loading integer data, e.g. land-sea mask
        nc_fill = get(v.attrib, "_FillValue", fill_value)
        # use === to include NaNs but also expect fill value to have same type as data
        data[data .=== nc_fill] .= fill_value
    end
    close(ds)
    return data
end

end # module
