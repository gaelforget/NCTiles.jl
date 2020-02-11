
"""
    getnsteps(d)

Helper function: Determines number of time steps in data d. Input d can be
BinData, NCData, TileData, or an Array.
"""
function getnsteps(d)
    if isa(d,BinData)
        if isa(d.fnames,Array)
            nsteps = length(d.fnames)
        else
            nsteps = 1
        end
    elseif isa(d,NCData)
        ds = Dataset(d.fname)
        dsvar = ds[d.varname]
        if hastimedim(ds,dsvar)
            tdim = dimnames(dsvar)[findtimedim(ds,dsvar)]
            nsteps = length(ds[tdim])
        else
            nsteps = 1
        end
    elseif isa(d,TileData)
        nsteps = getnsteps(d.vals)
    elseif isa(d,Array) && isa(d[1],Array)
        nsteps = length(d)
    else
        nsteps = 1
    end

    return nsteps
end

"""
    checkdims(v0::Array,var::NCvar)

Helper function: Checks that the size of the data about to be written to the file
matches the provided dimensions.
"""
function checkdims(v0::Array,var::NCvar)
    if isa(var.values,Array) && isa(var.values[1],Number)
        dimlist = getfield.(var.dims,:name)
    else
        dimlist = getfield.(var.dims[istimedim.(var.dims).==false],:name)
    end
    if length(size(v0)) != length(dimlist)
        dimlist = join(dimlist,", ")
        error("Size of $(var.name) $(size(v0)) does not match its dimension list: ($dimlist)")
    end

end

"""
    istimedim(d::Union{NCvar,NCDatasets.CFVariable})

Helper function: determines whether d is a time dimension.
"""
function istimedim(d::Union{NCvar,NCDatasets.CFVariable})
    if isa(d,NCvar)
        units = lowercase(d.units)
        longname = get(d.atts,"long_name","")
    else
        units = lowercase(d.attrib["units"])
        longname = get(d.attrib,"long_name","")
    end
    timeUnits = ["minutes","seconds","hours","days","minute","second","hour","day"]
    return any(occursin.(timeUnits,Ref(units))) || longname=="Time coordinate"
end

"""
    findtimedim(v::Array)

Helper function: finds which dimension is a time dimension, if any. Can be Array of
    NCvars or NCDatasets.CFVariables.
"""
function findtimedim(dims::Array)
    return findall(istimedim.(dims))[1]
end

findtimedim(v::NCvar) = findtimedim(v.dims)
findtimedim(ds::NCDatasets.Dataset,v::NCDatasets.CFVariable) = findtimedim([ds[d] for d in dimnames(v)])

"""
    hastimedim::Array)

Helper function: determines whether an array of dimensions has a time dimension. Can 
    be Array of NCvars or NCDatasets.CFVariables.
"""
function hastimedim(dims::Array)
    return any(istimedim.(dims))
end

"""
    hastimedim(v::NCvar)

Helper function: determines whether a variable has a time dimension.
"""
function hastimedim(v::NCvar)
    ncvardim = isa.(v.dims,NCvar)
    return any(ncvardim) && hastimedim(v.dims[ncvardim])
end

hastimedim(ds::NCDatasets.Dataset,v::NCDatasets.CFVariable) = hastimedim([ds[d] for d in dimnames(v)])

"""
    getdims(myflds::Dict)

Helper function: retrieves unique dimensions in dictionary of fields.
"""
function getdims(myflds::Dict)
    dims = unique(vcat([myflds[v].dims for v in keys(myflds)]...))
    dims = filter( d -> isa(d,NCvar),dims)
end

"""
    hastiledata(myfld::NCvar)

Helper function: returns whether given field or dict of fields has any tiled data.
"""
function hastiledata(myfld::NCvar)
    return isa(myfld.values,TileData)
end

hastiledata(myflds::Dict) = any([hastiledata(myflds[f]) for f in keys(myflds)])