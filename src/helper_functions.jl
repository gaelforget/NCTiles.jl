
"""
    getnsteps(d)

Helper function: Determines number of time steps in data d. Input d can be
BinData, NCData, TileData, or an Array. Not exported.
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
        tdim = dimnames(ds[d.varname])[end]
        timeUnits = ["minutes","seconds","hours","days","minute","second","hour","day"]
        if any(occursin.(timeUnits,Ref(lowercase(ds[tdim].attrib["units"]))))
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
matches the provided dimensions. Not exported.
"""
function checkdims(v0::Array,var::NCvar)
    dimlist = getfield.(var.dims[istimedim.(var.dims).==false],:name)
    if length(size(v0)) != length(dimlist)
        dimlist = join(dimlist,", ")
        error("Size of $(var.name) $(size(v0)) does not match its dimension list: ($dimlist)")
    end

end

"""
    istimedim(d::NCvar)

Helper function: determines whether d is a time dimension. Not exported.
"""
function istimedim(d::NCvar)

    timeUnits = ["minutes","seconds","hours","days","minute","second","hour","day"]
    return any(occursin.(timeUnits,Ref(lowercase(d.units))))

end

"""
    findtimedim(v::NCvar)

Helper function: finds which dimension is a time dimension, if any. Not exported.
"""
function findtimedim(v::NCvar)
    return findall(istimedim.(v.dims))[1]
end

"""
    hastimedim(v::NCvar)

Helper function: determines whether a variable has a time dimension. Not exported.
"""
function hastimedim(v::NCvar)
    ncvardim = isa.(v.dims,NCvar)
    return any(ncvardim) && any(istimedim.(v.dims[ncvardim]))
end

"""
    hastimedimdims::Array{NCvar})

Helper function: determines whether an array of dimensions has a time dimension. Not
    exported.
"""
function hastimedim(dims::Array{NCvar})
    return any(istimedim.(dims))
end