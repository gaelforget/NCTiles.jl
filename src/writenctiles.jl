using NCDatasets,NetCDF,Dates,MeshArrays,Printf

"""
    NCvar

Data structure containing information needed to write a NetCDF file. This
includes a list of filenames (see `Bindata`) if the data is not loaded into memory.
"""
struct NCvar
    name::String
    units::String
    dims
    values
    atts::Union{Dict,Nothing}
    backend::Module
end

function NCvar(name::String, units::String, dims::NCvar, values, atts::Union{Dict,Nothing}, backend::Module)
    return NCvar(name, units, [dims], values, atts, backend)
end

"""
    BinData

Data structure containing a string or an array of strings (NetCDF
    file names) as well as metadata needed to read a file.
"""
struct BinData # Pointer to data stored in binary files- contains info needed to read in
    fnames::Union{Array{String},String}
    precision::Type
    iosize::Tuple
    fldidx::Int
end

"""
    BinData(fnames::Union{Array{String},String},precision::Type,iosize::Tuple)

Construct a BinData struct for files that contain one field.
"""
function BinData(fnames::Union{Array{String},String},precision::Type,iosize::Tuple)
    return BinData(fnames,precision,iosize,1)
end

function BinData(fnames::Union{Array{String},String},precision::Type,iosize::Tuple,meta::Dict,fldname::String)
    return BinData(fnames,precision,iosize,findfirst(meta["fldList"] .== fldname)[2])
end

"""
    NCData

Data structure containing a string or an array of strings (file names) of
    NetCDF files as well as information needed to read a file.
"""
struct NCData
    fname::AbstractString
    varname::AbstractString
    backend::Module
    precision::Type
end


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
    TileData{T}

Data structure containing either a `MeshArray` struct or `BinData` struct (see `vals`),
    `MeshArray` structs describing the tile layout (`tileinfo`), and other information for
    reading/writing tile data.
"""
struct TileData{T}
    vals::T
    tileinfo::Dict
    tilesize::Tuple
    precision::Type
    numtiles::Int
end

"""
    TileData(vals,tilesize::Tuple)

Construct a TileData struct. First generate the tileinfo, precision, and numtiles attributes.
"""
function TileData(vals,tilesize::Tuple,grid::gcmgrid)
    ni,nj = tilesize
    tileinfo = findtiles(ni,nj,grid)
    if isa(vals,BinData)
        prec = vals.precision
    else
        prec = Float32
    end
    return TileData(vals,tileinfo,tilesize,Float32,Int(maximum(tileinfo["tileNo"])))
end
TileData(vals,tilesize::Tuple,grid::String="LLC90") = TileData(vals,tilesize::Tuple,GCMGridSpec(grid))


"""
    findidx(A,val)

Helper function for getting the indices for tiles. A is a gcmfaces struct, val is a
    numeric value to get the indices of. Returns a Dict of the indices for each face.
    Not currently exported. Maybe move to MeshArrays?
"""
function findidx(A,val)
    tileidx = Dict()
    if isa(A,MeshArrays.gcmfaces)
        nfaces = A.nFaces
    elseif isa(A,MeshArray)
        nfaces = A.grid.nFaces
    end
    for i = 1:nfaces
        idx = findall(x->x==val,A.f[i])
        tileidx[i] = [x for x in idx]
    end
    return tileidx
end

"""
    gettile(fldvals,tileinfo,tilesize,tilenum::Int)

Helper function for retrieving a tile from a gcmfaces struct as a numeric Array. Not
    currently exported.
"""
function gettile(fldvals,tileinfo,tilesize,tilenum::Int)
    tilidx = findidx(tileinfo["tileNo"],tilenum)

    if isa(fldvals,MeshArrays.gcmfaces)
        is3D = length(size(fldvals)) == 3
        if is3D; n3 = size(fldvals)[3]; end
    else
        is3D = length(size(fldvals)) == 2 && size(fldvals)[2] > 1
        if is3D; n3 = size(fldvals)[2]; end
    end

    if is3D # has depth
        tilfld = Array{Any,3}(nothing,0,tilesize[2],n3)
    else
        tilfld = Array{Any,2}(nothing,0,tilesize[2])
    end
    if isa(fldvals,MeshArrays.gcmfaces)
        nfaces = fldvals.nFaces
    elseif isa(fldvals,MeshArray)
        nfaces = fldvals.grid.nFaces
    end
    for iF = 1:nfaces
        if ~isempty(tilidx[iF])
            imin = minimum(tilidx[iF])[1]; imax = maximum(tilidx[iF])[1]
            jmin = minimum(tilidx[iF])[2]; jmax = maximum(tilidx[iF])[2]
            if is3D
                if isa(fldvals,MeshArray)
                    tilfld = [tilfld; cat([fldvals.f[iF,d][imin:imax,jmin:jmax] for d in 1:n3]...,dims=3)]
                else
                    tilfld = [tilfld; fldvals.f[iF][imin:imax,jmin:jmax,:]]
                end
            else
                tilfld = [tilfld; fldvals.f[iF][imin:imax,jmin:jmax]]
            end
        end
    end

    if is3D
        tilfld = reshape(tilfld,tilesize[1],tilesize[2],:)
    else
        tilfld = reshape(tilfld,tilesize[1],tilesize[2])
    end

    return tilfld
end

"""
    gettiles(fldvals,tilenum::Int)

Helper function for retrieving a tile from a gcmfaces struct as a numeric Array along
    with associated latitude and longitude. Not currently exported.
"""
function gettiles(tilfld,tilenum::Int)
    tilesize = tilfld.tilesize

    tilfld = gettile(tilfld.vals,tilfld.tileinfo,tilesize,tilenum)
    tillat = gettile(tilfld.tileinfo["XC"],tilfld.tileinfo,tilfld.tilesize,tilenum)
    tillon = gettile(tilfld.tileinfo["YC"],tilfld.tileinfo,tilfld.tilesize,tilenum)

    return tilfld,tillat,tillon
end

"""
    getindex(v::NCvar,i::Integer)

Retrive an NCvar at time step i.

Currently only works when values is a BinData, NCdata, or array of time steps.
"""
function getindex(v::NCvar,i::Integer)
    NCvar(v.name,v.units,v.dims,v.values[i],v.atts,v.backend)
end

"""
    getindex(d::BinData,i::Integer)

 Return a new BinData object with the ith file, this being time index i.
"""
function getindex(d::BinData,i::Integer) # Gets file at time index i
    if isa(d.fnames,Array)
        newfnames = [d.fnames[i]]
    elseif i == 1
        newfnames = d.fnames
    end

    BinData(newfnames,d.precision,d.iosize)
end

"""
    getindex(d::NCData,i)

Retrieve the ith time step from the NetCDF file.
"""
function getindex(d::NCData,i) # Gets data at time index i
    readncdata(d,i)
end

"""
    readbin(fname::String,prec::Type,iosize::Tuple,fldidx=1)

Read in a binary file as an array as previously done via `MeshArrays.read_bin`
"""
function readbin(fname::String,prec::Type,iosize::Tuple,fldidx=1)
    n3 = 1
    if length(iosize) == 3
        n1,n2,n3 = iosize
    else
        n1,n2 = iosize
    end

    if prec == Float64
        reclen = 8
    else
        reclen = 4
    end

    #if isempty(n3)
    #    n3=Int64(stat(fname).size/n1/n2/reclen)
    #end

    field = Array{prec,1}(undef,(n1*n2*n3))
    fid = open(fname)
    if fldidx != 1
        seek(fid,n1*n2*n3*(fldidx-1)*reclen)
    end
    read!(fid,field)
    close(fid)
    field = hton.(field)

    n3>1 ? s=(n1,n2,n3) : s=(n1,n2)
    return reshape(field,s)

end

"""
    readbin(flddata::BinData,tidx=1)

Read in a binary file as an array as previously done via `MeshArrays.read_bin`
"""
readbin(flddata::BinData,tidx=1) = readbin(flddata.fnames[tidx],flddata.precision,flddata.iosize,flddata.fldidx)
"""
    readncdata(var::NCData,i::Union{Colon,Integer}=:)

Read netcdf file as specified in `NCData` argument. Optional
argument `i` can be used to read a specific records / times.
"""
function readncdata(var::NCData,i::Union{Colon,Integer}=:)
    ds = Dataset(var.fname)
    ndims = length(size(ds[var.varname]))

    if ndims == 1
        values = ds[var.varname][i]
    elseif ndims == 2
        values = ds[var.varname][:,i]
    elseif ndims == 3
        values = ds[var.varname][:,:,i]
    elseif ndims == 4
        values = ds[var.varname][:,:,:,i]
    else
        println("Variables greater than four dimensions not currently supported")
    end
    close(ds)
    return values
end

function readdata(flddata,tidx=1)
    if isa(flddata,BinData)
        res = readbin(flddata,tidx)#readbin(flddata.fnames[tidx],flddata.precision,flddata.iosize,flddata.fldidx)
    elseif isa(flddata,NCData)
        res = readncdata(flddata,tidx)
    end
    return res
end


"""
    addDim(ds::NCDatasets.Dataset,dimvar::NCvar) # NCDatasets

Add a dimension to a `NCDatasets.Dataset`
"""
function addDim(ds::NCDatasets.Dataset,dimvar::NCvar)
    defDim(ds,dimvar.name,dimvar.dims[1])
end

"""
    addDim(dimvar::NCvar)

Add a dimension to a NetCDF file using `NetCDF.jl`
"""
function addDim(dimvar::NCvar) #NetCDF
    NcDim(dimvar.name,collect(dimvar.values),
    atts = merge(Dict(("units" =>dimvar.units)),dimvar.atts),
    unlimited = dimvar.dims==Inf)
end

"""
    addVar(ds::NCDatasets.Dataset,field::NCvar)

Add a variable to a NetCDF file using `NCDatasets.jl`
"""
function addVar(ds::NCDatasets.Dataset,field::NCvar)
    if ~isempty(field.units)
        atts = merge(Dict(("units" =>field.units)),field.atts)
    else
        atts = field.atts
    end
    if isa(field.values,Array)
        if isa(field.values[1],Array)
            prec = typeof(field.values[1][1])
        else
            prec = typeof(field.values[1])
        end
        fieldvar = defVar(ds,field.name,prec,
                            tuple([f.name for f in field.dims]...),attrib=atts)
    else
        fieldvar = defVar(ds,field.name,field.values.precision,
                            tuple([f.name for f in field.dims]...),attrib=atts)
    end
end

"""
    addVar(field::NCvar,dimlist::Array{NetCDF.NcDim})

Add a variable with dimensions dimlist to a NetCDF file using `NetCDF.jl`
"""
function addVar(field::NCvar,dimlist::Array{NetCDF.NcDim})
    if ~isempty(field.units)
        atts = merge(Dict(("units" =>field.units)),field.atts)
    else
        atts = field.atts
    end
    fieldvar = NcVar(field.name,dimlist,
                        atts = atts,
                        t = field.values.precision)
end

"""
    addVar(field::NCvar})

Add a variable and its dimensions to a NetCDF file using `NetCDF.jl`
"""
function addVar(field::NCvar)
    dimlist = addDim.(field.dims)
    if ~isempty(field.units)
        attributes = merge(Dict(("units" =>field.units)),field.atts)
    else
        attributes = field.atts
    end
    fieldvar = NcVar(field.name,dimlist,
        atts = merge(Dict(("units" =>field.units)),field.atts))#,
        #t = field.values.precision)
end

"""
    addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar},var::NCvar)

Fill variable with data in netcdf file using either `NCDatasets.jl`
or `NetCDF.jl`
"""
function addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar,Array},var::NCvar,startidx=1)
    isBinData = isa(var.values,BinData)
    isNCData = isa(var.values,NCData)
    isTileData = isa(var.values,TileData)
    if isTileData
        isBinData = isa(var.values.vals,BinData)
    end
    nsteps = getnsteps(var.values)
    ndims = length(var.dims)
    if hastimedim(var)
        ndims = ndims-1
    end
    if isBinData || isNCData || isTileData || isa(var.values[1],Array) # Binary files or array of timesteps

        if isBinData && ~ isTileData
            if isa(var.values.fnames,Array)
                fnames = var.values.fnames
            else
                fnames = [var.values.fnames]
            end
        end

        for i = startidx:nsteps
            if isTileData
                writetiles.(v,Ref(var),1:var.values.numtiles,Ref(i))
            else
                if isBinData || isNCData
                    v0 = readdata(var.values,i)
                else
                    v0 = var.values[i]
                end

                if ndims == 1
                    v[:,i] = v0
                elseif ndims == 2
                    v[:,:,i] = v0
                elseif ndims == 3
                    v[:,:,:,i] = v0
                end
            end
        end

    elseif isa(var.values[1],Number) # Single array of data- just insert it
        ndims = length(size(var.values))
        if ndims == 1
            v[:] = var.values
        elseif ndims == 2
            v[:,:] = var.values
        elseif ndims == 3
            v[:,:,:] = var.values
        elseif ndims == 4
            v[:,:,:,:] = var.values
        end
    else
        print("Unrecognized values")
    end

end


"""
    writetiles(v,var,tilenum,timeidx=1)

Helper function for writing a tile to a NetCDF file.
"""
function writetiles(v,var,tilenum,timeidx=1)
    if isa(v,Array)
        v = v[findfirst(isequal(var.name),name.(v))]
    end
    tileinfo = var.values.tileinfo; tilesize = var.values.tilesize; grid = tileinfo["XC"].grid
    if isa(var.values.vals,BinData)
        iosize = var.values.vals.iosize
        prec = var.values.precision
        if isa(var.values.vals.fnames,Array)
            fnames = var.values.vals.fnames
        else
            fnames = [var.values.vals.fnames]
        end
        if length(iosize) == 3
            f = Array{Array{prec,2},2}(undef,grid.nFaces,iosize[3])
            exarray = MeshArray(grid,f)
        else
            exarray = tileinfo["XC"]
        end
        v0 = read(readbin(fnames[timeidx],prec,iosize,
                                        var.values.vals.fldidx),exarray)
    else
        if isa(var.values.vals,MeshArray) || isa(var.values.vals,MeshArrays.gcmfaces)
            v0 = var.values.vals
        else
            v0 = var.values.vals[timeidx]
        end
    end

    v0 = gettile(v0,tileinfo,tilesize,tilenum)
    numdims = length(size(v0))
    if numdims == 1
        v[:,timeidx] = v0
    elseif numdims == 2
        v[:,:,timeidx] = v0
    elseif numdims == 3
        v[:,:,:,timeidx] = v0
    end

end


"""
    addDimData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar,Array},var::NCvar)

Add dimension data to predefined dimensions in a NetCDF file.
"""
function addDimData(ds,dimvar::NCvar)
    atts = merge(Dict(("units" =>dimvar.units)),dimvar.atts)
    defVar(ds,dimvar.name,dimvar.values,(dimvar.name,),attrib=atts)
end

"""
    createfile(filename, field::Union{NCvar,Dict{String,NCvar}}, README;
                fillval=NaN, missval=NaN, itile=1, ntile=1)

Create NetCDF file and add variable + dimension definitions
using either `NCDatasets.jl` or `NetCDF.jl`
"""
function createfile(filename, field::Union{NCvar,Dict}, README;
                    fillval=NaN, missval=NaN, itile=1, ntile=1)

    if isa(field,Dict)

        dims = unique(vcat([field[v].dims for v in keys(field)]...))
        dims = filter( d -> isa(d,NCvar),dims)
        field = collect(values(field))

        fieldnames = getfield.(field[findall(hastimedim.(field))],:name)
        backend = field[1].backend
    else
        dims = field.dims
        fieldnames = field.name
        backend = field.backend
        if backend == NCDatasets
            field = [field]
        end
    end

    if isa(fieldnames,Array)
        fieldnamestring = join(fieldnames,",")
    else
        fieldnamestring = fieldnames
    end

    file_atts = vcat(["date" => Dates.format(today(),"dd-u-yyyy"),
    "Conventions" => "CF-1.6",
    "description" => fieldnamestring*" -- "*README[1]],
    [string(Char(65+(i-2))) => README[i] for i in 2:length(README)],
    [#string(Char(65+length(README)-1)) => "file created using NCTiles.jl",
    "_FillValue" => fillval,
    "missing_value" => missval,
    "itile" => itile,
    "ntile" => ntile])

    if backend == NCDatasets
        ds =  Dataset(filename,"c",attrib=file_atts)
        dimlist = addDim.(Ref(ds),dims)
        fieldvar = addVar.(Ref(ds), field)
        if length(fieldvar) == 1
            fieldvar = fieldvar[1]
        end
        return ds,fieldvar,dimlist
    elseif backend == NetCDF # Not supporting multiple variables right now
        #= for i = 1:length(field)
            f = field[1];
            if isa(f.values,BinData) || isa(f.values,NCData)
                prec = f.values.precision
            else
                prec = eltype(f.values)
            end
            nccreate(filename,f.name,Tuple(vcat([[d.name; [collect(d.values)]; d.atts] for d in f.dims]...))...,atts = f.atts,t = prec)
        end
        #ncclose(filename)
        ds = NetCDF.open(filename)
        if length(field) == 1
            fieldvar = ds[field[1].name]
        else
            fieldvar = getkey.(Ref(ds),[f.name for f in field])
        end

        return ds,fieldvar,collect(values(ds.dim)) =#
        #fieldvar = Array{NcVar,1}
        dimlist = addDim.(dims)
        #for f in field
        #    fieldvar = vcat(fieldvar,addVar(f))
        #end
        field = field
        fieldvar = addVar(field)#[addVar(f) for f in field]
        file_atts = Dict(file_atts)
        #nccreate(filename,)
        return NetCDF.create(filename,fieldvar, gatts = file_atts),fieldvar,dimlist
    end

end

"""
    readncfile(fname,backend::Module=NCDatasets)

Read in a NetCDF file and return variables/dimensions as `NCvar` structs, and
    file attributes as `Dict`. Large variables/dimensions are not loaded into
    memory. This can use either `NCDatasets.jl` or `NetCDF.jl`
"""
function readncfile(fname,backend::Module=NCDatasets)

    ds = Dataset(fname)

    dims = Dict{AbstractString,NCvar}()
    vars = Dict{AbstractString,NCvar}()

    for k in keys(ds)
        if ~haskey(dims,k)
            k_units = get(ds[k].attrib,"units","")
            k_atts = Dict(ds[k].attrib)
            if length(dimnames(ds[k])) == 1 && k == dimnames(ds[k])[1] # this variable is also a dimension
                k_dims = length(ds[k])
                if length(size(ds[k])) == 1
                    k_values = ds[k][:]
                elseif length(size(ds[k])) == 2
                    k_values = ds[k][:,:]
                else
                    k_values = NCData(fname,k,backend,typeof(ds[k].var[:][1])) # Pointer to this variable in the file
                end
                if !isa(k_values,NCData) && !any(isa.(k_values,Ref(Missing)))
                    k_values = typeof(k_values[1]).(k_values)
                end
                dims[k] = NCvar(k,k_units,k_dims,k_values,k_atts,backend)
            else # Not a dimension
                k_dims = Array{NCvar,1}()
                for d in dimnames(ds[k])
                    #global k_dims,dims
                    if !haskey(dims,d) # Add the dimension
                        d_units = get(ds[d].attrib,"units","")
                        d_dims = length(ds[d])
                        d_values = NCData(fname,d,backend,typeof(ds[d][1])) # Pointer to this variable in the file
                        d_atts = Dict(ds[d].attrib)
                        dims[d] = NCvar(d,d_units,d_dims,d_values,d_atts,backend)
                    end
                    k_dims = cat(k_dims,dims[d],dims=1)
                end
                hasTimeDim = hastimedim(k_dims)
                if  ~hasTimeDim && length(k_dims) == 1# extra variable, not main field, just load the data now
                    k_values = ds[k][:]
                elseif ~hasTimeDim && length(k_dims) == 2
                    k_values = ds[k][:,:]
                elseif ~hasTimeDim && length(k_dims) == 3
                    k_values = ds[k][:,:,:]
                elseif hasTimeDim
                    k_values = NCData(fname,k,backend,eltype(ds[k].var[:])) # Pointer to this variable in the file
                end

                if !hasTimeDim && !any(isa.(k_values,Ref(Missing)))
                    k_values = typeof(k_values[1]).(k_values)
                end

                vars[k] = NCvar(k,k_units,k_dims,k_values,k_atts,backend)
            end
        end
    end

    atts = ["_FillValue","missing_value","itile","ntile"]
    fileatts = Dict()
    for a in atts
        if haskey(ds.attrib,a)
            fileatts[a] = ds.attrib[a]
        end
    end
    close(ds)
    return vars,dims,fileatts
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

"""
    parsemeta(metafile)

Parse out an `MITgcm` metadata file and return a `Dict` of fields in the file.
"""
function parsemeta(metafile)

    meta = read(metafile,String)
    meta = split(meta,";\n")
    meta = meta[isempty.(meta).==false]
    meta = replace.(meta,Ref(",\n"=>";"))
    meta = replace.(meta,Ref("\n"=>""))
    meta = replace.(meta,Ref("}"=>"]"))
    meta = replace.(meta,Ref("{"=>"["))
    #meta = replace.(meta,Ref(" "=>""))
    meta = replace.(meta,Ref("'"=>"\""))
    meta = replace.(meta,Ref(";]"=>"]"))
    meta = replace.(meta,Ref(","=>" "))

    meta = split.(meta,"=")
    meta = [[replace(x[1]," "=>"") x[2]] for x in meta]

    metaDict = Dict{String,Any}(m[1] => m[2] for m in meta)

    for k in keys(metaDict)
        val = eval(Meta.parse(metaDict[k]))
        if isa(val[1],String)
            val = replace.(val,Ref(" "=>""))
        end
        if length(val) == 1
            val = val[1]
        end
        metaDict[k] = val
    end
    metaDict["dataprec"] = titlecase(metaDict["dataprec"])
    return metaDict

end

"""
    readAvailDiagnosticsLog(fname,fldname)

Get the information for a particular field from the `available_diagnostics.log`
    file (`MITgcm` output).
"""
function readAvailDiagnosticsLog(fname,fldname)
    availdiags = readlines(fname)
    line = availdiags[findall(occursin.(@sprintf("%-8s",fldname),availdiags))[1]]

    line = split(line,'|')
    line = lstrip.(rstrip.(line))

    diagInfo = Dict([
    "diagNum" => parse(Int,line[1]),
    "fldname" => line[2],
    "levs" => parse(Int,line[3]),
    "mate" => line[4],
    "code" => line[5],
    "units" => line[6],
    "title" => line[7]
    ])

end
