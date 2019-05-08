using NCDatasets,NetCDF,Dates

struct NCvar
    name::String
    units::String
    dims
    values
    atts::Union{Dict,Nothing}
    backend::Module
end

struct Bindata # Pointer to data stored in binary files- contains info needed to read in
    fnames::Union{Array{String},String}
    precision::Type
    iosize::Tuple
end

struct NCData
    fname::AbstractString
    varname::AbstractString
    backend::Module
    precision::Type
end

# Functions to get single time step
function getindex(v::NCvar,i::Integer)
    NCvar(v.name,v.units,v.dims,v.values[i],v.atts,v.backend)
end

function getindex(d::Bindata,i) # Gets file at time index i
    if isa(d.fnames,Array)
        newfnames = d.fnames[i]
    else
        newfnames = d.fnames
    end

    Bindata(newfnames,d.precision,d.iosize)
end

function getindex(d::NCData,i) # Gets data at time index i
    readncdata(d,i)
end

function readbin(fname::String,prec::Type,iosize::Tuple)
    n1,n2 = iosize

    if prec == Float64
        reclen = 6
    else
        reclen=4
    end

    n3=Int64(stat(fname).size/n1/n2/reclen)

    fid = open(fname)
    field = Array{prec,1}(undef,(n1*n2*n3))
    read!(fid,field)
    field = hton.(field)

    n3>1 ? s=(n1,n2,n3) : s=(n1,n2)
    return reshape(field,s)

end

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

function addDim(ds::NCDatasets.Dataset,dimvar::NCvar) # NCDatasets
    defDim(ds,dimvar.name,dimvar.dims[1])
end

function addDim(dimvar::NCvar) #NetCDF
    NcDim(dimvar.name,collect(dimvar.values),
        atts = merge(Dict(("units" =>dimvar.units)),dimvar.atts),
        unlimited = dimvar.dims==Inf)
end

function addVar(ds::NCDatasets.Dataset,field::NCvar)
    atts = merge(Dict(("units" =>field.units)),field.atts)
    if isa(field.values,Array)
        fieldvar = defVar(ds,field.name,field.values,tuple([f.name for f in field.dims]...),attrib=atts)
    else
        fieldvar = defVar(ds,field.name,field.values.precision,tuple([f.name for f in field.dims]...),attrib=atts)
    end
end

function addVar(field::NCvar,dimlist::Array{NetCDF.NcDim})
    fieldvar = NcVar(field.name,dimlist,
        atts = merge(Dict(("units" =>field.units)),field.atts), 
        t = field.values.precision)
end

function addVar(field::NCvar)
    dimlist = addDim.(field.dims)
    fieldvar = NcVar(field.name,dimlist,
        atts = merge(Dict(("units" =>field.units)),field.atts))#, 
        #t = field.values.precision)
end

function addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar},var::NCvar)
    isBinData = isa(var.values,Bindata)
    isNCData = isa(var.values,NCData)

    if isBinData || isNCData ||  isa(var.values[1],Array) # Binary files or array of timesteps

        if isBinData
            ndims = length(var.values.iosize)
            if isa(var.values.fnames,Array)
                nsteps = length(var.values.fnames)
                fnames = var.values.fnames
            else
                nsteps = 1
                fnames = [var.values.fnames]
            end
        elseif isNCData
            ndims = length(var.dims)
            if hastimedim(var)
                ndims = ndims-1
                timedimidx = findtimedim(var)
                nsteps = var.dims[timedimidx].dims[1]
            else
                nsteps = 1
            end
            
        else
            ndims = length(size(var.values[1]))
            nsteps = length(var.values)
        end
        
        for i = 1:nsteps
            if isBinData
                v0 = readbin(fnames[i],var.values.precision,var.values.iosize)
            elseif isNCData
                v0 = readncdata(var.values,i)
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

function addDimData(ds,dimvar::NCvar)
    atts = merge(Dict(("units" =>dimvar.units)),dimvar.atts)
    defVar(ds,dimvar.name,dimvar.values,(dimvar.name,),attrib=atts)
end

function createfile(filename, field::Union{NCvar,Dict{AbstractString,NCvar}}, README; fillval=NaN, missval=NaN, ff=1, ntile=1)
    
    if isa(field,Dict{AbstractString,NCvar})
        
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
    
    file_atts = vcat(["date" => Dates.format(today(),"dd-u-yyyy"),
                "Conventions" => "CF-1.6",
                "description" => join(fieldnames,",")*" -- "*README[1]],
                [string(Char(65+(i-2))) => README[i] for i in 2:length(README)],
                [string(Char(65+length(README)-1)) => "file created using NCTiles.jl",
                "_FillValue" => fillval,
                "missing_value" => missval,
                "itile" => ff,
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
            if isa(f.values,Bindata) || isa(f.values,NCData)
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

# Helpers for finding time dimension
function istimedim(d::NCvar)

    timeUnits = ["minutes","seconds","hours","days","minute","second","hour","day"]
    return any(occursin.(timeUnits,Ref(lowercase(d.units))))

end

function findtimedim(v::NCvar)
    return findall(istimedim.(v.dims))[1]
end

function hastimedim(v::NCvar)
    ncvardim = isa.(v.dims,NCvar)
    return any(ncvardim) && any(istimedim.(v.dims[ncvardim]))
end

function hastimedim(dims::Array{NCvar})
    return any(istimedim.(dims))
end

function parsemeta(metafile)

    meta = read(metafile,String)
    meta = split(meta,";\n")
    meta = meta[isempty.(meta).==false]
    meta = replace.(meta,Ref(",\n"=>";"))
    meta = replace.(meta,Ref("\n"=>""))
    meta = replace.(meta,Ref("}"=>"]"))
    meta = replace.(meta,Ref("{"=>"["))
    meta = replace.(meta,Ref(" "=>""))
    meta = replace.(meta,Ref("'"=>"\""))
    meta = replace.(meta,Ref(";]"=>"]"))
    meta = replace.(meta,Ref(","=>" "))

    meta = split.(meta,"=")

    metaDict = Dict{String,Any}(m[1] => m[2] for m in meta)

    for k in keys(metaDict)
        val = eval(Meta.parse(metaDict[k]))
        if length(val) == 1
            val = val[1]
        end
        metaDict[k] = val
    end
    metaDict["dataprec"] = titlecase(metaDict["dataprec"])
    return metaDict

end

using Printf
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
