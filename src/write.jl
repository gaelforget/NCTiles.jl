
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
    replacevalues(vals,ncvar::NCvar)

Helper function: Replaces the "values" attribute of an ncvar without changing
anything else. Used to apply land mask.
"""
replacevalues(vals,ncvar::NCvar) = NCvar(ncvar.name,ncvar.units,ncvar.dims,vals,ncvar.atts,ncvar.backend)

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
    if isa(vals,BinData)
        prec = vals.precision
    else
        prec = eltype(vals)
    end
    τ=Tiles(grid,ni,nj); TileMap=MeshArray(grid,prec);
    [TileMap[τ[ii]["face"]][τ[ii]["i"],τ[ii]["j"]].=ii for ii in 1:length(τ)]
    tileinfo = Dict("tileNo" => TileMap)
    gridvars=GridLoad(grid)
    tileinfo["XC"]=gridvars["XC"]
    tileinfo["YC"]=gridvars["YC"]
    return TileData(vals,tileinfo,tilesize,prec,Int(maximum(tileinfo["tileNo"])))
end
TileData(vals,tilesize::Tuple,grid::String="LatLonCap") = TileData(vals,tilesize::Tuple,GridSpec(grid))
"""
    replacevalues(vals,td::TileData)

Helper function: Replaces the "values" attribute of a TilData struct without changing
anything else. Used to apply land mask.
"""
replacevalues(vals,td::TileData) = TileData(vals,td.tileinfo,td.tilesize,td.precision,td.numtiles)

import Base: getindex

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
    dimvals = dimvar.values[:]
    NcDim(dimvar.name,collect(dimvals),
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
    if isa(field.values,Array)
        prec = eltype(field.values)
    else
        prec = field.values.precision
    end
    fieldvar = NcVar(field.name,dimlist,
                        atts = atts,
                        t = prec)
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
        atts = merge(Dict(("units" =>field.units)),field.atts))
end

"""
    addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar},var::NCvar)

Fill variable with data in netcdf file using either `NCDatasets.jl`
or `NetCDF.jl`
"""
function addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar,Array},var::NCvar;startidx=1,land_mask=nothing)
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

    if isTileData
        if isnothing(land_mask)
            gridvars = GridLoad(var.values.tileinfo["XC"].grid)
            if any(occursin.("_s",getfield.(var.dims,:name)))
                land_mask = gridvars["hFacS"]
            elseif any(occursin.("_w",getfield.(var.dims,:name)))
                land_mask = gridvars["hFacW"]
            else
                land_mask = gridvars["hFacC"]
            end
            for f in land_mask.fIndex
                for d in 1:size(land_mask,2)
                    land_mask[f,d][land_mask[f,d].==0] .= NaN
                    land_mask[f,d][land_mask[f,d].>0] .= 1
                end
            end
        end
    end

    if (isa(var.values,Array) && isa(var.values[1],Number)) || ~hastimedim(var)
        if isa(var.values,Array)
            v0 = var.values
        else
            v0 = readdata(var.values,:)
        end

        
        if ~isnothing(land_mask)
            if length(size(v0)) == 1 && length(size(land_mask)) == 2
                v0 = v0 .* land_mask[:,1]
            elseif length(size(v0)) == 2 && length(size(land_mask)) == 3
                v0 = v0 .* land_mask[:,:,1]
            else
                v0 = v0 .* land_mask
            end
        end
        if isTileData
            tmpvar = replacevalues(replacevalues(v0,var.values),var)
            writetiles.(v,Ref(tmpvar),1:var.values.numtiles,Ref(:),Ref(land_mask))
        else
            ndims = length(var.dims)
            checkdims(v0,var::NCvar)
            if ndims == 1
                v[:] = v0
            elseif ndims == 2
                v[:,:] = v0
            elseif ndims == 3
                v[:,:,:] = v0
            elseif ndims == 4
                v[:,:,:,:] = v0
            end
        end

    elseif isBinData || isNCData || isTileData || isa(var.values[1],Array) # Binary files or array of timesteps

        if isBinData && ~ isTileData
            if isa(var.values.fnames,Array)
                fnames = var.values.fnames
            else
                fnames = [var.values.fnames]
            end
        end

        if ~hastimedim(var) #if no time steps, just read in and insert the whole thing
            v0 = readdata(var.values,:)
            v[:] = v0
        else
            for i = startidx:nsteps
                if isBinData || isNCData || isTileData
                    v0 = readdata(var.values,i)
                else
                    v0 = var.values[i]
                end

                if ~isnothing(land_mask)
                    if length(size(v0)) == 1 && length(size(land_mask)) == 2
                        v0 = v0 .* land_mask[:,1]
                    elseif length(size(v0)) == 2 && length(size(land_mask)) == 3
                        v0 = v0 .* land_mask[:,:,1]
                    else
                        v0 = v0 .* land_mask
                    end
                end
                if isTileData
                    tmpvar = replacevalues(replacevalues(v0,var.values),var)
                    writetiles.(v,Ref(tmpvar),1:var.values.numtiles,Ref(i),Ref(land_mask))
                else
                    checkdims(v0,var)
                    if ndims == 0
                        v[i] = v0
                    elseif ndims == 1
                        v[:,i] = v0
                    elseif ndims == 2
                        v[:,:,i] = v0
                    elseif ndims == 3
                        v[:,:,:,i] = v0
                    end
                end
            end
        end

    else
        print("Unrecognized values")
    end
    return nothing
end


"""
    writetiles(v,var,tilenum,timeidx=1)

Helper function for writing a tile to a NetCDF file.
"""
function writetiles(v,var,tilenum,timeidx=1,land_mask=nothing)
    
    if isa(v,Array)
        if var.backend == NCDatasets
            v = v[findfirst(isequal(var.name),name.(v))]
        else
            varnames = [vtmp.name for vtmp in v]
            v = v[findfirst(isequal(var.name),varnames)]
        end
    end
    tileinfo = var.values.tileinfo; tilesize = var.values.tilesize
    if isa(var.values.vals,MeshArray) || isa(var.values.vals,MeshArrays.gcmfaces)
        v0 = var.values.vals
    else
        v0 = readdata(var.values,timeidx)
    end

    v0 = gettile(v0,tileinfo,tilesize,tilenum)
    checkdims(v0,var::NCvar)
    if hastimedim(var)
        numdims = ndims(v0)
        if numdims == 1
            v[:,timeidx] = v0
        elseif numdims == 2
            v[:,:,timeidx] = v0
        elseif numdims == 3
            v[:,:,:,timeidx] = v0
        end
    else
        v[:] = v0
    end

end

"""
    addDimData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar,Array},var::NCvar)

Add dimension data to predefined dimensions in a NetCDF file.
"""
function addDimData(ds,dimvar::NCvar)
    atts = merge(Dict(("units" =>dimvar.units)),dimvar.atts)
    defVar(ds,dimvar.name,dimvar.values[:],(dimvar.name,),attrib=atts)
end

"""
    createfile(filename, field::Union{NCvar,Dict{String,NCvar}}, README;
                fillval=NaN, missval=NaN, itile=1, ntile=1)

Create NetCDF file and add variable + dimension definitions
using either `NCDatasets.jl` or `NetCDF.jl`
"""
function createfile(filename, field::Union{NCvar,Dict}, rdm="";
                    fillval=NaN, missval=NaN, itile=1, ntile=1, attribs=nothing)

    if isa(field,Dict)

        dims = getdims(field::Dict)
        #dims = unique(vcat([field[v].dims for v in keys(field)]...))
        #dims = filter( d -> isa(d,NCvar),dims)
        field = collect(values(field))

        fldnames = getfield.(field[findall(hastimedim.(field))],:name)
        backend = field[1].backend
    else
        dims = field.dims
        fldnames = field.name
        backend = field.backend
        if backend == NCDatasets
            field = [field]
        end
    end

    if isa(fldnames,Array)
        fldnamestring = join(filter(x -> x !== "climatology_bounds", fldnames),",")
    else
        fldnamestring = fldnames
    end

    README=deepcopy(rdm)
    if ~isempty(README) && ~occursin(fldnamestring*" -- ",README[1])
        README[1] = fldnamestring*" -- "*README[1]
    end

    if ~isnothing(attribs)
        fillval = pop!(attribs,"_FillValue",NaN)
        missingval = pop!(attribs,"missing_value",NaN)
        haskey(attribs,"itile") ? itile = pop!(attribs,"itile") : nothing
        haskey(attribs,"ntile") ? ntile = pop!(attribs,"ntile") : nothing
        description = [pop!(attribs,k,"") for k in ["description","A","B","C","D","E","F","G","H","I","J"]]
        description = description[.~isempty.(description)]
        if isempty(README)
            README = description
        elseif !(description == README)
            README = vcat(README,description)
        end
        pop!(attribs,"date",nothing)
        if isempty(attribs); attribs=nothing; end
    end

    if ~isempty(README)
        if isa(README,Array) && length(README) > 1
            description = vcat("description" => README[1],
            [string(Char(65+(i-2))) => README[i] for i in 2:length(README)])
        else
            description = ["description" => isa(README,Array) ? README[1] : README]
        end
    else
        description = ""
    end

    file_atts = vcat(["date" => Dates.format(today(),"dd-u-yyyy"),
    "Conventions" => "CF-1.6"], "NCTiles" => NCTiles.version())

    if ~isempty(description)
        file_atts = vcat(file_atts,description)
    end
    file_atts = vcat(file_atts,
    [#string(Char(65+length(README)-1)) => "file created using NCTiles.jl",
    "_FillValue" => fillval,
    "missing_value" => missval,
    "itile" => itile,
    "ntile" => ntile])

    if ~isnothing(attribs)
        file_atts = vcat(file_atts,[k => attribs[k] for k in keys(attribs)])
    end

    if backend == NCDatasets
        ds =  Dataset(filename,"c",attrib=file_atts)
        dimlist = addDim.(Ref(ds),dims)
        fieldvar = addVar.(Ref(ds), field)
        if length(fieldvar) == 1
            fieldvar = fieldvar[1]
        end
        return ds,fieldvar,dimlist
    elseif backend == NetCDF
        dimlist = addDim.(dims)
        if isa(field,Array)
            dimnames = getfield.(dims,:name)
            field = filter(f ->~any(f.name.==dimnames),field) # need to get dimlist for each variable
            nflds = length(field)
            fieldvar = Array{NcVar}(undef,nflds)
            for i in 1:nflds
                f = field[i]
                dimidx = [findfirst(d -> d .== dim.name,dimnames) for dim in f.dims]
                fieldvar[i] = addVar(f,dimlist[dimidx])
            end
        else
            fieldvar = addVar(field,dimlist)
        end
        file_atts = Dict(file_atts)
        return NetCDF.create(filename,fieldvar, gatts = file_atts),fieldvar,dimlist
    end

end

import Base: write
using Pkg
"""
    write(myfld::NCvar,savename::String;README="",globalattribs=Dict())

Creates NetCDF file and writes myfld and all its dimensions to the file.
"""
function write(myfld::NCvar,savename::String;README="",globalattribs=Dict())
    if hastiledata(myfld) # Create one file for each tile

        numtiles = myfld.values.numtiles
        savenames = savename*".".*lpad.(string.(1:numtiles),4,"0").*".nc"

        datasets = [createfile(savenames[tidx],myfld,README, itile = tidx, ntile = length(savenames), attribs = globalattribs) for tidx in 1:length(savenames)]

        ds = [x[1] for x in datasets]
        fldvars = [x[2] for x in datasets]

        addData(fldvars,myfld)

        if myfld.backend == NCDatasets
            dims = myfld.dims

            for dim in dims
                addDimData.(ds,Ref(dim))
            end
            
            close.(ds)
        else
            if Pkg.installed()["NetCDF"] < v"0.10"
                NetCDF.close.(ds)
            else
                finalize.(ds)
            end
        end
    else # Create one file
        ds,fldvar,dimlist = createfile(savename,myfld,README,attribs=globalattribs)

        # Add field and dimension data
        addData(fldvar,myfld)

        if myfld.backend == NCDatasets
            addDimData.(Ref(ds),myfld.dims)

            # Close the file
            close(ds)
        else
            if Pkg.installed()["NetCDF"] < v"0.10"
                NetCDF.close(ds)
            else
                finalize(ds)
            end
        end
    end
end

"""
    write(myflds::Dict,savename::String;README="",globalattribs=Dict())

Creates NetCDF file and writes myflds and all their dimensions to the file.
"""
function write(myflds::Dict{AbstractString,NCvar},savename::String;README="",globalattribs=Dict())

    if hastiledata(myflds)
        fldnames = collect(keys(myflds))
        tilefld = myflds[fldnames[findfirst([hastiledata(myflds[f]) for f in fldnames])]]
        numtiles = tilefld.values.numtiles

        landidx = findfirst(get.([myflds[f].atts for f in fldnames],"standard_name","none").=="land_binary_mask")
        if ~isnothing(landidx); land_mask = myflds[fldnames[landidx]].values; else land_mask = nothing; end
        if isa(land_mask,TileData); land_mask = land_mask.vals; end
        savenames = savename*".".*lpad.(string.(1:numtiles),4,"0").*".nc"

        datasets = [createfile(savenames[tidx],myflds,README, itile = tidx, ntile = length(savenames), attribs = globalattribs) for tidx in 1:length(savenames)]

        ds = [x[1] for x in datasets]
        fldvars = [x[2] for x in datasets]

        for k in keys(myflds)
            if isa(myflds[k].values,TileData)
                addData(fldvars,myflds[k],land_mask = land_mask)
            else
                if myflds[k].backend == NCDatasets
                    tmpfldvars = [fv[findfirst(isequal(k),name.(fv))] for fv in fldvars]
                    #tmpnames = [name.(fv) for fv in fldvars]
                else
                    tmpfldvars = [fv[findfirst(isequal(k),getproperty.(fv,:name))] for fv in fldvars]
                    #tmpnames = [getproperty.(fv,:name) for fv in fldvars]
                end
                
                addData.(tmpfldvars,Ref(myflds[k]))
            end
        end
        if isa(ds[1],NCDatasets.Dataset)
            dims = unique(vcat([myflds[v].dims for v in keys(myflds)]...))
            dims = filter( d -> isa(d,NCvar),dims)

            for dim in dims
                addDimData.(ds,Ref(dim))
            end

            close.(ds)
        else
            if Pkg.installed()["NetCDF"] < v"0.10"
                NetCDF.close.(ds)
            else
                finalize.(ds)
            end
        end

    else
        dims = getdims(myflds)

        ## Here down to move to other branch- high level API for writing
        ds,fldvar,dimlist = createfile(savename,myflds,README,attribs=globalattribs)

        # Add field data
        for k in keys(myflds)
            addData(ds[k],myflds[k])
        end

        if isa(ds,NCDatasets.Dataset)
            # Only insert data for dims with data
            dims = filter( d -> ~isempty(d.values[:]),dims)

            # Add dimension data
            addDimData.(Ref(ds),dims)

            # Close the file
            close(ds)
        else
            if Pkg.installed()["NetCDF"] < v"0.10"
                NetCDF.close(ds)
            else
                finalize(ds)
            end
        end
    end
end

write(myflds::Dict,savename::String;README="",globalattribs=Dict()) = write(Dict{AbstractString,NCvar}(myflds),savename;README="",globalattribs=Dict())