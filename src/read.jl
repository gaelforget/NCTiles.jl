
"""
    readbin(fname::String,prec::Type,iosize::Tuple,fldidx=1)

Read in a binary file to an Array.
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

Read in a binary file as an array as specified by BinData
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

"""
    readdata(flddata,tidx=1)

Generic wrapper function to read data from a file regardless of file/data type.
"""
function readdata(flddata,tidx=1)
    if isa(flddata,BinData)
        res = readbin(flddata,tidx)
    elseif isa(flddata,NCData)
        res = readncdata(flddata,tidx)
    elseif isa(flddata,TileData)
        if isa(flddata.vals,BinData)
            iosize = flddata.vals.iosize
            grid = flddata.tileinfo["XC"].grid
            if length(iosize) == 3
                f = Array{Array{flddata.vals.precision,2},2}(undef,grid.nFaces,iosize[3])
                exarray = MeshArray(grid,f)
            else
                exarray = flddata.tileinfo["XC"]
            end
            res = read(readbin(flddata.vals,tidx),exarray)
        else
            if isa(flddata.vals,MeshArray) || isa(flddata.vals,MeshArrays.gcmfaces)
                res = flddata.vals
            else
                res = flddata.vals[tidx]
            end
        end
    end
    return res
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
                    if !haskey(dims,d) 
                        if haskey(ds,d) # Add the dimension
                            d_units = get(ds[d].attrib,"units","")
                            d_dims = length(ds[d])
                            d_values = NCData(fname,d,backend,typeof(ds[d][1])) # Pointer to this variable in the file
                            d_atts = Dict(ds[d].attrib)
                            dims[d] = NCvar(d,d_units,d_dims,d_values,d_atts,backend)
                        else
                            dims[d] = NCvar(d,"",ds.dim[d],[],Dict(),backend)
                        end

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
