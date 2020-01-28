
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
argument `i` can be used to read a specific records / times. Not exported.
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

Generic wrapper function to read data from a file regardless of file/data type. Not exported.
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
