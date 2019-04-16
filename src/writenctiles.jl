using NCDatasets,Dates

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

struct NCvar
    name::String
    longname::String
    units::String
    dims
    values
    atts::Union{Dict,Nothing}
end

struct Bindata # Pointer to data stored in binary files- contains info needed to read in
    fnames::Union{Array{String},String}
    precision::Type
    iosize::Tuple
end

function addDim(ds,dimvar::NCvar)
    defDim(ds,dimvar.name,dimvar.dims[1])
end

function addData(v::NCDatasets.CFVariable,var::NCvar)

    if isa(var.values,Bindata) ||  isa(var.values[1],Array) # Binary files or array of timesteps

        isBinData = isa(var.values,Bindata)

        if isBinData
            ndims = length(var.values.iosize)
            nsteps = length(var.values.fnames)
        else
            ndims = length(size(var.values[1]))
            nsteps = length(var.values)
        end

        for i = 1:nsteps
            if isBinData
                v0 = readbin(var.values.fnames[i],var.values.precision,var.values.iosize)
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
    atts = merge(Dict(("long_name" => dimvar.longname,"units" =>dimvar.units)),dimvar.atts)
    defVar(ds,dimvar.name,dimvar.values,(dimvar.name,),attrib=atts)
end

function createfile(filename, fldname, README; fillval=NaN, missval=NaN, ff=1, ntile=1)

    file_atts = vcat(["date" => Dates.format(today(),"dd-u-yyyy"),
                "Conventions" => "CF-1.6",
                "description" => fldname*" -- "*README[1]],
                [string(Char(65+(i-2))) => README[i] for i in 2:length(README)],
                [string(Char(65+length(README)-1)) => "file created using NCTiles.jl",
                "_FillValue" => fillval,
                "missing_value" => missval,
                "itile" => ff,
                "ntile" => ntile])
                
    return Dataset("test.nc","c",attrib=file_atts)

end