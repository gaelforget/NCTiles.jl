using NCTiles, NCDatasets

function writetestfile(fname,field,package)

    if isfile(fname)
        rm(fname)
    end

    README = ["This file is written as a test for NCTiles.jl."]
    ds,fldvar,dimlist = createfile(fname,field,README)

    addData(fldvar,field)
    if package == NCDatasets
        addDimData.(Ref(ds),field.dims)
    end
    # Close the file
    close(ds)
end

function testfile(fname,checkfld)
    pass = true
    ds = Dataset(fname)
    varcheck = Dict()
    varcheck["dims"] = Dict()
    for d in checkfld.dims
        fildim = ds[d.name]
        dimcheck = Dict()
        # check size
        dimcheck["size"] = length(fildim) == length(d.values)
        dimcheck["values"] = d.values == fildim
        pass = pass && dimcheck["size"] && dimcheck["values"]
        # check attributes
        for k in keys(fildim.attrib)
            if k == "units"
                dimcheck["units"] = fildim.attrib[k] == d.units
                pass = pass && dimcheck["units"]
            else
                dimcheck[k] = fildim.attrib[k] == d.atts[k]
                pass = pass && dimcheck[k]
            end
        end

        varcheck["dims"][d.name] = dimcheck
    end

    dsvar = ds[checkfld.name]

    for k in keys(dsvar.attrib)
        if k == "units"
            varcheck["units"] = dsvar.attrib[k] == checkfld.units
            pass = pass && varcheck["units"]
        else
            varcheck[k] = dsvar.attrib[k] == checkfld.atts[k]
            pass = pass && varcheck[k]
        end
    end

    # Check values of test data
    varcheck["values"] = false
    if isa(checkfld.values,Array)
        if isa(checkfld.values[1],Array)
            varcheck["values"] = true
            for t in 1:length(checkfld.values)
                if length(size(dsvar)) == 3
                    varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,:,t]
                elseif length(size(dsvar)) == 4
                    varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,:,:,t]
                elseif length(size(dsvar)) == 2
                    varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,t]
                end
            end
        else
            varcheck["values"] = varcheck["values"] && checkfld.values == dsvar
        end
    elseif isa(checkfld.values,BinData)
        varcheck["values"] = true
        for t in 1:length(checkfld.values.fnames)
            testdata = readbin(checkfld.values,t)
            if length(size(dsvar)) == 3
                varcheck["values"] = varcheck["values"] && testdata == dsvar[:,:,t]
            elseif length(size(dsvar)) == 4
                varcheck["values"] = varcheck["values"] && testdata == dsvar[:,:,:,t]
            elseif length(size(dsvar)) == 2
                varcheck["values"] = varcheck["values"] && testdata == dsvar[:,t]
            end
        end
    elseif isa(checkfld.values,NCData)
        varcheck["values"] = true
        testfld = Dataset(checkfld.values.fname)[checkfld.values.varname]
        for t in 1:size(testfld)[end]
            if length(size(dsvar)) == 3
                varcheck["values"] = varcheck["values"] && testfld[:,:,t] == dsvar[:,:,t]
            elseif length(size(dsvar)) == 4
                varcheck["values"] = varcheck["values"] && testfld[:,:,:,t] == dsvar[:,:,:,t]
            elseif length(size(dsvar)) == 2
                varcheck["values"] = varcheck["values"] && testfld[:,t] == dsvar[:,t]
            end
        end
    elseif isa(checkfld.values,TileData)
        println("No value check for TileData")
    end
    pass = pass && varcheck["values"]

    if ~pass
        missmatchstring = ""
        for k in keys(varcheck)
            if k == "dims"
                for d in keys(varcheck[k])
                    for a in keys(varcheck[k][d])
                        if ~varcheck[k][d][a]
                            missmatchstring = join([missmatchstring,join([k,d,a],'-')],", ")
                        end
                    end
                end
            else
                if ~varcheck[k]
                    missmatchstring = join([missmatchstring,k],", ")
                end
            end
        end
        println("The following file elements do not match the field elements:\n"*missmatchstring[3:end])
    end
    return pass
end

function maketestdata()
    lon=-180:20:180; lat=-90:20:90;
    depth = 5:10:100
    n1,n2,n3 = (length(lon),length(lat),length(depth))
    
    timeinterval = 1
    nsteps = 5
    time_steps = timeinterval/2:timeinterval:timeinterval*nsteps
    time_units = "days from 1992-01-01"
    
    data2d = [ones(Float32,n1,n2) for t in 1:nsteps]
    data3d = [ones(Float32,n1,n2,n3) for t in 1:nsteps]
    units = "test_units"
    longname2d = "Two Dimensional Test Data"
    longname3d = "Three Dimensional Test Data"
    
    dims = [NCvar("lon_c","degrees_east",n1,lon,Dict("long_name" => "longitude"),NCDatasets),
            NCvar("lat_c","degrees_north",n2,lat,Dict("long_name" => "latitude"),NCDatasets),
            NCvar("dep_c","m",n3,depth,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
            NCvar("time","days from 1992-01-01",Inf,collect(time_steps),Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)
            ]


    #fnames2d = ["testdata/data2d."*lpad(t,2,"0")*".data" for t in 1:length(data2d)]
    #fnames3d = ["testdata/data3d."*lpad(t,2,"0")*".data" for t in 1:length(data3d)]

    fnames2d = [mktemp()[1] for t in 1:length(data2d)]
    fnames3d = [mktemp()[1] for t in 1:length(data3d)]

    for t in 1:length(data2d)
        write(fnames2d[t],hton.(data2d[t]))
    end
    
    for t in 1:length(data3d)
        write(fnames3d[t],hton.(data3d[t]))
    end

    return Dict(["data2d" => data2d, "data3d" => data3d,
                 "units" => units, "dims" => dims,
                 "longname2d" => longname2d,"longname3d" => longname3d,
                 "fnames2d" => fnames2d, "fnames3d" => fnames3d])

end