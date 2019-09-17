using NCTiles, NCDatasets, MeshArrays, Dates

function writetestfile(fname,field,package)

    if isfile(fname)
        rm(fname)
    end

    README = ["This file is written as a test for NCTiles.jl."]
    if isa(field,Dict) # TileData Example
        savenames = joinpath.(ncfiltile2d*".".*lpad.(string.(1:field["data2d"].values.numtiles),4,"0").*".nc")
        dims = unique(vcat([field[v].dims for v in keys(field)]...))
        datasets = [createfile(savenames[tidx],field,README,ntile = length(savenames), itile = tidx) for tidx in 1:length(savenames)]
        
        ds = [x[1] for x in datasets]
        fldvars = [x[2] for x in datasets]
        #dims = [x[3] for x in datasets]

        for k in keys(field)
            if isa(field[k].values,TileData)
                addData(fldvars,field[k])
            else
                tmpfldvars = [fv[findfirst(isequal(k),name.(fv))] for fv in fldvars]
                addData.(tmpfldvars,Ref(field[k]))
            end
        end

        for dim in dims
            addDimData.(ds,Ref(dim))
        end
        close.(ds)

    else
        ds,fldvar,dimlist = createfile(fname,field,README)

        addData(fldvar,field)
        if package == NCDatasets
            addDimData.(Ref(ds),field.dims)
        end
        # Close the file
        close(ds)
    end
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
        if isa(fildim[1],DateTime)
            dimvals = [Hour.(ds["tim"][t] .- DateTime(replace(d.units,"days since "=>""),"yyyy-mm-dd HH:MM:SS")).value for t in 1:length(ds["tim"])]/24
            dimcheck["values"] = d.values == dimvals
        else
            dimcheck["values"] = d.values == fildim
            if length(size(d.values)) > 1
                dimcheck["values"] = d.values[:,1] == fildim
            end
        end
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
        varcheck["values"] = true
        if isa(checkfld.values[1],Array)
            if NCTiles.hastimedim(checkfld)
                for t in 1:length(checkfld.values)
                    if length(size(dsvar)) == 3
                        varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,:,t]
                    elseif length(size(dsvar)) == 4
                        varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,:,:,t]
                    elseif length(size(dsvar)) == 2
                        varcheck["values"] = varcheck["values"] && checkfld.values[t] == dsvar[:,t]
                    end
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
        itile = Int(ds.attrib["itile"])
        varcheck["values"] = true
        if isa(checkfld.values.vals,BinData)
            for t in 1:length(checkfld.values.vals.fnames)
                data = read(readbin(checkfld.values.vals.fnames[t],
                                    checkfld.values.vals.precision,checkfld.values.vals.iosize,
                                    checkfld.values.vals.fldidx),checkfld.values.tileinfo["XC"])
                testdata = NCTiles.gettile(data,checkfld.values.tileinfo,checkfld.values.tilesize,itile)
                if length(size(dsvar)) == 3
                    varcheck["values"] = varcheck["values"] && testdata == dsvar[:,:,t]
                elseif length(size(dsvar)) == 4
                    varcheck["values"] = varcheck["values"] && testdata == dsvar[:,:,:,t]
                elseif length(size(dsvar)) == 2
                    varcheck["values"] = varcheck["values"] && testdata == dsvar[:,t]
                end
            end
        else
            testdata = NCTiles.gettile(checkfld.values.vals,checkfld.values.tileinfo,checkfld.values.tilesize,itile)
            varcheck["values"] = varcheck["values"] && ((testdata .== dsvar) + isnan.(dsvar)) == ones(Int,size(dsvar)...)
        end
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
    fnamestile2d = [mktemp()[1] for t in 1:length(data2d)]

    for t in 1:length(data2d)
        write(fnames2d[t],hton.(data2d[t]))
    end
    
    for t in 1:length(data3d)
        write(fnames3d[t],hton.(data3d[t]))
    end

    tilesize = (16,16)
    grid = GridSpec("CS32",getgrid()*"/")
    for f in fnamestile2d
        cp(joinpath(grid.path,"Depth.data"),f,force=true)
    end
    dep = -1 .* collect([-5 -15 -25 -35 -45 -55 -65 -75.0049972534180 -85.0250015258789 -95.0950012207031 -105.309997558594 -115.870002746582 -127.150001525879 -139.740005493164 -154.470001220703 -172.399993896484 -194.735000610352 -222.710006713867 -257.470001220703 -299.929992675781 -350.679992675781 -409.929992675781 -477.470001220703 -552.710021972656 -634.734985351563 -722.400024414063 -814.469970703125 -909.739990234375 -1007.15502929688 -1105.90502929688 -1205.53503417969 -1306.20495605469 -1409.15002441406  -1517.09497070313 -1634.17504882813 -1765.13500976563 -1914.15002441406 -2084.03491210938 -2276.22509765625 -2491.25000000000 -2729.25000000000 -2990.25000000000 -3274.25000000000 -3581.25000000000 -3911.25000000000 -4264.25000000000 -4640.25000000000 -5039.25000000000 -5461.25000000000 -5906.25000000000]')
    tile_ex = Dict(["tilesize" => tilesize,
                    "grid" => grid,
                    "dims" =>   [NCvar("i_c","1",tilesize[1],1:tilesize[1],Dict("long_name" => "Cartesian coordinate 1"),NCDatasets),
                                NCvar("j_c","1",tilesize[2],1:tilesize[2],Dict("long_name" => "Cartesian coordinate 2"),NCDatasets),
                                NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
                                NCvar("tim","days since 1992-1-1 0:0:0",Inf,time_steps,Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)],
                    "fnamestile2d" => fnamestile2d
                    ])

    return Dict(["data2d" => data2d, "data3d" => data3d,
                 "units" => units, "dims" => dims,
                 "longname2d" => longname2d,"longname3d" => longname3d,
                 "fnames2d" => fnames2d, "fnames3d" => fnames3d, "tile_ex" => tile_ex])

end

function getgrid()

    testdir = abspath(joinpath(dirname(pathof(NCTiles)),"..","test"))
    griddir = joinpath(testdir,"GRID_CS32")
    if ~ispath(griddir)
        try
            file = joinpath(testdir,"GRID_CS32.tar.gz")
            download("https://github.com/gaelforget/GRID_CS32/archive/master.tar.gz",file)
            run(`tar -xzf $file -C $testdir`)
            olddir = joinpath(testdir,"GRID_CS32-master")
            run(`mv $olddir $griddir`)
            run(`rm $file`)
        catch e
            println("Could not download and extract grid. To fix, get the grid from https://github.com/gaelforget/GRID_CS32/ and put it in the "*testdir*" directory.")
        end
    end
        
    return testdir
end