using NCTiles,NCDatasets,NetCDF,Dates,Printf

interpdir = "data/cs510/diags_interp/"
availdiagsfname = "data/available_diagnostics.log"
ncdir = "data/cs510/interp_nctiles/"
fnameprefix = ncdir*"/darwin_v0.2_cs510"
groups = ["Nutrients","Bulk_Ecosystem_Characteristcs","Phytoplankton_Functional_Types","Ocean_Color"]
README = readlines("data/cs510/README")

if !isdir(ncdir)
    mkpath(ncdir)
end

# Hard coding dims for now
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
nsteps = 2
timeinterval = 3
timesteps = timeinterval/2:timeinterval:timeinterval*nsteps
timeunits = "days from 1992-01-01"
timestart = DateTime(1992,01,01)
dep = -1 .* collect([-5 -15 -25 -35 -45 -55 -65 -75.0049972534180 -85.0250015258789 -95.0950012207031 -105.309997558594 -115.870002746582 -127.150001525879 -139.740005493164 -154.470001220703 -172.399993896484 -194.735000610352 -222.710006713867 -257.470001220703 -299.929992675781 -350.679992675781 -409.929992675781 -477.470001220703 -552.710021972656 -634.734985351563 -722.400024414063 -814.469970703125 -909.739990234375 -1007.15502929688 -1105.90502929688 -1205.53503417969 -1306.20495605469 -1409.15002441406  -1517.09497070313 -1634.17504882813 -1765.13500976563 -1914.15002441406 -2084.03491210938 -2276.22509765625 -2491.25000000000 -2729.25000000000 -2990.25000000000 -3274.25000000000 -3581.25000000000 -3911.25000000000 -4264.25000000000 -4640.25000000000 -5039.25000000000 -5461.25000000000 -5906.25000000000]')
n1,n2,n3 = (length(lon),length(lat),length(dep))

dimdict = Dict{AbstractString,NCvar}([
    "lon_c" => NCvar("lon_c","degrees east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
    "lat_c" => NCvar("lat_c","degrees north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
    "dep_c" => NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
    "time" => NCvar("time",timeunits,Inf,timesteps,Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)
])

dims = collect(values(dimdict))

for group in groups
    println(group)

    selectfields = filter(x -> isdir(interpdir*group*"/"*x),readdir(interpdir*group))

    flds = Dict{AbstractString,NCvar}()

    for fldname in selectfields
        
        fullpath = interpdir*group*'/'*fldname*'/'
        fnames = fullpath.*filter(x -> occursin(".data",x), readdir(fullpath))
        metafile = fullpath*filter(x -> occursin(".meta",x), readdir(fullpath))[1]
        metafile = parsemeta(metafile)
        diaginfo = readAvailDiagnosticsLog(availdiagsfname,fldname)


        prec = eval(Meta.parse(metafile["dataprec"]))
        units = diaginfo["units"]
        longname = diaginfo["title"]

        mydims = Array{NCvar,1}()
        for i in 1:metafile["nDims"]
            dimidx = getindex.(getfield.(dims,:dims),1) .== metafile["dimList"][i,1]
            mydims = cat(mydims,dims[dimidx],dims=1)
        end

        iosize = Tuple([d.dims[1] for d in mydims])

        mydims = cat(mydims,dimdict["time"],dims=1)

        fielddata = Bindata(fnames,prec,iosize)
        flds[fldname] = NCvar(fldname,units,mydims,fielddata,
            Dict("long_name" => longname),NCDatasets)

    end

    # Set up NCvars
    filename = ncdir*"/"*group*".nc"
    
    # Create the NetCDF file and populate with dimension and field info
    ds,fieldvars,dimlist = createfile(filename,flds,README)

    # Add field and dimension data and close the file
    for k in keys(flds)
        addData(ds[k],flds[k])
    end

    addDimData.(Ref(ds),dims)
    close(ds)

    # Create NetCDF file for each time step
    for i in 1:length(timesteps)
        global timestart,timesteps,fnameprefix,README
        thistime = timestart + Dates.Day(floor(timesteps[i])-1)

        yyyyddd = @sprintf("%i%03i",year(thistime),dayofyear(thistime))

        thisfilename = join((fnameprefix,lowercase(group),yyyyddd),"_")*".nc"

        ds,fieldvars,dimlist = createfile(thisfilename,flds,README)

        # Add field and dimension data and close the file
        for k in keys(flds)
            addData(ds[k],flds[k][i])
        end

        if haskey(dimdict,"time")
            hastimedim = true
            timedim = pop!(dimdict,"time")
        else
            hastimedim = false
        end

        addDimData.(Ref(ds),collect(values(dimdict)))
        if hastimedim
            addDimData(ds,timedim[i])
            dimdict["time"] = timedim
        end

        close(ds)

    end
end
