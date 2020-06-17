using NCTiles,NCDatasets,NetCDF,Dates,Printf,MITgcmTools

examplesdir = joinpath("data","ex3")
indir = joinpath(examplesdir,"data","cs510","diags_interp")
availdiagsfname = joinpath(examplesdir,"data","available_diagnostics.log")
savedir = joinpath(examplesdir,"interp_nctiles")

if ~ispath(savedir); mkpath(savedir); end
README = readlines(joinpath(examplesdir,"README"))

exdir = "../../../run/run_cs510/"
fnameprefix = "darwin_v0.2_cs510"
groups = ["Nutrients","Phytoplankton_Functional_Types","Ocean_Color"]

rename = Dict(["TRAC05" => "PO4",
"TRAC06" => "SiO2",
"TRAC07" => "FeT", 
"TRAC19" => "O2"])

unitsrename = Dict(["mmol C/" => "mmol N/m^3",
                    "mmol O/" => "mmol O/m^3",
                    "mmol P/" => "mmol P/m^3",
                    "mmol Si" => "mmol Si/m^3",
                    "mmol Fe" => "mmol Fe/m^3"])


# Hard coding dims for now
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;

timeinterval = 3
tstart = 732; tend = 8766;
#tstart = 6; tend = 6;
timesteps = ((tstart - timeinterval/2):timeinterval:(tend - timeinterval/2))[1:4]
nsteps = length(timesteps)
timeinfname = lpad.(string.(Integer.(collect((timesteps.+timeinterval/2).*24))),10,'0')
timeunits = "days since 1992-1-1 0:0:0"
timestart = DateTime(1992,01,01)
dep = -1 .* collect([-5 -15 -25 -35 -45 -55 -65 -75.0049972534180 -85.0250015258789 -95.0950012207031 -105.309997558594 -115.870002746582 -127.150001525879 -139.740005493164 -154.470001220703 -172.399993896484 -194.735000610352 -222.710006713867 -257.470001220703 -299.929992675781 -350.679992675781 -409.929992675781 -477.470001220703 -552.710021972656 -634.734985351563 -722.400024414063 -814.469970703125 -909.739990234375 -1007.15502929688 -1105.90502929688 -1205.53503417969 -1306.20495605469 -1409.15002441406  -1517.09497070313 -1634.17504882813 -1765.13500976563 -1914.15002441406 -2084.03491210938 -2276.22509765625 -2491.25000000000 -2729.25000000000 -2990.25000000000 -3274.25000000000 -3581.25000000000 -3911.25000000000 -4264.25000000000 -4640.25000000000 -5039.25000000000 -5461.25000000000 -5906.25000000000]')
n1,n2,n3 = (length(lon),length(lat),length(dep))

dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
        NCvar("time",timeunits,Inf,timesteps,Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)]

dims_NetCDF = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NetCDF),
                NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NetCDF),
                NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NetCDF),
                NCvar("time",timeunits,Inf,timesteps,Dict(("long_name" => "time","standard_name" => "time")),NetCDF)]

for group in groups
    println(group)
    groupindir = joinpath(indir,group)
    groupsavedir = joinpath(savedir,group)
    if !isdir(groupsavedir); mkpath(groupsavedir); end

    selectfields = filter(x -> isdir(joinpath(groupindir,x)),readdir(groupindir))

    flds = Dict{AbstractString,NCvar}()
    flds_NetCDF = Dict{AbstractString,NCvar}()

    timestepidx = []

    for fldname in selectfields
        fullpath = joinpath(groupindir,fldname)
        fnames_timerange = fldname*'.'.*timeinfname.*".data"
        fnames_avail = filter(x -> occursin(".data",x), readdir(fullpath))
        fnames = joinpath.(Ref(fullpath),intersect(fnames_timerange,fnames_avail))
        timestepidx = [findfirst(occursin.(f,fnames_timerange)) for f in fnames_avail]
        timestepidx = timestepidx[.!isa.(timestepidx,Nothing)]
        metafile = joinpath.(Ref(fullpath),fldname*'.'.*timeinfname.*".meta")
        #filter(x -> occursin(".meta",x), readdir(fullpath))[1]
        metafile = parsemeta(metafile[1])
        diaginfo = readAvailDiagnosticsLog(availdiagsfname,fldname)

        #timesteps = parse.(Int,getindex.(split.(fnames,'.'),2))./24 .- timeinterval/2

        # Rename if needed
        if any(occursin.(fldname,collect(keys(rename))))
            fldname = rename[fldname]
        end

        prec = eval(Meta.parse(metafile["dataprec"]))
        units = get(unitsrename,diaginfo["units"],diaginfo["units"])
        longname = diaginfo["title"]

        if metafile["nDims"] == 3
            mydims = dims
            mydims_NetCDF = dims_NetCDF
        else
            mydims = dims[[1; 2; 4]]
            mydims_NetCDF = dims_NetCDF[[1; 2; 4]]
        end

        iosize = Tuple([d.dims[1] for d in mydims[1:end-1]])

        fielddata = BinData(fnames,prec,iosize)
        flds[fldname] = NCvar(fldname,units,mydims,fielddata,
            Dict("long_name" => longname),NCDatasets)
        flds_NetCDF[fldname] = NCvar(fldname,units,mydims_NetCDF,fielddata,
            Dict("long_name" => longname),NetCDF)

    end


    # Set up NCvars
    filename_NCDatasets = joinpath(savedir,join([fnameprefix,lowercase(group)],"_")*"_NCDatasets.nc")
    filename_NetCDF = joinpath(savedir,join([fnameprefix,lowercase(group)],"_")*"_NetCDF.nc")
    
    rm(filename_NCDatasets,force=true)
    rm(filename_NetCDF,force=true)
    println(filename_NCDatasets)
    write(flds,filename_NCDatasets,README=README)
    write(flds_NetCDF,filename_NetCDF,README=README)
    
#= Line above is shorthand for:
    # Create the NetCDF file and populate with dimension and field info
    ds,fieldvars,dimlist = createfile(filename,flds,README)

    # Add field and dimension data and close the file
    for k in keys(flds)
        addData(ds[k],flds[k])
    end

    addDimData.(Ref(ds),dims)
    close(ds)
=#
end







