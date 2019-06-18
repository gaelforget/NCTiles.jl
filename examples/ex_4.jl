using MeshArrays, NCDatasets, NCTiles
GCMGridSpec("LLC90","../../../grids/")

datadir = "../../../DarwinModelOutputSamples/sample1/output/"
saveloc = "../../../DarwinModelOutputSamples/sample1/tiles/"
fnames = datadir.*filter(x->occursin("data",x),readdir(datadir))

metaname = "../../../DarwinModelOutputSamples/sample1/output/ptr_3d_set1.0000000732.meta"
metafile = parsemeta(metaname)
README = readlines("README")

timeunits = "days since 1992-1-1 0:0:0"
nsteps = 12
timeinterval = 30
time_steps = timeinterval/2:timeinterval:timeinterval*nsteps
dep = -1 .* collect([-5 -15 -25 -35 -45 -55 -65 -75.0049972534180 -85.0250015258789 -95.0950012207031 -105.309997558594 -115.870002746582 -127.150001525879 -139.740005493164 -154.470001220703 -172.399993896484 -194.735000610352 -222.710006713867 -257.470001220703 -299.929992675781 -350.679992675781 -409.929992675781 -477.470001220703 -552.710021972656 -634.734985351563 -722.400024414063 -814.469970703125 -909.739990234375 -1007.15502929688 -1105.90502929688 -1205.53503417969 -1306.20495605469 -1409.15002441406  -1517.09497070313 -1634.17504882813 -1765.13500976563 -1914.15002441406 -2084.03491210938 -2276.22509765625 -2491.25000000000 -2729.25000000000 -2990.25000000000 -3274.25000000000 -3581.25000000000 -3911.25000000000 -4264.25000000000 -4640.25000000000 -5039.25000000000 -5461.25000000000 -5906.25000000000]')
prec = Float32
iosize = (90,1170,50)
tilesize = (90,90)
dims = [
    NCvar("i_c","1",tilesize[1],1:tilesize[1],Dict("long_name" => "Cartesian coordinate 1"),NCDatasets),
    NCvar("j_c","1",tilesize[2],1:tilesize[2],Dict("long_name" => "Cartesian coordinate 2"),NCDatasets),
    NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
    NCvar("time",timeunits,Inf,time_steps,Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)
]

fldidx = 1:106
for fidx in fldidx
    @time begin
        global dims
        fldname = metafile["fldList"][fidx]
        println("Processing "*fldname)
        diaginfo = readAvailDiagnosticsLog("available_diagnostics.log",fldname)
        fielddata = BinData(fnames,prec,iosize,fidx)
        tilfld = TileData(fielddata,tilesize)
        tillat = TileData(tilfld.tileinfo["YC"],tilfld.tileinfo,tilfld.tilesize,tilfld.precision,tilfld.numtiles)
        tillon = TileData(tilfld.tileinfo["XC"],tilfld.tileinfo,tilfld.tilesize,tilfld.precision,tilfld.numtiles)
        flds = Dict([fldname => NCvar(fldname,diaginfo["units"],dims,tilfld,Dict(),NCDatasets),
                    "lon_c" => NCvar("lon_c","degrees east",dims[1:2],tillon,Dict("long_name" => "longitude"),NCDatasets),
                    "lat_c" => NCvar("lat_c","degrees north",dims[1:2],tillat,Dict("long_name" => "latitude"),NCDatasets)
        ]) 
        #myfld = NCvar(fldname,diaginfo["units"],dims,TileData(flds[1],tilesize),Dict(),NCDatasets)
        savepath = saveloc*fldname
        if !isdir(savepath)
            mkpath(savepath)
        end

        numtiles = flds[fldname].values.numtiles
        savenames = savepath*"/"*fldname*".".*lpad.(string.(1:numtiles),ndigits(numtiles),"0").*".nc"

        #savenames = savepath*"/ex2_Tiles_".*string.(1:flds[fldname].values.numtiles).*".nc"

        datasets = createfile.(savenames,Ref(flds),Ref(README))

        ds = [x[1] for x in datasets]
        fldvars = [x[2] for x in datasets]
        #dims = [x[3] for x in datasets]

        addData(fldvars,flds[fldname])
        addData(fldvars,flds["lat_c"])
        addData(fldvars,flds["lon_c"])

        for dim in dims
            addDimData.(ds,Ref(dim))
        end

        close.(ds)
    end
end