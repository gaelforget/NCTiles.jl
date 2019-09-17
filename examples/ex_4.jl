using MeshArrays, NCDatasets, NCTiles
grid = GridSpec("LLC90","grids/")

exampledir = joinpath("data","ex4")
datadir = joinpath(exampledir,"model_output")
saveloc = joinpath(exampledir,"tiles")
fnames = joinpath.(Ref(datadir),filter(x->occursin("data",x),readdir(datadir)))

metafile = parsemeta(joinpath(datadir,"ptr_3d_set1.0000000732.meta"))
README = readlines(joinpath(datadir,"README"))

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
    NCvar("tim",timeunits,Inf,time_steps,Dict(("long_name" => "time","standard_name" => "time")),NCDatasets)
]

gridvars = GridLoad(grid)

land = gridvars["hFacC"]
for f in land.fIndex
    for d in 1:size(land,2)
        land[f,d][land[f,d].==0] .= NaN
        land[f,d][land[f,d].>0] .= 1
    end
end

tilarea = TileData(gridvars["RAC"],tilesize,grid)
tilland = TileData(land,tilesize,grid)
thic = gridvars["DRC"]

fldidx = 1:106
for fidx in fldidx
    @time begin
        fldname = metafile["fldList"][fidx]
        println("Processing "*fldname)
        diaginfo = readAvailDiagnosticsLog(joinpath(exampledir,"available_diagnostics.log"),fldname)
        fielddata = BinData(fnames,prec,iosize,fidx)
        tilfld = TileData(fielddata,tilesize,grid)
        tillat = TileData(gridvars["YC"],tilfld.tileinfo,tilfld.tilesize,tilfld.precision,tilfld.numtiles)
        tillon = TileData(gridvars["XC"],tilfld.tileinfo,tilfld.tilesize,tilfld.precision,tilfld.numtiles)
        flds = Dict([fldname => NCvar(fldname,diaginfo["units"],dims,tilfld,Dict(),NCDatasets),
                    "lon" => NCvar("lon","degrees_east",dims[1:2],tillon,Dict("long_name" => "longitude"),NCDatasets),
                    "lat" => NCvar("lat","degrees_north",dims[1:2],tillat,Dict("long_name" => "latitude"),NCDatasets),
                    "area" => NCvar("area","m^2",dims[1:2],tilarea,Dict(["long_name" => "grid cell area", "standard_name" => "cell_area"]),NCDatasets),
                    "land" => NCvar("land","1",dims[1:3],tilland,Dict(["long_name" => "land mask", "standard_name" => "land_binary_mask"]),NCDatasets),
                    "thic" => NCvar("thic","m",dims[3],thic,Dict("standard_name" => "cell_thickness"),NCDatasets)
        ]) 
        savepath = joinpath(saveloc,fldname)
        if !isdir(savepath); mkpath(savepath); end

        numtiles = flds[fldname].values.numtiles
        savenames = joinpath.(Ref(savepath),fldname*".".*lpad.(string.(1:numtiles),4,"0").*".nc")
        rm.(savenames, force=true)

        datasets = [createfile(savenames[tidx],flds,README, itile = tidx, ntile = length(savenames)) for tidx in 1:length(savenames)]

        ds = [x[1] for x in datasets]
        fldvars = [x[2] for x in datasets]
        #dims = [x[3] for x in datasets]

        for k in keys(flds)
            if isa(flds[k].values,TileData)
                addData(fldvars,flds[k])
            else
                tmpfldvars = [fv[findfirst(isequal(k),name.(fv))] for fv in fldvars]
                addData.(tmpfldvars,Ref(flds[k]))
            end
        end

        for dim in dims
            addDimData.(ds,Ref(dim))
        end

        close.(ds)
    end
end