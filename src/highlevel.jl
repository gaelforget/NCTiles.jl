using NCTiles,NCDatasets,ClimateTools


"""
writeNetCDFtiles(flds::Dict,savenamebase::String,README::Array)

Function to write out tiled NetCDF files. Flds should be a Dict of NCVars, 
savenamebase should be the prefix of the filenames to which the tile 
number and file exension is added, including full path to the save 
location, and README should be an Array of strings containing the
description to write into the files.
"""
function writeNetCDFtiles(flds::Dict,savenamebase::String,README::Array)

    fldnames = collect(keys(flds))
    tilefld = flds[fldnames[findfirst(isa.([flds[f].values for f in fldnames],TileData))]]
    numtiles = tilefld.values.numtiles

    landidx = findfirst(get.([flds[f].atts for f in fldnames],"standard_name","none").=="land_binary_mask")
    if ~isnothing(landidx); land_mask = flds[fldnames[landidx]].values; else land_mask = nothing; end
    if isa(land_mask,TileData); land_mask = land_mask.vals; end
    savenames = savenamebase*".".*lpad.(string.(1:numtiles),4,"0").*".nc"

    datasets = [createfile(savenames[tidx],flds,README, itile = tidx, ntile = length(savenames)) for tidx in 1:length(savenames)]

    ds = [x[1] for x in datasets]
    fldvars = [x[2] for x in datasets]

    for k in keys(flds)
        if isa(flds[k].values,TileData)
            addData(fldvars,flds[k],land_mask = land_mask)
        else
            tmpfldvars = [fv[findfirst(isequal(k),name.(fv))] for fv in fldvars]
            addData.(tmpfldvars,Ref(flds[k]))
        end
    end

    dims = unique(vcat([flds[v].dims for v in keys(flds)]...))
    ims = filter( d -> isa(d,NCvar),dims)

    for dim in dims
        addDimData.(ds,Ref(dim))
    end

    close.(ds);

    return nothing

end
