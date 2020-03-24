using ClimateTools,NCDatasets,NCTiles
# ClimateTools requires that lat,lon,time requires specific field names
# latitude: lat, latitude, rlat, y, yc
# longitude: lon, longitude, rlon, x, xc
# time: time
examplesdir = joinpath("data","ex5")
fldname = "Chl050"
indir = joinpath(examplesdir,"infiles")
savedir = joinpath(examplesdir,"outfiles")
if ~ispath(savedir); mkpath(savedir); end

"""
        climgridtoncvar(C::ClimGrid)

Creates an NCvar struct from a ClimGrid object. NCvar struct can then be written
to a NetCDF file using write().

Ex: C = load(fname,fldname)
writefld = climgridtoncvar(C)
write(writefld,"myfile.nc")
"""
function climgridtoncvar(C::ClimGrid)
        x, y, timevec = ClimateTools.getdims(C) # may need to check number of dims first
        timevec = NCDatasets.timeencode(timevec, C.timeattrib["units"], get(C.timeattrib,"calendar","standard"))
        
        
        dims = [NCvar(C.dimension_dict["lon"],C.lonunits,size(C.data)[1],x,Dict("long_name" => "longitude"),NCDatasets),
                NCvar(C.dimension_dict["lat"],C.latunits,size(C.data)[2],y,Dict("long_name" => "latitude"),NCDatasets),
                NCvar("time",C.timeattrib["units"],Inf,timevec,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
                ]
        
        return NCvar(fldname,C.dataunits,dims,C.data.data,C.varattribs,NCDatasets)
end   

# First test one of our files
fname = joinpath(indir,"Chl050.nc")
README = readlines(joinpath(examplesdir,"README"))
savename = joinpath(savedir,"ex5.nc") 

C = load(fname,fldname)
writefld = climgridtoncvar(C)
write(writefld,savename)

# Files from ClimateTools Example
gcmfiles =["tasmax_day_MIROC5_historical_r1i1p1_19800101-19891231.nc",
"tasmax_day_MIROC5_historical_r1i1p1_19900101-19991231.nc",
"tasmax_day_MIROC5_historical_r1i1p1_20000101-20091231.nc"]
fldname = "tasmax"
savename = joinpath(savedir,gcmfiles[1])
C = load(joinpath(indir,gcmfiles[1]),fldname)
writefld = climgridtoncvar(C)
write(writefld,savename)

# Following ClimateTools Example
# Extraction
poly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]
model = load(joinpath.(Ref(indir),gcmfiles), "tasmax", poly=poly_reg)
writefld = climgridtoncvar(model)
write(writefld,joinpath(savedir,"ex5_extraction.nc"),globalattribs=model.globalattribs)
