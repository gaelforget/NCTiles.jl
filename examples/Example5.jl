# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# # Example 5
#
# A `ClimGrid` struct from the `ClimateTools.jl` package, including metadata read from file, is written to a netcdf file via `NCTiles.jl`. 
#
# _Note: `ClimateTools.jl` requires that requires specific dimension names as follows_
#
# - latitude: lat, latitude, rlat, y, yc
# - longitude: lon, longitude, rlon, x, xc
# - time: time

# +
using ClimateBase, NCDatasets, NCTiles

# File Paths
inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
if ~ispath(outputs); mkpath(outputs); end

savedir = joinpath(outputs,"ex5")
if ~ispath(savedir); mkpath(savedir); end

"""
        climgridtoncvar(C::ClimGrid,N::String)

Creates an NCvar struct from a ClimArray object. NCvar struct can then be written
to a NetCDF file using write().

Ex: C = ClimateBase.ncread(fil, "ETAN")
writefld = ClimArray_to_NCvar(C,"ETAN")
write(writefld,"myfile.nc")
"""
function ClimArray_to_NCvar(C::ClimArray,N::String)
        x = C.dims[1][:]
        y = C.dims[2][:]
        timevec = DateTime.(C.dims[3][:])
        timeunit = "days since 1992-01-01"
        timevec = NCDatasets.timeencode(timevec, timeunit)
        
        dims = [NCvar("lon","degrees_east",size(C.data)[1],x,Dict("long_name" => "longitude"),NCDatasets),
                NCvar("lat","degrees_north",size(C.data)[2],y,Dict("long_name" => "latitude"),NCDatasets),
                NCvar("time",timeunit,Inf,timevec,Dict(("long_name" => "Ti","standard_name" => "time")),NCDatasets)
                ]
        
        return NCvar(N,C.attrib["units"],dims,C.data,C.attrib,NCDatasets)
end   

## Read via NCTiles.jl and write via ClimateBase.jl

fil=joinpath(outputs,"ex1/ex1_NetCDF.nc")
ncvars,ncdims,fileatts = readncfile(fil)

#filout=joinpath(savedir,"ex5_NCTiles.nc")
#write(ncvars["ETAN"],filout,globalattribs=fileatts)

lons=ncdims["lon_c"].values[:];
lats=ncdims["lat_c"].values[:];
t=ncdims["tim"].values[:];
dimensions = (Lon(lons), Lat(lats), ClimateBase.Ti(t))

data=ncvars["ETAN"].values[:];
units=ncvars["ETAN"].units;
A = ClimArray(data, dimensions; name = "ETAN", attrib = Dict("units" => units))

filout=joinpath(savedir,"ex5_ClimateBase.nc")
ncwrite(filout, A)

## Read via ClimateBase.jl and write via NCTiles.jl

ETAN = ClimateBase.ncread(fil, "ETAN")
ETAN = ClimArray_to_NCvar(ETAN,"ETAN")
filout=joinpath(savedir,"ex5_ClimArray_to_NCvar.nc")
write(ETAN,filout)
