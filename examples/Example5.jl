# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
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
using ClimateTools,NCDatasets,NCTiles

# File Paths
inputs = "input/"
NCTiles.get_testcases_if_needed(inputs)

outputs = "output/"
if ~ispath(outputs); mkpath(outputs); end

savedir = joinpath(outputs,"ex5")
if ~ispath(savedir); mkpath(savedir); end

#original climate model output file 
file_in="tas_day_MIROC5_piControl_r1i1p1_20000101-20091231.nc"
field_name = "tas"

if ~isfile(joinpath(inputs,file_in))
    run(`wget http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC5/piControl/day/atmos/day/r1i1p1/v20161012/tas/tas_day_MIROC5_piControl_r1i1p1_20000101-20091231.nc`)
    run(`mv tas_day_MIROC5_piControl_r1i1p1_20000101-20091231.nc ../inputs/`)
end
# -


"""
        climgridtoncvar(C::ClimGrid,N::String)

Creates an NCvar struct from a ClimGrid object. NCvar struct can then be written
to a NetCDF file using write().

Ex: C = load(fname,fldname)
writefld = climgridtoncvar(C,fldname)
write(writefld,"myfile.nc")
"""
function climgridtoncvar(C::ClimGrid,N::String)
        x, y, timevec = ClimateTools.getdims(C) # may need to check number of dims first
        timevec = NCDatasets.timeencode(timevec, C.timeattrib["units"], get(C.timeattrib,"calendar","standard"))
        
        
        dims = [NCvar(C.dimension_dict["lon"],C.lonunits,size(C.data)[1],x,Dict("long_name" => "longitude"),NCDatasets),
                NCvar(C.dimension_dict["lat"],C.latunits,size(C.data)[2],y,Dict("long_name" => "latitude"),NCDatasets),
                NCvar("time",C.timeattrib["units"],Inf,timevec,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
                ]
        
        return NCvar(N,C.dataunits,dims,C.data.data,C.varattribs,NCDatasets)
end   

# Whole Data Set
C = load(joinpath(inputs,file_in),field_name)
writefld = climgridtoncvar(C,field_name)
write(writefld,joinpath(savedir,"ex5_whole.nc"))

# Extracted Sub-Set
poly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]
E = load(joinpath(inputs,file_in),field_name, poly=poly_reg)
writefld = climgridtoncvar(E,field_name)
write(writefld,joinpath(savedir,"ex5_extract.nc"),
    globalattribs=E.globalattribs)


