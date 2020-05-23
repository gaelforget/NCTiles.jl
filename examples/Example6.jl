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

# # Example 6

# +
using NCTiles,NCDatasets

# File Paths
inputs = "input/"
NCTiles.get_testcases_if_needed(inputs)

outputs = "output/"
if ~ispath(outputs); mkpath(outputs); end

savedir = joinpath(outputs,"ex6")
if ~ispath(savedir); mkpath(savedir); end

file_in="FeT.0001.nc"
field_name = "FeT"
README = [field_name*" -- Source: Gael Forget; version: alpha."]

#srv="http://engaging-opendap.mit.edu:8080/thredds/dodsC/las/"
#fil="id-2e0ea5ca2c/data_usr_local_tomcat_content_cbiomes_20200206_17_Nutrients_FeT.0001.nc.jnl"
#ncdata=Dataset(srv*fil)

ncdata=Dataset(joinpath(inputs*file_in))

# +
#ncvars,ncdims,fileatts = readncfile(joinpath(inputs*file_in))
#rm(joinpath(savedir,"ex6a.nc"),force=true)
#write(ncvars,joinpath(savedir,"ex6a.nc"),README=README)

# +
lon=ncdata["lon_c"][:]
lat=ncdata["lat_c"][:]
dep_c=ncdata["dep_c"][:]
tim=ncdata["tim"][:]
tim_units=get(ncdata["tim"].attrib,"units","")

dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("dep_c","m",size(dep_c),dep_c,Dict(["long_name" => "depth","positive"=>"down","standard_name"=>"depth"]),NCDatasets),
        NCvar("tim",tim_units,Inf,tim,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]

v=ncdata[field_name]
u=get(v.attrib,"units","")
long_name=get(v.attrib,"long_name","")

vv=v[:,:,:,:]
vv[findall(vv.<-1.0e33)].=NaN
field = NCvar(field_name,u,dims,vv,Dict("long_name" => long_name),NCDatasets)

rm(joinpath(savedir,"ex6a.nc"),force=true)
write(field,joinpath(savedir,"ex6a.nc"),README=README)
# + {}
ncvars,ncdims,fileatts = readncfile(joinpath(inputs*"FeT.0062.nc"))
rm(joinpath(savedir,"ex6b.nc"),force=true)
write(ncvars,joinpath(savedir,"ex6b.nc"),README=README)

climatology_bounds = NCvar(ncvars["climatology_bounds"].name,
                                        ncvars["climatology_bounds"].units,
                                        [ncdims["tcb"],dims[end]],
                                        ncvars["climatology_bounds"].values,
                                        ncvars["climatology_bounds"].atts,
                                        ncvars["climatology_bounds"].backend)
# -


