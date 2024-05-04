
# # Example 1 : Regular Grid
#
# - two-dimensional fields on a regular `longitude,latitude` grid are retrieved from binary files
# - they are then written to a new netcdf file using `NCTiles.write()`
# 
# This example illustrates the use of either `NCDatasets.jl` or `NetCDF.jl` as the backend.

using NCTiles
import NCTiles: NCDatasets, NetCDF, Dates

# ### File Paths
inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
if ~ispath(outputs); mkpath(outputs); end

selectfields = ["ETAN"]
indir = joinpath(inputs,"diags_interp",selectfields[1])
fnames = joinpath.(Ref(indir),filter(x -> occursin(".data",x), readdir(indir)))

savedir = joinpath(outputs,"ex1")
if ~ispath(savedir); mkpath(savedir); end

# ### Grid specifications & Metadata

prec = Float32
README = ["File created by","example 1 of NCTiles.jl","on "*string(Dates.now())]
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
n1,n2 = (length(lon),length(lat))
tim_units = "days since 1992-01-01"
tim = vec([14.0 45.0 74.0])
longname = "Surface Height Anomaly"
units = "m"

# ### 1. Using `NCDatasets.jl` As The Backend

# Define Dimensions as `NCvar`s
dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("tim",tim_units,Inf,tim,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]

# Define Variable Array as `BinData`
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NCDatasets)

# Create the NetCDF file and populate with dimensions and variable
write(field,joinpath(savedir,"ex1_NCDatasets.nc"),README=README)

# ### 2. Using `NetCDF.jl` As The Backend

# Define Dimensions as `NCvar`s
dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NetCDF),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NetCDF),
        NCvar("tim",tim_units,Inf,tim,Dict(("long_name" => "time","standard_name" => "time")),NetCDF)
        ]

# Define Variable Array as `BinData`
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NetCDF)

write(field,joinpath(savedir,"ex1_NetCDF.nc"),README=README)

# _Note: the `write` function is a shorthand for_
#
# ```
# # Create the NetCDF file and populate with dimension and field info, as well as dimension data
# ncfile,fldvar,dimlist = createfile(joinpath(savedir,"ex1_NetCDF.nc"),field,README)
#
# # Add field data
# addData(fldvar,field)
#
# # Close the file-  only needed for NetCDF v < 0.10
# finalize(ncfile)
# ```



