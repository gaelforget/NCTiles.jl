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

# # Example 1
#
# Two-dimensional field on a regular ("lat-lon") grid gets written to netcdf file using either `NCDatasets.jl` or `NetCDF.jl` as the backend.

using NCTiles,NCDatasets,NetCDF

# ### File Paths & Names

examplesdir = joinpath("data","ex1")
selectfields = ["Chl050"]
indir = joinpath(examplesdir,"diags_interp",selectfields[1])
savedir = joinpath(examplesdir,"interp_ncfiles")
if ~ispath(savedir); mkpath(savedir); end
fnames = joinpath.(Ref(indir),filter(x -> occursin(".data",x), readdir(indir)))

# ### Grid specifications & Metadata
#

prec = Float32
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
n1,n2 = (length(lon),length(lat))
nsteps = 2
timeinterval = 3
time_steps = timeinterval/2:timeinterval:timeinterval*nsteps
longname = "Average chlorophyll concentration (top 50m)"
units = "mg Chl"
README = readlines(joinpath(examplesdir,"README"))

# ### 1. Using `NCDatasets.jl` For The I/O Backend
#
# The write function used below is a shorthand for:
#
# ```
# # Create the NetCDF file and populate with dimension and field info
# ds,fldvar,dimlist = createfile(joinpath(savedir,"ex1_NCDatasets.nc"),field,README)
#
# # Add field and dimension data
# addData(fldvar,field)
# addDimData.(Ref(ds),field.dims)
#
# # Close file
# close(ds)
# ```

# +
# Define Dimensions as `NCvar`s
dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("tim","days since 1992-01-01",Inf,time_steps,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]

# Define Variable Array as `BinData`
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NCDatasets)
# -

# Create the NetCDF file and populate with dimensions and variable
write(field,joinpath(savedir,"ex1_NCDatasets.nc"),README=README)

# ### 2. Using `NetCDF.jl` For The I/O Backend

# +
# Define Dimensions as `NCvar`s
dims = [NCvar("lon_c","degrees east",size(lon),lon,Dict("long_name" => "longitude"),NetCDF),
        NCvar("lat_c","degrees north",size(lat),lat,Dict("long_name" => "latitude"),NetCDF),
        NCvar("tim","days from 1992-01-01",Inf,collect(time_steps),Dict(("long_name" => "time","standard_name" => "time")),NetCDF)
        ]

# Define Variable Array as `BinData`
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NetCDF)

# Create the NetCDF file and populate with dimension and field info, as well as dimension data
ncfile,fldvar,dimlist = createfile(joinpath(savedir,"ex1_NetCDF.nc"),field,README)

# Add field data
addData(fldvar,field)

# Close the file
#NetCDF.close(ncfile)


# -


