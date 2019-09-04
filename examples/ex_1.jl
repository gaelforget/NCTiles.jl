using NCTiles,NCDatasets,NetCDF

# Point to interpolated 2D Data
examplesdir = joinpath("data","ex1")
selectfields = ["Chl050"]
indir = joinpath(examplesdir,"diags_interp",selectfields[1])
savedir = joinpath(examplesdir,"interp_ncfiles")
if ~ispath(savedir); mkpath(savedir); end
fnames = joinpath.(Ref(indir),filter(x -> occursin(".data",x), readdir(indir)))

# Hard coding these for now
prec = Float32
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
n1,n2 = (length(lon),length(lat))
nsteps = 2
timeinterval = 3
time_steps = timeinterval/2:timeinterval:timeinterval*nsteps
longname = "Average chlorophyll concentration (top 50m)"
units = "mg Chl"

README = readlines(joinpath(examplesdir,"README"))

# Using NCDatasets

# Set up data for writing to NCfiles
# Define dimensions
dims = [NCvar("lon_c","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat_c","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("tim","days since 1992-01-01",Inf,time_steps,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]

# Define field- BinData contains the filenames where the data sits so it's only loaded when needed
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NCDatasets)

# Create the NetCDF file and populate with dimension and field info
ds,fldvar,dimlist = createfile(joinpath(savedir,"ex1_NCDatasets.nc"),field,README)

# Add field and dimension data
addData(fldvar,field)
addDimData.(Ref(ds),field.dims)

# Close the file
close(ds)

# Using NetCDF

# Set up data for writing to NCfiles
# Define dimensions
dims = [NCvar("lon_c","degrees east",size(lon),lon,Dict("long_name" => "longitude"),NetCDF),
        NCvar("lat_c","degrees north",size(lat),lat,Dict("long_name" => "latitude"),NetCDF),
        NCvar("tim","days from 1992-01-01",Inf,collect(time_steps),Dict(("long_name" => "time","standard_name" => "time")),NetCDF)
        ]

# Define field- BinData contains the filenames where the data sits so it's only loaded when needed
fielddata = BinData(fnames,prec,(n1,n2))
field = NCvar(selectfields[1],units,dims,fielddata,Dict("long_name" => longname),NetCDF)

using NetCDF

# Create the NetCDF file and populate with dimension and field info, as well as dimension data
ncfile,fldvar,dimlist = createfile(joinpath(savedir,"ex1_NetCDF.nc"),field,README)

# Add field data
addData(fldvar,field)

# Close the file
NetCDF.close(ncfile)
