using NCTiles,NCDatasets

# Point to NetCDF file
examplesdir = joinpath("data","ex2")
indir = joinpath(examplesdir,"nctiles_in")
savedir = joinpath(examplesdir,"nctiles_out")
selectfields = ["THETA"]
fname = joinpath(indir,"THETA.0009.nc")
if ~ispath(savedir); mkpath(savedir); end
README = readlines(joinpath(examplesdir,"README"))

# Get all the metadata from the file and set up NCvars
ncvars,ncdims,fileatts = readncfile(fname)

# Rewrite to a file
write(ncvars,joinpath(savedir,"ex2.nc"),README=README,globalattribs=fileatts)

#= Line above is shorthand for:
# Create the NetCDF file and populate with dimension and field info
ds,fieldvars,dimlist = createfile(filename,ncvars,README,fillval = f_fillval,missval=f_missval,itile=f_ff,ntile=f_ntile)

# Add field and dimension data and close the file
addData(ds["THETA"],ncvars["THETA"])
addDimData.(Ref(ds),collect(values(ncdims)))
close(ds)
=#