# # Example 2 : Lazily Read
#
# Two-dimensional fields are read from the netcdf file generated in `Example1`, and then re-written to a new netcdf file.

using NCTiles,NetCDF

inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
if ~ispath(outputs); mkpath(outputs); end

# Using NCDatasets backend
file_in=outputs*"ex1/ex1_NCDatasets.nc"
~isfile(file_in) ? error("Running Example1 first is needed to run Example2") : nothing

file_out=outputs*"ex2/ex2_NCDatasets.nc"
if ~ispath(outputs*"ex2/"); mkpath(outputs*"ex2/"); end

README = readlines(joinpath(inputs,"README"))

# Get all the metadata from the file and set up NCvars
ncvars,ncdims,fileatts = readncfile(file_in)

# Rewrite to a file
write(ncvars,file_out,README=README,globalattribs=fileatts)

# Using NetCDF backend
file_out=outputs*"ex2/ex2_NetCDF.nc"
ncvars,ncdims,fileatts = readncfile(file_in,NetCDF)
write(ncvars,file_out,README=README,globalattribs=fileatts)



