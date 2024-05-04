# # Example 4 : Group Variables
#
# Two variables (three-dimensional ones) are read from the netcdf files generated in `Example3`, combined into a single data structure, and then re-written together into a new netcdf file.

using NCTiles
import NCTiles: NCDatasets, Dates

inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()
outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
nt="0003"

file_in1=outputs*"ex3/THETA/THETA.$nt.nc"
file_in2=outputs*"ex3/SALT/SALT.$nt.nc"
~isfile(file_in1) ? error("Running Example3 first is needed to run Example4") : nothing

file_out=outputs*"ex4/TS.$nt.nc"
if ~ispath(outputs*"ex4/"); mkpath(outputs*"ex4/"); end

README = ["File created by","example 4 of NCTiles.jl","on "*string(Dates.now())]

# Get the first 3D variable
ncvars,ncdims,fileatts = readncfile(file_in1)

# Add a second 3D variable
tmp,_,_ = readncfile(file_in2)
T=ncvars["THETA"]
S=tmp["SALT"]
ncvars["SALT"]=NCvar(S.name,S.units,T.dims,S.values,S.atts,T.backend)

# Rewrite to a new file
write(ncvars,file_out,README=README,globalattribs=fileatts)

