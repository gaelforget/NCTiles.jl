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

# # Example 4
#
# Two three-dimensional variables are read from the netcdf files generated in `Example3`, combined into a single data structure, and then re-written together into a new netcdf file.

# +
using NCTiles,NCDatasets

inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()
outputs = "output/"
nt="0003"

file_in1=outputs*"ex3/THETA/THETA.$nt.nc"
file_in2=outputs*"ex3/SALT/SALT.$nt.nc"
~isfile(file_in1) ? error("Running Example3 first is needed to run Example4") : nothing

file_out=outputs*"ex4/TS.$nt.nc"
if ~ispath(outputs*"ex4/"); mkpath(outputs*"ex4/"); end

README = readlines(joinpath(inputs,"README"))

# Get the first 3D variable
ncvars,ncdims,fileatts = readncfile(file_in1)

# Add a second 3D variable
tmp,_,_ = readncfile(file_in2)
T=ncvars["THETA"]
S=tmp["SALT"]
ncvars["SALT"]=NCvar(S.name,S.units,T.dims,S.values,S.atts,T.backend)

# Rewrite to a new file
write(ncvars,file_out,README=README,globalattribs=fileatts)
# -

