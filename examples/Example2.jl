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

# # Example 2
#
# Two-dimensional fields are read from the netcdf file generated in `Example1`, and then re-written to a new netcdf file.

# +
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
file_in=outputs*"ex1/ex1_NetCDF.nc"
~isfile(file_in) ? error("Running Example1 first is needed to run Example2") : nothing

file_out=outputs*"ex2/ex2_NetCDF.nc"
if ~ispath(outputs*"ex2/"); mkpath(outputs*"ex2/"); end

README = readlines(joinpath(inputs,"README"))

# Get all the metadata from the file and set up NCvars
ncvars,ncdims,fileatts = readncfile(file_in,NetCDF)

# Rewrite to a file
write(ncvars,file_out,README=README,globalattribs=fileatts)

# -


