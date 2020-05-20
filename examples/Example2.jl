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

# # Example 2
#
# Two-dimensional fields are read from the netcdf file generated in `Example1`, and then re-written to a new netcdf file.

# +
using NCTiles,NCDatasets

inputs = "input/"
outputs = "output/"

file_in=outputs*"ex1/ex1_NCDatasets.nc"
~isfile(file_in) ? error("Running Example1 first is needed to run Example2") : nothing

file_out=outputs*"ex2/ex2_NCDatasets.nc"
if ~ispath(outputs*"ex2/"); mkpath(outputs*"ex2/"); end

README = readlines(joinpath(inputs,"README"))

# Get all the metadata from the file and set up NCvars
ncvars,ncdims,fileatts = readncfile(file_in)

# Rewrite to a file
write(ncvars,file_out,README=README,globalattribs=fileatts)
# -


