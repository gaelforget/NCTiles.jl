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

# # Example 6
#
# An example of one tile with a climatology time axis from the global ocean domain available in https://github.com/gaelforget/nctiles-testcases
#

# +
using NCTiles,NCDatasets

# File Paths
inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
if ~ispath(outputs); mkpath(outputs); end

savedir = joinpath(outputs,"ex6")
if ~ispath(savedir); mkpath(savedir); end

field_name = "FeT"
README = [field_name*" -- Source: Gael Forget; version: alpha."];
# -

# ### One tile example

ncvars,ncdims,fileatts = readncfile(joinpath(inputs,"diags_nctiles/FeT.0062.nc"))
rm(joinpath(savedir,"ex6a.nc"),force=true)
write(ncvars,joinpath(savedir,"ex6a.nc"),README=README)


# ### Define Dimensions
#
# For reference:
#
#
# ```
# - i_c: NCvar("i_c", "1", 30, 1.0:30.0, Dict("units" => "1","long_name" => "Cartesian coordinate 1"), NCDatasets) 
# - j_c: NCvar("j_c", "1", 30, 1.0:30.0, Dict("units" => "1","long_name" => "Cartesian coordinate 2"), NCDatasets)
# - k_c: NCvar("k_c", "1", 50, 1.0:50.0, Dict("units" => "1","long_name" => "Cartesian coordinate 3"), NCDatasets)
# - t:   NCvar("t", "1", 12, 1.0:12.0, Dict("units" => "1","long_name" => "Time coordinate"), NCDatasets)
# - tcb:  NCvar("tcb", "", 2, Any[], Dict{Any,Any}(), NCDatasets)
# ```
#
#
# Note : tcb is unitles, only has a dimension, no values

FeT_dims = [ncdims["i_c"],
                ncdims["j_c"],
                ncdims["k_c"],
                ncdims["t"]]
clim_dims = [ncdims["tcb"],
                ncdims["t"]]

# ### Define Variables
#
# Note : lat, lon, dep, tim, thic, area, land defined same as non-climatology example

# +
FeT = NCvar("FeT", # name
                "mmol Fe", # units
                FeT_dims, # dimensions
                NCData("input/diags_nctiles/FeT.0062.nc", "FeT", NCDatasets, Float32), # values- to be read from file
                Dict("coordinates" => "lon lat dep tim","long_name" => "FeT concentration"), # attributes
                NCDatasets) # backend

climatology_bounds = NCvar("climatology_bounds", # name
                                "days since 1992-1-1 0:0:0", # units
                                clim_dims, # dimensions
                                NCData("input/diags_nctiles/FeT.0062.nc", "climatology_bounds", NCDatasets, Float32), # values- to be read from file
                                Dict("long_name" => "climatology_bounds"), # attributes
                                ncvars["climatology_bounds"].backend) # backend
# -

# Write to file

writevars = Dict(["FeT" => FeT,
                "lon" => ncvars["lon"],
                "lat" => ncvars["lat"],
                "dep" => ncvars["dep"],
                "tim" => ncvars["tim"],
                "thic" => ncvars["thic"],
                "area" => ncvars["area"],
                "land" => ncvars["land"],
                "climatology_bounds" => climatology_bounds])
rm(joinpath(savedir,"ex6b.nc"),force=true)
write(writevars,joinpath(savedir,"ex6b.nc"),README=README)


