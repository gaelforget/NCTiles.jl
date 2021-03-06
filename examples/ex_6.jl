using NCTiles,NCDatasets

# Climatology example: This example has two parts. First, we take a NetCDF file with climatology_bounds
# and rewrite it as is (6a). Then we are write a new NetCDF file with updated dimensions (6b).
# This example is based on the processing that had to be done to put this data into
# LAS/OpenDap. In order to display LAS prefers the data have dimensions lon,lat(,dep),time and that each
# of these be one dimensional (linear), whereas the data is originally formatted to have dimensions
# i,j(,k),t with 2D latitude and longitude added as variables (curvelinear). Additionally, the attribute
# "positive" needed to be added to the depth dimension. This is climatology data and includes
# the climatology_bounds.

# Point to NetCDF file
examplesdir = joinpath("data","ex6")
indir = joinpath(examplesdir,"input")
savedir = joinpath(examplesdir,"output")
if ~ispath(savedir); mkpath(savedir); end
fldname = "FeT"
fname = joinpath(indir,fldname*".0001.nc")
README = [fldname*" -- Source: Gael Forget; version: alpha."]

# Read in current NetCDF File
ncvars,ncdims,fileatts = readncfile(fname)

# Ex 6a: Rewrite as-is
rm(joinpath(savedir,"ex6a.nc"),force=true)
write(ncvars,joinpath(savedir,"ex6a.nc"),README=README)

# Ex 6b: Rewrite NetCDF file with updated dimensions

# Reworking lat,lon,dep Dimensions
dims = []
dimdict = Dict(["i_" => "lon", "j_" => "lat", "k_" => "dep"])
for d in ncvars[fldname].dims
    global dims
    if haskey(dimdict,d.name[1:end-1])
        newdim = dimdict[d.name[1:end-1]]
        newdimname = join([newdim, d.name[end]],"_")
        attribs = ncvars[newdim].atts
        if newdim == "dep" # Add "positive" => "down" to depth
            attribs["positive"] = "down"
            attribs["standard_name"] = "depth"
        end
        # Grab one-dimensional values for lat/lon
        vals = newdim == "lat" ? ncvars[newdim].values[1,:] : ncvars[newdim].values[:,1]

        # Add updated dimension to dims array
        dims = vcat(dims,NCvar(newdimname,
                                ncvars[newdim].units,
                                ncdims[d.name].dims,
                                vals,
                                attribs, 
                                NCDatasets))
        pop!(ncvars,newdim) #remove extra dimension variable
    end
end

# Reworking Time Dimension
dims = vcat(dims,NCvar("tim",
                        ncvars["tim"].units,
                        ncdims["t"].dims,
                        Array{DateTime}(ncvars["tim"].values[:]),
                        ncvars["tim"].atts,
                        NCDatasets))
pop!(ncvars,"tim")

# Replacing dimensions in the NCvars
ncvars[fldname] = NCvar(ncvars[fldname].name,
                                ncvars[fldname].units,
                                dims,
                                ncvars[fldname].values,
                                ncvars[fldname].atts,
                                ncvars[fldname].backend)
ncvars["climatology_bounds"] = NCvar(ncvars["climatology_bounds"].name,
                                        ncvars["climatology_bounds"].units,
                                        [ncdims["tcb"],dims[end]],
                                        ncvars["climatology_bounds"].values,
                                        ncvars["climatology_bounds"].atts,
                                        ncvars["climatology_bounds"].backend)

    
# Write result to a new NetCDF file
rm(joinpath(savedir,"ex6b.nc"),force=true)
write(ncvars,joinpath(savedir,"ex6b.nc"),README=README)