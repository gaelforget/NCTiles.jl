# # Example 5 : Array Types
#
# A `ClimArray` (struct from `ClimateBase.jl`), which includes metadata read from nectdf file, is written back to a netcdf file via `NCTiles.jl`. Then the reverse is done to illustrate consistent workflows.

using ClimateBase, NCTiles, Unitful, NCDatasets, Dates

# File Paths
inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/")
if ~ispath(outputs); mkpath(outputs); end

savedir = joinpath(outputs,"ex5")
if ~ispath(savedir); mkpath(savedir); end

# ## 1. Read via NCTiles.jl and write via ClimateBase.jl

# ### Helper function

"""
        NCvar_to_ClimArray(ncvar::NCvar,ncdims::Dict)

Creates a ClimArray struct, from a NCvar struct, which can then be written
to a NetCDF file using `ClimateBase.ncwrite()`.
"""
function NCvar_to_ClimArray(ncvar::NCvar,ncdims::Dict)
    lons=ncdims["lon_c"].values[:]
    lats=ncdims["lat_c"].values[:]
    t=ncdims["tim"].values[:]
    dimensions = (Lon(lons), Lat(lats), ClimateBase.Ti(t))

    name=ncvar.name
    data=ncvar.values[:]
    units=ncvar.units
    return ClimArray(data, dimensions; name = name, attrib = Dict("units" => units))
end

# ### Main call sequence : 

fil=joinpath(outputs,"ex1/ex1_NetCDF.nc")
ncvars,ncdims,fileatts = NCTiles.readncfile(fil)

#filout=joinpath(savedir,"ex5_NCTiles.nc")
#write(ncvars["ETAN"],filout,globalattribs=fileatts)

filout=joinpath(savedir,"ex5_ClimateBase.nc")
A=NCvar_to_ClimArray(ncvars["ETAN"],ncdims)
ClimateBase.ncwrite(filout, A)

# ## 2. Read via ClimateBase.jl and write via NCTiles.jl
#
# ### Helper function

"""
        ClimArray_to_NCvar(C::ClimGrid,N::String)

Creates an NCvar struct, from a ClimArray struct, which can then be written
to a NetCDF file using `NCTiles.write()`.

```
C = ClimateBase.ncread(fil, "ETAN")
writefld = ClimArray_to_NCvar(C,"ETAN")
NCTiles.write(writefld,"myfile.nc")
```
"""
function ClimArray_to_NCvar(C::ClimArray,N::String)
        x = C.dims[1][:]
        y = C.dims[2][:]
        timevec = DateTime.(C.dims[3][:])
        timeunit = "days since 1992-01-01"
        timevec = NCDatasets.timeencode(timevec, timeunit)
        
        dims = [NCvar("lon","degrees_east",size(C.data)[1],x,Dict("long_name" => "longitude"),NCDatasets),
                NCvar("lat","degrees_north",size(C.data)[2],y,Dict("long_name" => "latitude"),NCDatasets),
                NCvar("time",timeunit,Inf,timevec,Dict(("long_name" => "Ti","standard_name" => "time")),NCDatasets)
                ]
        
        return NCvar(N,C.attrib["units"],dims,C.data,C.attrib,NCDatasets)
end   

# ## Main call sequence : 
#
# _Note: `ClimateBase.jl` relies on specific dimension names to identify space and time dimensions (as follows). This sometimes lead to warnings when files use different names._
#
# - latitude: lat, latitude, rlat, y, yc
# - longitude: lon, longitude, rlon, x, xc
# - time: time

ETAN = ClimateBase.ncread(fil, "ETAN")
ETAN = ClimArray_to_NCvar(ETAN,"ETAN")
filout=joinpath(savedir,"ex5_ClimArray_to_NCvar.nc")
NCTiles.write(ETAN,filout)

##

using MeshArrays

function ClimArray_to_MeshArray(C::ClimArray)
        u=uparse(C.attrib["units"])
        n=string(C.name)
        ln=C.attrib["long_name"]
        tim=DateTime.(collect(tmp.dims[3][:]))        
        m=varmeta(u,fill(0.5,3),tim,n,ln)
        #MeshArray(C.data;meta=m)

        nlon=length(tmp.dims[1][:])
        nlat=length(tmp.dims[2][:])
        XC = MeshArray(collect(tmp.dims[1][:])*ones(1,nlon))
        YC = MeshArray(ones(nlat,1)*collect(tmp.dims[2][:])')
        Γ = (XC=XC,YC=YC)

        MeshArray(C.data;meta=m),Γ
end

function NCvar_to_MeshArray(ncvar::NCvar)
        u=uparse(ncvar.units)
        n=string(ncvar.name)
        ln=ncvar.atts["long_name"]
        tim=ncvar.dims[3].values[:]
        m=varmeta(u,fill(0.5,3),tim,n,ln)
        #MeshArray(ncvar.values[:];meta=m)

        nlon=length(ncvar.dims[1].values[:])
        nlat=length(ncvar.dims[2].values[:])
        XC = MeshArray(ncvar.dims[1].values[:]*ones(1,nlon))
        YC = MeshArray(ones(nlat,1)*ncvar.dims[2].values[:]')
        Γ = (XC=XC,YC=YC)

        MeshArray(ncvar.values[:];meta=m),Γ
end
    
fil=joinpath(outputs,"ex1/ex1_NetCDF.nc")
tmp = ClimateBase.ncread(fil, "ETAN")
ETAN,Γ = ClimArray_to_MeshArray(tmp)

ncvars,ncdims,fileatts = NCTiles.readncfile(fil)
ETAN,Γ = NCvar_to_MeshArray(ncvars["ETAN"])
