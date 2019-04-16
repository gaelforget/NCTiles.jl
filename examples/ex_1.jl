using NCTiles

interpdir = "../../../DarwinModelOutputSamples/ptr/diags_interp_20190326_1256/"

selectfields = ["Chl050","Plk050"]

fnames = interpdir*selectfields[1]*'/'.*filter(x -> occursin(".data",x), readdir(interpdir*selectfields[1]))

prec = Float32
lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
n1,n2 = (length(lon),length(lat))
nsteps = 2
timeinterval = 3
time_steps = timeinterval/2:timeinterval:timeinterval*nsteps
longname = "Average chlorophyll concentration (top 50m)"
units = "mg Chl"

README = readlines("README")

# Set up data for writing to NCfiles
dims = [NCvar("lon_c","longitude","degrees east",size(lon),lon,Dict()),
        NCvar("lat_c","latitude","degrees north",size(lat),lat,Dict()),
        NCvar("time","time","days from 1992-01-01",Inf,time_steps,Dict(("standard_name" => "time")))
        ]

fielddata = Bindata(fnames,prec,(n1,n2))

field = NCvar(selectfields[1],longname,units,dims,fielddata,Dict())


# Try out NCDatasets
# Only supports 1D dims??
using NCDatasets

ds = createfile("test.nc",field.name,README)

addDim.(Ref(ds),field.dims)

atts = merge(Dict(("long_name" => field.longname,"units" =>field.units)),field.atts)
fieldvar = defVar(ds,field.name,field.values.precision,(field.dims[1].name,field.dims[2].name,field.dims[3].name),attrib=atts)

addData(fieldvar,field)
addDimData.(Ref(ds),field.dims)

close(ds)
