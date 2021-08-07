# User Guide

`NCTiles.jl` first goes through lazy operations, on data structure, as it obtains information about variables etc. The second phase calls `write` function to instantiate and write files. _Note:_ some of the included functions are interfaces to `MITgcm` output.

Data structures:

- `NCvar` contains information needed to write a NetCDF file which can include a list of filenames (see `BinData`) if the data is not loaded into memory.
- `NCData` contains a string (NetCDF file name) + metadata to read files.
- `BinData` is a container for one field.
- `TileData` contains a `MeshArray` or `BinData` struct in `vals`,
    information about the tile layout in `tileinfo`, and metadata needed to
    read/write tile data.

As an example:

```julia
struct TileData{T}
    vals::T
    tileinfo::Dict
    tilesize::Tuple
    precision::Type
    numtiles::Int
end
```

## Data Structures & Functions

Higher-level APIs, which are practical for automated or distributed workflows which can be called upon e.g. as a model runs forward in time, are readily documented in the examples. Here we take a more detailed look at a basic one to document the underlying / core data structures and functionalities.

The core functionality of NCTiles comes from a series of data structures that contain the information needed write to NetCDF files. This includes the information and methods needed to read from source files. The data structure used for writing a variable is `NCvar`, which includes that variable's data and metadata. The data itself can be in memory and included directly in the `NCvar` struct, or can be described in another data structure (`NCData` / `BinData` / `TileData`) :

- `BinData` for data in binary files
- `NCData` for data in NetCDF files
- `TileData` for e.g. tiled model output when subdomains is often written out in distributed fashion across file collections (see [MeshArrays](https://juliaclimate.github.io/MeshArrays.jl/dev) for suitable Earth domain decomposition examples).

## Basic Example

Here we show how to write a metadata-rich `NetCDF` file from a series of binary data files, which represents output from a climate model (`MITgcm` output in this example). We try here to document the `metadata` specification in detail as one of the main goals of `NCTiles.jl` is to facilitate the production of metadata-rich data sets that are easily reuseable. 

### Define Dimensions

The first step for creating a NetCDF file is to define your dimensions. Each dimension is specified by an `NCvar`. Dimensions should be in an `Array` in the order corresponding to your variable data (if your data dimensions are lon x lat x time, dimensions should be in that order as well). In this example we have a regular half-degree lat-lon grid with 10 time steps (as in `ex_1.jl`). This is how we define the dimensions:

```julia
lon = -179.75:0.5:179.75
lat = -89.75:0.5:89.75
time = 1:10

dims = [NCvar("lon","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("time","days since 1992-01-01",Inf,time,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]
```

Let's go through the `NCvar` constructor. Here is the struct definition for reference:

```julia
struct NCvar
    name::String
    units::String
    dims
    values
    atts::Union{Dict,Nothing}
    backend::Module
end
```

The first attribute, `name`, should be a `String` and is what you want to call the variable in the file. The second are the units, which should also be a `String`. We then specify the dimensions, `dims`. For Dimension variables `dims` should be of length 1 (calling `size` on your dimension values like above if sufficient). Next you specify the actual dimension values. For a Dimension variable, this must be a 1 dimensional array, like above. After the values you can specify any additional attributes that you want to add to the variable as a dictionary. The last attribute is the backend, which allows you to choose between `NCDatasets.jl` and `NetCDF.jl`. We have some support for `NetCDF.jl` and full support for `NCDatasets.jl`. Note that in creating these `NCvar` structs we do not do any CF Compliance checks, it is the user's responsibility to provide CF-compliant units.

### Define the Data Source

Once you've created the dimensions for your NetCDF file you can create `NCvar` for your variable. Here we are going to create one pointing to data that is stored in multiple Binary files, one for each time step. The first step is to create this pointer to the data, which is the `BinData` struct. For example:

```julia
precision = Float32
datapath = "data/binfiles/"
fnames = joinpath.(Ref(datapath),
                    ["Chl050.001.data","Chl050.002.data","Chl050.003.data","Chl050.004.data","Chl050.005.data","Chl050.006.data","Chl050.007.data","Chl050.008.data","Chl050.009.data","Chl050.010.data"])
vardata = BinData(fnames,precision,(length(lon),length(lat)))
# or: vardata = BinData(fnames,precision,(length(lon),length(lat)),1)
```

And for reference, the struct definition for `BinData`:

```julia
struct BinData # Pointer to data stored in binary files- contains info needed to read in
    fnames::Union{Array{String},String}
    precision::Type
    iosize::Tuple
    fldidx::Int
end
```

In order to read data from a binary file, we need to know where the files are and their names, the precision that the data is written in, and the dimensions of the data. The first argument, `fnames`, should be a single file path `String` or an `Array` of file paths `String`s. The second should be the precision that the data is written in the file, here our data is `Float32`. Precision should be a `Type`. Finally we need to know the size of the data that we are reading from the file, this should be specified as a `Tuple`. If we have multiple variables written in the same file, we can additionally specify the index of that variable, say if it's the 10th variable in the file. In this example there is only one variable in the file, so we can specify 1 or leave it out and it will be assumed to be 1.

### Create the NCvar

Now we can create the NCvar for the variable we want to write to the file.

```julia
varname = "Chl050"
units = "mg Chl"
longname = "Average chlorophyll concentration (top 50m)"
myvar = NCvar(varname,units,dims,vardata,Dict("long_name" => longname),NCDatasets)
```

Creating the final `NCvar` for our variable is similar to creating a dimension `NCvar`. We specify the name we want to use in the file and the units. Here we use the `dims` array and the `vardata` struct we created above. We specify a `long_name` attribute, and finally indicate that we want to use `NCDatasets` in the backend.

### Writing to the NetCDF File

Assuming you've created the above structs as expected, executing the `write` function is as simple as:

```julia
README = "A useful README that describes the data in the file."
attributes = Dict(["_FillValue"=>NaN, "missing_value" => NaN])
write(myvar,"data/mydata.nc",README=README,globalattribs=attributes)
```

The `write` function requires at a minimum an `NCvar` and the output file path. It writes the `NCvar` to the file with default global attributes. Additionally you can specify a `README` and global attributes, by passing a `String` or Array of Strings to the `README` keyword argument or by providing a `Dict` to the `globalattribs` keyword argument, as shown above.

If you would like to write multiple variables to the same file, you can pass a `Dict{String,NCvar}` into the `write` function:

```julia
myvars = Dict(["myvar1" => myvar1,
                "myvar2" => myvar2])
write(myvar,"data/mydata.nc")
```

Where the keys of the `Dict` should match the `name` attributes of the `NCvar` struct values.

## Other Data Structures

In the example above we wrote a NetCDF file with data sourced from Binary Files, specified by the `BinData` struct. We have a few other structs for different kinds of data:

- `NCData`: for data sourced from a NetCDF file
- `TileData`: for data to be written into separate tile files

#### NCData

`NCData` structs contain the necessary information to read data from a NetCDF file.

```julia
struct NCData
    fname::AbstractString
    varname::AbstractString
    backend::Module
    precision::Type
end
```

For example, if you wanted to use the NetCDF file created before as a data source, you would use the `NCData` constructor:

```julia
myvardata = NCData("data/mydata.nc","Chl050",NCDatasets)
```

Where the arguments are: file path; variable name; backend.

Alternatively, we provide the function `readncfile` which creates `NCvar`s containing the `NCData` structs for all the variables in the file:

```julia
ncvars,ncdims,fileatts = readncfile("data/mydata.nc")
```

Here, the `ncvars` dictionary contains `NCvar`s of all the variables in the file. Each `NCvar` has `NCData` structs in the `values` attribute, which avoids reading in all the data from the file. In this case the `NCData` can be accessed as `myvardata = ncvars["Chl050"].values`.

To re-write this exact file run:

```julia
write(ncvars,joinpath("data/mydata2.nc"),globalattribs=fileatts)
```

You can see this process demonstrated in `ex_2.jl`.

### TileData

The `TileData` struct is used to chunk up data and write to separate files. We do this using the `MeshArrays` package. This is demonstrated in more detail in `ex_4.jl`. First, specify your grid and read in the grid variables:

```julia
grid = GridSpec("LatLonCap","grids/GRID_LLC90/")
gridvars = GridLoad(grid)
```

Where `GridSpec()` and `GridLoad()` are from the `MeshArrays` package (you can refer to the `MeshArrays` documentation for more information about these functions and grids).

The next step is to specify the tile, or chunk, size as a `tuple`. The data is chunked in the horizontal dimension, so tile sizes should be two dimensional `tuple`. If the data is three dimensional, say its full dimension is `NxMx10` and the tile size is `nxm`, the chunks will be `nxmx10`. Here we set the tile size to `90x90`:

```julia
tilesize = (90,90)
```

When defining dimensions for `TileData` variables, the horizontal dimensions should be the size of the tiles, and their values integers `1:n` or `1:m` for an `nxm` tile:

```julia
time = 1:10
dep = gridvars["RC"]
dims = [
    NCvar("i_c","1",tilesize[1],1:tilesize[1],Dict("long_name" => "Cartesian coordinate 1"),NCDatasets),
    NCvar("j_c","1",tilesize[2],1:tilesize[2],Dict("long_name" => "Cartesian coordinate 2"),NCDatasets),
    NCvar("dep_c","m",size(dep),dep,Dict("long_name" => "depth","standard_name" => "depth","positive" => "down"),NCDatasets),
    NCvar("time","days since 1992-01-01",Inf,time,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
]
```

The latitude and longitude variables will be written to the file separately, their data specified by `TileData` structs:

```julia
tillat = TileData(gridvars["YC"],tilesize,grid)
varlat = NCvar("lat","degrees_north",dims[1:2],tillat,Dict("long_name" => "latitude"),NCDatasets)
tillon = TileData(gridvars["XC"],tilesize,grid)
varlon = NCvar("lon","degrees_east",dims[1:2],tillon,Dict("long_name" => "longitude"),NCDatasets)
```

Since the data for latitude and longitude are held in memory (in the `gridvars` dictionary), we can specify it directly. At construction, the TileData struct will create the mapping for which indices of the latitude and longitude data should be put in each tile. When a file is written, `NCTiles` will use this mapping to extract the chunk for that file. The dimensions for the corresponding `NCvar`s should have the dimensions `dims[1:2]`, corresponding to `i_c` and `j_c`.

The variable we want to write is in a binary data file, so we can use a `BinData` struct in the `TileData` for our variable:

```julia
vardata = TileData(BinData(fnames,prec,iosize),
                    tilesize,
                    grid)
myvar = NCvar(varname,"myunits",dims,vardata,Dict(),NCDatasets)
```

The final step is to create the `NCvar`s and write them to the `NetCDF` files:

```julia
vars = Dict([varname => myvar,
                    "lon" => varlon,
                    "lat" => varlat
            ])

savenamebase = "data/mytiledata"
write(vars,savenamebase)
```

The `write` function will create one file for each tile, using `savenamebase` as a prefix for the file path. It will append a zero-padded number to the end of the filename, along with the extension `.nc`. For this example we would have the files `data/mytiledata.0001.nc`, `data/mytiledata.0002.nc`, ..., `data/mytiledata.0013.nc`.

