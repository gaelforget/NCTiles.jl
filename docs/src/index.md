# `NCTiles.jl` Package Documentation

**NCTiles.jl** reads and writes [NetCDF files](https://en.wikipedia.org/wiki/NetCDF) that represent e.g. subdivisions of Earth's surface (`tiles`). Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and generation of [CF-compliant](http://cfconventions.org) files are key goals of this package. 

## Contents

```@contents
```

_`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015))._

## Main Features

`NCTiles.jl` first goes through lazy operations, on data stucture, as it obtains information about variables etc. The second phase calls `createfile` and then the `add*` functions to instantiate files. _Note:_ some of the included functions are interfaces to `MITgcm` output.

Data structures:

- `NCData` contains a string or an array of strings (NetCDF file names) + metadata to read files. 
- `NCvar` contains information needed to write a NetCDF file which can include a list of filenames (see `Bindata`) if the data is not loaded into memory.
- `BinData` is a container for one field.
- `TileData ` contains a `MeshArray` or `BinData` struct in `vals`,
    information about the tile layout in `tileinfo`, and metadata needed for
    reading/writing tile data.

As an example:
    
```
struct TileData{T}
    vals::T
    tileinfo::Dict
    tilesize::Tuple
    precision::Type
    numtiles::Int
end
```    

## Use examples

`DataStructures/06_nctiles.ipynb` in this [GlobalOceanNotebooks repo](https://github.com/gaelforget/GlobalOceanNotebooks/) provides a series of examples that overlap somewhat with the ones found in `Examples/ex*.jl`:

- `ex_1.jl` reads a `binary` file containing one interpolated 2D field on a regular grid. It then writes that array to a `NetCDF`/`NCTiles` file.
- `ex_2.jl` reads data from a `NetCDF` file containing one `tile` of model output. It then writes it to a new `NetCDF`/`NCTiles` file. This uses 3D data on a non-regular grid for one ocean subdivision (`tile`).
- `ex_3.jl` is an example of interpolated model output processing in `CBIOMES` where several variables are included in the same `NetCDF`/`NCTiles` file.
- `ex_4.jl` generates a tiled netcdf output (i.e., a `nctiles` output) for a global 2D field on the non-regular `LLC90` grid (see `MeshArrays.jl`). Since the tile width is set to 90, this creates 13 files.

## Using NCTiles

The core functionality of NCTiles comes from a series of data structures that contain the information needed write NetCDF files, including the information need to read from a variety of types of source files. The main data structure used for writing is `NCvar`, which includes all the information needed to write a variable to a NetCDF file. The data itself can be in memory and included directly in the `NCvar` struct, or can be described in another class of data structures, with names ending in `Data`. These included `BinData`, for data in binary files, `NCData`, for data in NetCDF files, and `TileData` for data to be written out in separate tiles.

### Basic Example

In this basic example we will show how to write a NetCDF file from a series of Binary data files.

#### Define Dimensions

The first step for creating a NetCDF file would be to define your dimensions. Each dimension is specified by an `NCvar`. Dimensions should be in an `Array` in the order corresponding to your variable data (if your data dimensions are lon x lat x time, dimensions should be in that order as well). In this example we have a square half-degree lat-lon grid with 10 time steps (this example is based loosly on `ex_1.jl` above). We define the dimensions like so:

```julia
lon = -179.75:0.5:179.75
lat = -89.75:0.5:89.75
time = 1:10

dims = [NCvar("lon","degrees_east",size(lon),lon,Dict("long_name" => "longitude"),NCDatasets),
        NCvar("lat","degrees_north",size(lat),lat,Dict("long_name" => "latitude"),NCDatasets),
        NCvar("time","days since 1992-01-01",Inf,time,Dict(("long_name" => "tim","standard_name" => "time")),NCDatasets)
        ]
```

Let's go through the NCvar constructor. Here is the struct definition for reference:

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

The first attribute, `name`, should be a `String` and is what you want to call the variable in the file. The second are the units, which should also be a `String`. We then specify the dimensions, `dims`. For Dimension variables `dims` should be of length 1 (calling `size` on your dimension values like above if sufficnt). Next you specify the actual dimension values. For a Dimension variable, this must be a 1 dimensional array, like above. After the values you can specify any additional attributes that you want to add to the variable as a dictionary. The last attribute is the backend, which allows you to choose between `NCDatasets.jl` and `NetCDF.jl`. We have some support for `NetCDF.jl` and full support for `NCDatsets.jl`. Note that in creating these `NCvar` structs we do not do any CF Compliance checks, it is up to you to give CF compliant units and attributes.

#### Define the Data Source

Once you've created the dimensions for your NetCDF file you can create `NCvar` for your variable. Here we are going to create one pointing to data that is stored in multiple Binary files, one for each time step. The first step is to create this pointer to the data, which is the `BinData` struct:

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

#### Create the NCvar

Now we can create the NCvar for the variable we want to write to the file.

```julia
varname = "Chl050"
units = "mg Chl"
longname = "Average chlorophyll concentration (top 50m)"
myvar = NCvar(varname,units,dims,vardata,Dict("long_name" => longname),NCDatasets)
```

Creating the final `NCvar` for our variable is similar to creating the Dimension `NCvar`s. We specify the name we want to use in the file and the units. Here we use the `dims` array we created above, followed by the `vardata` `BinData` struct we created. We specify a long_name attribute, and finally indicate that we want to use NCDatasets in the backend.

#### Writing to the NetCDF File

Assuming you've created the above structs properly, writing is fairly straightforward:

```julia
README = "A useful README that describes the data in the file."
attributes = Dict(["_FillValue"=>NaN, "missing_value" => NaN])
write(myvar,"data/mydata.nc",README=README,globalattribs=attributes)
```

The `write` function requires at a minimum an `NCvar` and the path of a file to write to. It will write the `NCvar` to the file with default global attributes. Additionally you can specify a `README` and global attributes, by passing a `String` or Array of Strings to the `README` keyword argument or by providing a `Dict` to the `globalattribs` keyword argument, as shown above.

If you would like to write multiple variables to the same file, you can pass a `Dict{String,NCvar}` into the `write` function:

```julia
myvars = Dict(["myvar1" => myvar1,
                "myvar2" => myvar2])
write(myvar,"data/mydata.nc")
```

Where the keys of the `Dict` should match the `name` attributes of the `NCvar` struct values.

### Other Data Structures

In the example above we wrote a NetCDF file with data sourced from Binary Files, specified by the `BinData` struct. We have a few other structs for different kinds of data:

- `NCData`: for data sourced from a NetCDF file
- `TileData`: for data to be written into separated tile files

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

For example, if I wanted to use the NetCDF file I wrote above as a data source, I'd create:

```julia
myvardata = NCData("data/mydata.nc","Chl050",NCDatasets)
```

Where I specify the name of the file, the name of the variable in the file, and the backend.

Alternatively, we proved the function `readncfile` which creates `NCvar`s containing the `NCData` structs for all the variables in the file:

```julia
ncvars,ncdims,fileatts = readncfile("data/mydata.nc")
myvardata = ncvars["Chl050"].values
```

With the command above, ncvars will be a dictionary containing `NCvar`s of all the variables in the file, containing `NCData` structs rather than reading in all the data in the file.

#### TileData

## Index

```@index
```

## API / Functions

```@autodocs
Modules = [NCTiles]
Order   = [:type,:function]
Private = false
```
