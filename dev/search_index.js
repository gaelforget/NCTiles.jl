var documenterSearchIndex = {"docs":
[{"location":"#NCTiles.jl-Package-Documentation-1","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCTiles.jl reads and writes NetCDF files that represent e.g. subdivisions of Earth's surface (tiles). Inter-operability with popular climate model grids and MeshArrays.jl and generation of CF-compliant files are key goals of this package. ","category":"page"},{"location":"#Contents-1","page":"NCTiles.jl Package Documentation","title":"Contents","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCTiles.jl derives from the earlier nctiles implementation in gcmfaces (Forget et al. 2015).","category":"page"},{"location":"#Main-Features-1","page":"NCTiles.jl Package Documentation","title":"Main Features","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCTiles.jl first goes through lazy operations, on data structure, as it obtains information about variables etc. The second phase calls write function to instantiate and write files. Note: some of the included functions are interfaces to MITgcm output.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Data structures:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCvar contains information needed to write a NetCDF file which can include a list of filenames (see BinData) if the data is not loaded into memory.\nNCData contains a string (NetCDF file name) + metadata to read files.\nBinData is a container for one field.\nTileData contains a MeshArray or BinData struct in vals,   information about the tile layout in tileinfo, and metadata needed to   read/write tile data.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"As an example:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"struct TileData{T}\n    vals::T\n    tileinfo::Dict\n    tilesize::Tuple\n    precision::Type\n    numtiles::Int\nend","category":"page"},{"location":"#Use-examples-1","page":"NCTiles.jl Package Documentation","title":"Use examples","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"DataStructures/06_nctiles.ipynb in this GlobalOceanNotebooks repo provides a series of examples. Additional examples found in the examples/ folder include:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Example1.jl reads two-dimensional fields on a regular grid (\"lat-lon\") read from binary files, and then writes them to a netcdf file. This example illustrates the use of either NCDatasets.jl or NetCDF.jl as the backend.\nExample2.jl reads two-dimensional fields from the netcdf file generated in Example1, and then re-writes them to a new netcdf file.\nex_3.jl is an example of interpolated model output processing in CBIOMES where several variables are included in the same NetCDF/NCTiles file.\nex_4.jl generates a tiled NetCDF output (i.e., a nctiles output) for a global 2D field on the non-regular LLC90 grid (see MeshArrays.jl). Since the tile width is set to 90, this creates 13 files.\nex_5.jl shows how to write a ClimGrid struct from the ClimateTools package to a NetCDF/NCTiles file using NCTiles.","category":"page"},{"location":"#Using-NCTiles-1","page":"NCTiles.jl Package Documentation","title":"Using NCTiles","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The core functionality of NCTiles comes from a series of data structures that contain the information needed write to NetCDF files. This includes the information and methods needed to read from source files. The data structure used for writing a variable is NCvar, which includes that variable's data and metadata. The data itself can be in memory and included directly in the NCvar struct, or can be described in another class of data structures, with names ending in Data. These included BinData, for data in binary files, NCData, for data in NetCDF files, and TileData for data to be written out in chunks.","category":"page"},{"location":"#Basic-Example-1","page":"NCTiles.jl Package Documentation","title":"Basic Example","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Here we show how to write a NetCDF file from a series of Binary data files.","category":"page"},{"location":"#Define-Dimensions-1","page":"NCTiles.jl Package Documentation","title":"Define Dimensions","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The first step for creating a NetCDF file is to define your dimensions. Each dimension is specified by an NCvar. Dimensions should be in an Array in the order corresponding to your variable data (if your data dimensions are lon x lat x time, dimensions should be in that order as well). In this example we have a regular half-degree lat-lon grid with 10 time steps (as in ex_1.jl). This is how we define the dimensions:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"lon = -179.75:0.5:179.75\nlat = -89.75:0.5:89.75\ntime = 1:10\n\ndims = [NCvar(\"lon\",\"degrees_east\",size(lon),lon,Dict(\"long_name\" => \"longitude\"),NCDatasets),\n        NCvar(\"lat\",\"degrees_north\",size(lat),lat,Dict(\"long_name\" => \"latitude\"),NCDatasets),\n        NCvar(\"time\",\"days since 1992-01-01\",Inf,time,Dict((\"long_name\" => \"tim\",\"standard_name\" => \"time\")),NCDatasets)\n        ]","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Let's go through the NCvar constructor. Here is the struct definition for reference:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"struct NCvar\n    name::String\n    units::String\n    dims\n    values\n    atts::Union{Dict,Nothing}\n    backend::Module\nend","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The first attribute, name, should be a String and is what you want to call the variable in the file. The second are the units, which should also be a String. We then specify the dimensions, dims. For Dimension variables dims should be of length 1 (calling size on your dimension values like above if sufficient). Next you specify the actual dimension values. For a Dimension variable, this must be a 1 dimensional array, like above. After the values you can specify any additional attributes that you want to add to the variable as a dictionary. The last attribute is the backend, which allows you to choose between NCDatasets.jl and NetCDF.jl. We have some support for NetCDF.jl and full support for NCDatasets.jl. Note that in creating these NCvar structs we do not do any CF Compliance checks, it is the user's responsibility to provide CF-compliant units.","category":"page"},{"location":"#Define-the-Data-Source-1","page":"NCTiles.jl Package Documentation","title":"Define the Data Source","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Once you've created the dimensions for your NetCDF file you can create NCvar for your variable. Here we are going to create one pointing to data that is stored in multiple Binary files, one for each time step. The first step is to create this pointer to the data, which is the BinData struct. For example:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"precision = Float32\ndatapath = \"data/binfiles/\"\nfnames = joinpath.(Ref(datapath),\n                    [\"Chl050.001.data\",\"Chl050.002.data\",\"Chl050.003.data\",\"Chl050.004.data\",\"Chl050.005.data\",\"Chl050.006.data\",\"Chl050.007.data\",\"Chl050.008.data\",\"Chl050.009.data\",\"Chl050.010.data\"])\nvardata = BinData(fnames,precision,(length(lon),length(lat)))\n# or: vardata = BinData(fnames,precision,(length(lon),length(lat)),1)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"And for reference, the struct definition for BinData:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"struct BinData # Pointer to data stored in binary files- contains info needed to read in\n    fnames::Union{Array{String},String}\n    precision::Type\n    iosize::Tuple\n    fldidx::Int\nend","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"In order to read data from a binary file, we need to know where the files are and their names, the precision that the data is written in, and the dimensions of the data. The first argument, fnames, should be a single file path String or an Array of file paths Strings. The second should be the precision that the data is written in the file, here our data is Float32. Precision should be a Type. Finally we need to know the size of the data that we are reading from the file, this should be specified as a Tuple. If we have multiple variables written in the same file, we can additionally specify the index of that variable, say if it's the 10th variable in the file. In this example there is only one variable in the file, so we can specify 1 or leave it out and it will be assumed to be 1.","category":"page"},{"location":"#Create-the-NCvar-1","page":"NCTiles.jl Package Documentation","title":"Create the NCvar","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Now we can create the NCvar for the variable we want to write to the file.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"varname = \"Chl050\"\nunits = \"mg Chl\"\nlongname = \"Average chlorophyll concentration (top 50m)\"\nmyvar = NCvar(varname,units,dims,vardata,Dict(\"long_name\" => longname),NCDatasets)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Creating the final NCvar for our variable is similar to creating a dimension NCvar. We specify the name we want to use in the file and the units. Here we use the dims array and the vardata struct we created above. We specify a long_name attribute, and finally indicate that we want to use NCDatasets in the backend.","category":"page"},{"location":"#Writing-to-the-NetCDF-File-1","page":"NCTiles.jl Package Documentation","title":"Writing to the NetCDF File","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Assuming you've created the above structs as expected, executing the write function is as simple as:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"README = \"A useful README that describes the data in the file.\"\nattributes = Dict([\"_FillValue\"=>NaN, \"missing_value\" => NaN])\nwrite(myvar,\"data/mydata.nc\",README=README,globalattribs=attributes)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The write function requires at a minimum an NCvar and the output file path. It writes the NCvar to the file with default global attributes. Additionally you can specify a README and global attributes, by passing a String or Array of Strings to the README keyword argument or by providing a Dict to the globalattribs keyword argument, as shown above.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"If you would like to write multiple variables to the same file, you can pass a Dict{String,NCvar} into the write function:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"myvars = Dict([\"myvar1\" => myvar1,\n                \"myvar2\" => myvar2])\nwrite(myvar,\"data/mydata.nc\")","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Where the keys of the Dict should match the name attributes of the NCvar struct values.","category":"page"},{"location":"#Other-Data-Structures-1","page":"NCTiles.jl Package Documentation","title":"Other Data Structures","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"In the example above we wrote a NetCDF file with data sourced from Binary Files, specified by the BinData struct. We have a few other structs for different kinds of data:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCData: for data sourced from a NetCDF file\nTileData: for data to be written into separate tile files","category":"page"},{"location":"#NCData-1","page":"NCTiles.jl Package Documentation","title":"NCData","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"NCData structs contain the necessary information to read data from a NetCDF file.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"struct NCData\n    fname::AbstractString\n    varname::AbstractString\n    backend::Module\n    precision::Type\nend","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"For example, if you wanted to use the NetCDF file created before as a data source, you would use the NCData constructor:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"myvardata = NCData(\"data/mydata.nc\",\"Chl050\",NCDatasets)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Where the arguments are: file path; variable name; backend.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Alternatively, we provide the function readncfile which creates NCvars containing the NCData structs for all the variables in the file:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"ncvars,ncdims,fileatts = readncfile(\"data/mydata.nc\")","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Here, the ncvars dictionary contains NCvars of all the variables in the file. Each NCvar has NCData structs in the values attribute, which avoids reading in all the data from the file. In this case the NCData can be accessed as myvardata = ncvars[\"Chl050\"].values.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"To re-write this exact file run:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"write(ncvars,joinpath(\"data/mydata2.nc\"),globalattribs=fileatts)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"You can see this process demonstrated in ex_2.jl.","category":"page"},{"location":"#TileData-1","page":"NCTiles.jl Package Documentation","title":"TileData","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The TileData struct is used to chunk up data and write to separate files. We do this using the MeshArrays package. This is demonstrated in more detail in ex_4.jl. First, specify your grid and read in the grid variables:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"grid = GridSpec(\"LatLonCap\",\"grids/GRID_LLC90/\")\ngridvars = GridLoad(grid)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Where GridSpec() and GridLoad() are from the MeshArrays package (you can refer to the MeshArrays documentation for more information about these functions and grids).","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The next step is to specify the tile, or chunk, size as a tuple. The data is chunked in the horizontal dimension, so tile sizes should be two dimensional tuple. If the data is three dimensional, say its full dimension is NxMx10 and the tile size is nxm, the chunks will be nxmx10. Here we set the tile size to 90x90:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"tilesize = (90,90)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"When defining dimensions for TileData variables, the horizontal dimensions should be the size of the tiles, and their values integers 1:n or 1:m for an nxm tile:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"time = 1:10\ndep = gridvars[\"RC\"]\ndims = [\n    NCvar(\"i_c\",\"1\",tilesize[1],1:tilesize[1],Dict(\"long_name\" => \"Cartesian coordinate 1\"),NCDatasets),\n    NCvar(\"j_c\",\"1\",tilesize[2],1:tilesize[2],Dict(\"long_name\" => \"Cartesian coordinate 2\"),NCDatasets),\n    NCvar(\"dep_c\",\"m\",size(dep),dep,Dict(\"long_name\" => \"depth\",\"standard_name\" => \"depth\",\"positive\" => \"down\"),NCDatasets),\n    NCvar(\"time\",\"days since 1992-01-01\",Inf,time,Dict((\"long_name\" => \"tim\",\"standard_name\" => \"time\")),NCDatasets)\n]","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The latitude and longitude variables will be written to the file separately, their data specified by TileData structs:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"tillat = TileData(gridvars[\"YC\"],tilesize,grid)\nvarlat = NCvar(\"lat\",\"degrees_north\",dims[1:2],tillat,Dict(\"long_name\" => \"latitude\"),NCDatasets)\ntillon = TileData(gridvars[\"XC\"],tilesize,grid)\nvarlon = NCvar(\"lon\",\"degrees_east\",dims[1:2],tillon,Dict(\"long_name\" => \"longitude\"),NCDatasets)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Since the data for latitude and longitude are held in memory (in the gridvars dictionary), we can specify it directly. At construction, the TileData struct will create the mapping for which indices of the latitude and longitude data should be put in each tile. When a file is written, NCTiles will use this mapping to extract the chunk for that file. The dimensions for the corresponding NCvars should have the dimensions dims[1:2], corresponding to i_c and j_c.","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The variable we want to write is in a binary data file, so we can use a BinData struct in the TileData for our variable:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"vardata = TileData(BinData(fnames,prec,iosize),\n                    tilesize,\n                    grid)\nmyvar = NCvar(varname,\"myunits\",dims,vardata,Dict(),NCDatasets)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The final step is to create the NCvars and write them to the NetCDF files:","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"vars = Dict([varname => myvar,\n                    \"lon\" => varlon,\n                    \"lat\" => varlat\n            ])\n\nsavenamebase = \"data/mytiledata\"\nwrite(vars,savenamebase)","category":"page"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"The write function will create one file for each tile, using savenamebase as a prefix for the file path. It will append a zero-padded number to the end of the filename, along with the extension .nc. For this example we would have the files data/mytiledata.0001.nc, data/mytiledata.0002.nc, ..., data/mytiledata.0013.nc.","category":"page"},{"location":"#Index-1","page":"NCTiles.jl Package Documentation","title":"Index","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"","category":"page"},{"location":"#API-/-Functions-1","page":"NCTiles.jl Package Documentation","title":"API / Functions","text":"","category":"section"},{"location":"#","page":"NCTiles.jl Package Documentation","title":"NCTiles.jl Package Documentation","text":"Modules = [NCTiles]\nOrder   = [:type,:function]\nPrivate = false","category":"page"},{"location":"#NCTiles.BinData","page":"NCTiles.jl Package Documentation","title":"NCTiles.BinData","text":"BinData\n\nData structure containing a string or an array of strings (NetCDF     file names) as well as metadata needed to read a file.\n\n\n\n\n\n","category":"type"},{"location":"#NCTiles.BinData-Tuple{Union{String, Array{String,N} where N},Type,Tuple}","page":"NCTiles.jl Package Documentation","title":"NCTiles.BinData","text":"BinData(fnames::Union{Array{String},String},precision::Type,iosize::Tuple)\n\nConstruct a BinData struct for files that contain one field.\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.NCData","page":"NCTiles.jl Package Documentation","title":"NCTiles.NCData","text":"NCData\n\nData structure containing a string or an array of strings (file names) of     NetCDF files as well as information needed to read a file.\n\n\n\n\n\n","category":"type"},{"location":"#NCTiles.NCvar","page":"NCTiles.jl Package Documentation","title":"NCTiles.NCvar","text":"NCvar\n\nData structure containing information needed to write a NetCDF file. This includes a list of filenames (see Bindata) if the data is not loaded into memory.\n\n\n\n\n\n","category":"type"},{"location":"#NCTiles.TileData","page":"NCTiles.jl Package Documentation","title":"NCTiles.TileData","text":"TileData{T}\n\nData structure containing either a MeshArray struct or BinData struct (see vals),     MeshArray structs describing the tile layout (tileinfo), and other information for     reading/writing tile data.\n\n\n\n\n\n","category":"type"},{"location":"#NCTiles.TileData-Tuple{Any,Tuple,MeshArrays.gcmgrid}","page":"NCTiles.jl Package Documentation","title":"NCTiles.TileData","text":"TileData(vals,tilesize::Tuple)\n\nConstruct a TileData struct. First generate the tileinfo, precision, and numtiles attributes.\n\n\n\n\n\n","category":"method"},{"location":"#Base.write-Tuple{Dict{AbstractString,NCvar},String}","page":"NCTiles.jl Package Documentation","title":"Base.write","text":"write(myflds::Dict,savename::String;README=\"\",globalattribs=Dict())\n\nCreates NetCDF file and writes myflds and all their dimensions to the file.\n\n\n\n\n\n","category":"method"},{"location":"#Base.write-Tuple{NCvar,String}","page":"NCTiles.jl Package Documentation","title":"Base.write","text":"write(myfld::NCvar,savename::String;README=\"\",globalattribs=Dict())\n\nCreates NetCDF file and writes myfld and all its dimensions to the file.\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addData-Tuple{Union{Array, NCDatasets.CFVariable, NetCDF.NcVar},NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addData","text":"addData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar},var::NCvar)\n\nFill variable with data in netcdf file using either NCDatasets.jl or NetCDF.jl\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addDim-Tuple{NCDatasets.NCDataset,NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addDim","text":"addDim(ds::NCDatasets.Dataset,dimvar::NCvar) # NCDatasets\n\nAdd a dimension to a NCDatasets.Dataset\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addDim-Tuple{NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addDim","text":"addDim(dimvar::NCvar)\n\nAdd a dimension to a NetCDF file using NetCDF.jl\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addDimData-Tuple{Any,NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addDimData","text":"addDimData(v::Union{NCDatasets.CFVariable,NetCDF.NcVar,Array},var::NCvar)\n\nAdd dimension data to predefined dimensions in a NetCDF file.\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addVar-Tuple{NCDatasets.NCDataset,NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addVar","text":"addVar(ds::NCDatasets.Dataset,field::NCvar)\n\nAdd a variable to a NetCDF file using NCDatasets.jl\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addVar-Tuple{NCvar,Array{NetCDF.NcDim,N} where N}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addVar","text":"addVar(field::NCvar,dimlist::Array{NetCDF.NcDim})\n\nAdd a variable with dimensions dimlist to a NetCDF file using NetCDF.jl\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.addVar-Tuple{NCvar}","page":"NCTiles.jl Package Documentation","title":"NCTiles.addVar","text":"addVar(field::NCvar})\n\nAdd a variable and its dimensions to a NetCDF file using NetCDF.jl\n\n\n\n\n\n","category":"method"},{"location":"#NCTiles.createfile","page":"NCTiles.jl Package Documentation","title":"NCTiles.createfile","text":"createfile(filename, field::Union{NCvar,Dict{String,NCvar}}, README;\n            fillval=NaN, missval=NaN, itile=1, ntile=1)\n\nCreate NetCDF file and add variable + dimension definitions using either NCDatasets.jl or NetCDF.jl\n\n\n\n\n\n","category":"function"},{"location":"#NCTiles.readbin","page":"NCTiles.jl Package Documentation","title":"NCTiles.readbin","text":"readbin(fname::String,prec::Type,iosize::Tuple,fldidx=1)\n\nRead in a binary file to an Array.\n\n\n\n\n\n","category":"function"},{"location":"#NCTiles.readbin","page":"NCTiles.jl Package Documentation","title":"NCTiles.readbin","text":"readbin(flddata::BinData,tidx=1)\n\nRead in a binary file as an array as specified by BinData\n\n\n\n\n\n","category":"function"},{"location":"#NCTiles.readncfile","page":"NCTiles.jl Package Documentation","title":"NCTiles.readncfile","text":"readncfile(fname,backend::Module=NCDatasets)\n\nRead in a NetCDF file and return variables/dimensions as NCvar structs, and     file attributes as Dict. Large variables/dimensions are not loaded into     memory. This can use either NCDatasets.jl or NetCDF.jl\n\n\n\n\n\n","category":"function"}]
}
