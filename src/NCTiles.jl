module NCTiles

using NCDatasets,NetCDF

greet() = print("Hello World!")

include("writenctiles.jl")

export NCvar, Bindata, NCData, readbin
export createfile, addDim, addVar, addData, addDimData, readncfile

end # module
