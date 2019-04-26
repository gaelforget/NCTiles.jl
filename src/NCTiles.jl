module NCTiles

using NetCDF, MeshArrays, Printf, NCDatasets

include("utilities.jl")
export read_nctiles

include("writenctiles.jl")

export NCvar, Bindata, NCData, readbin
export createfile, addDim, addVar, addData, addDimData, readncfile

end # module
