module NCTiles

using NCDatasets,NetCDF,Dates,MeshArrays,Printf
#using MITgcmTools


include("write.jl")
include("read.jl")
include("tile_support.jl")
include("helper_functions.jl")

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData


end # module
