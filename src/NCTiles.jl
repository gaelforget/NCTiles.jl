module NCTiles

using NCDatasets,NetCDF,Dates,MeshArrays,Printf
#using MITgcmTools

include("Write.jl")
include("Read.jl")
include("TileSupport.jl")
include("HelperFunctions.jl")

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData


end # module
