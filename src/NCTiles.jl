module NCTiles

using NetCDF, MeshArrays, Printf, NCDatasets
#using MITgcmTools

include("writenctiles.jl")
include("utilities.jl")

import Base: getindex

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData


end # module
