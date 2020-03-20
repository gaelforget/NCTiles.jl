module NCTiles

using NCDatasets,NetCDF,MeshArrays
using Dates,Printf,Pkg
#using MITgcmTools

version()=Pkg.TOML.parsefile(joinpath(dirname(pathof(NCTiles)), "..", "Project.toml"))["version"]

include("write.jl")
include("read.jl")
include("tile_support.jl")
include("helper_functions.jl")

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData
export write


end # module
