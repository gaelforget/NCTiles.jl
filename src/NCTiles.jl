module NCTiles

using NCDatasets,NetCDF,MeshArrays
using Dates,Printf,Pkg, Pkg.Artifacts
#using MITgcmTools

version()=Pkg.TOML.parsefile(joinpath(dirname(pathof(NCTiles)), "..", "Project.toml"))["version"]

p=dirname(pathof(NCTiles))
artifact_toml = joinpath(p, "../Artifacts.toml")
NCTILES_TESTCASES_hash = artifact_hash("NCTILES_TESTCASES", artifact_toml)
NCTILES_TESTCASES = joinpath(artifact_path(NCTILES_TESTCASES_hash)*"/","nctiles-testcases-0.1/")

include("write.jl")
include("read.jl")
include("tile_support.jl")
include("helper_functions.jl")

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData
export write


end # module
