module NCTiles

using NCDatasets,NetCDF,MeshArrays
using Dates,Printf,Pkg, Pkg.Artifacts
#using MITgcmTools

version()=Pkg.TOML.parsefile(joinpath(dirname(pathof(NCTiles)), "..", "Project.toml"))["version"]

p=dirname(pathof(NCTiles))
artifact_toml = joinpath(p, "../Artifacts.toml")
NCTILES_TESTCASES_hash = artifact_hash("NCTILES_TESTCASES", artifact_toml)
NCTILES_TESTCASES = joinpath(artifact_path(NCTILES_TESTCASES_hash)*"/","nctiles-testcases-0.1/")

"""
    pkg_version(name::String)

Get the version number (in current environment) for a package name 
"""
function pkg_version(name::String)
    tmp1=Pkg.dependencies()
    tmp2=Vector{Union{Missing, VersionNumber}}(missing, 1)
    for i in collect(keys(tmp1))
        if tmp1[i].name==name
            tmp2[1]=tmp1[i].version
        end
    end
    return tmp2[1]
end

include("write.jl")
include("read.jl")
include("tile_support.jl")
include("helper_functions.jl")

export NCvar, BinData, NCData, TileData, readbin, readncfile
export createfile, addDim, addVar, addData, addDimData
export write


end # module
