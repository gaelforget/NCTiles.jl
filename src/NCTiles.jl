module NCTiles

using NetCDF, MeshArrays, Printf, NCDatasets

include("utilities.jl")
export read_nctiles

import Base: getindex

export NCvar, Bindata, NCData, TileData, readbin
export createfile, addDim, addVar, addData, addDimData,
        readncfile, parsemeta, readAvailDiagnosticsLog,
        findtiles
        
include("writenctiles.jl")

end # module
