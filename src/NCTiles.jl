module NCTiles

using NetCDF, MeshArrays, Printf, NCDatasets

include("utilities.jl")
export read_nctiles

import Base: getindex

export NCvar, Bindata, NCData, readbin
export createfile, addDim, addVar, addData, addDimData,
        readncfile, parsemeta, readAvailDiagnosticsLog
        
include("writenctiles.jl")

end # module
