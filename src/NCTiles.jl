module NCTiles

using NCDatasets

greet() = print("Hello World!")

include("writenctiles.jl")

export NCvar, Bindata, readbin
export createfile, addDim, addData, addDimData

end # module
