# NCTiles.jl

[NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) mainly aims to make it easier to write (and read) [NetCDF](https://en.wikipedia.org/wiki/NetCDF) files that represent e.g. the whole Earth surface or a subdomain (`tiles`) based on the [CF conventions](http://cfconventions.org). 

Goals of [NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) include (1) inter-operability with climate model grids via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl), [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl), and [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl); (2) generation of [CF-compliant](http://cfconventions.org) files from suitable array formats. 

```@contents
Pages = [
    "maindocs.md",
    "examples.md",
]
Depth = 2
```

_`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015))._
