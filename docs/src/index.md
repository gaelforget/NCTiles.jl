# `NCTiles.jl` Package Documentation

[NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) mainly reads and writes [NetCDF](https://en.wikipedia.org/wiki/NetCDF) files (`nc`) that represent e.g. the whole Earth surface or a subdomain (`tiles`) based on the [CF coventions](http://cfconventions.org). 

Goals of [NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) include (1) inter-operability with popular climate model grids via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl); (2) making it easy to generate [CF-compliant](http://cfconventions.org) files from suitable Array formats. 

```@contents
Pages = [
    "maindocs.md",
    "examples.md",
    "API.md",
]
Depth = 2
```

_`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015))._
