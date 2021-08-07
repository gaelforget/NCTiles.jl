# `NCTiles.jl` Package Documentation

**NCTiles.jl** reads and writes [NetCDF files](https://en.wikipedia.org/wiki/NetCDF) that represent e.g. subdivisions of Earth's surface (`tiles`). Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and generation of [CF-compliant](http://cfconventions.org) files are key goals of this package. 

```@contents
Pages = [
    "maindocs.md",
    "examples.md",
    "API.md",
]
Depth = 2
```

_`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015))._
