# NCTiles.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/NCTiles.jl/dev)
[![DOI](https://zenodo.org/badge/179139682.svg)](https://zenodo.org/badge/latestdoi/179139682)
[![CI](https://github.com/gaelforget/NCTiles.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/gaelforget/NCTiles.jl/actions/workflows/ci.yml)

[NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) mainly aims to make it easier to write (and read) [NetCDF](https://en.wikipedia.org/wiki/NetCDF) files that represent e.g. the whole Earth surface or a subdomain (`tiles`) based on the [CF conventions](http://cfconventions.org). 

Goals of [NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) include (1) inter-operability with climate model grids via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl), [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl), and [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl); (2) generation of [CF-compliant](http://cfconventions.org) files from suitable array formats. 

