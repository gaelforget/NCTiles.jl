# NCTiles.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/NCTiles.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/NCTiles.jl/dev)
[![Travis Build Status](https://travis-ci.org/gaelforget/NCTiles.jl.svg?branch=master)](https://travis-ci.org/gaelforget/NCTiles.jl)
[![DOI](https://zenodo.org/badge/179139682.svg)](https://zenodo.org/badge/latestdoi/179139682)

[NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) mainly reads and writes [NetCDF](https://en.wikipedia.org/wiki/NetCDF) files (`nc`) that represent e.g. the whole Earth surface or a subdomain (`tiles`) based on the [CF coventions](http://cfconventions.org). Goals of [NCTiles.jl](https://github.com/gaelforget/NCTiles.jl) include (1) inter-operability with popular climate model grids via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl); (2) making it easy to generate [CF-compliant](http://cfconventions.org) files from suitable Array formats. 

