# NCTiles.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/NCTiles.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/NCTiles.jl/dev)
[![Travis Build Status](https://travis-ci.org/gaelforget/NCTiles.jl.svg?branch=master)](https://travis-ci.org/gaelforget/NCTiles.jl)
[![DOI](https://zenodo.org/badge/179139682.svg)](https://zenodo.org/badge/latestdoi/179139682)

**NCTiles.jl** reads and writes [NetCDF files](https://en.wikipedia.org/wiki/NetCDF) that represent e.g. subdivisions of Earth's surface (`tiles`). Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and generation of [CF-compliant](http://cfconventions.org) files are key goals of this package. 

## `Notebooks/`

- See `DataStructures/06_nctiles.ipynb` in this [GlobalOceanNotebooks repo](https://github.com/gaelforget/GlobalOceanNotebooks/)

## `Examples/`

- Several examples illustrate the package functionalities in detail (see `examples/ex*.jl`)

## Background

`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).


