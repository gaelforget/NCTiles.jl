# NCTiles.jl

[![Travis Build Status](https://travis-ci.org/gaelforget/NCTiles.jl.svg?branch=master)](https://travis-ci.org/gaelforget/NCTiles.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/NCTiles.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/NCTiles.jl/dev)

**NCTiles.jl** reads and writes tiled `Netcdf` files that represent e.g. subdivisions of Earth's surface. Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/gaelforget/MeshArrays.jl) is an important prospect. This package derives from the `nctiles` design previously implemented in [gcmfaces](https://github.com/gaelforget/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).

Two examples (`ex_1.jl` and `ex_2.jl`) are currently included. Each demonstrates a potential workflow and tests some of the package functionalities. `ex_1.jl` reads a `binary` file containing interpolated 2D data on a regular, global grid and writes to a `NetCDF` file. `ex_2.jl` reads data from a `NetCDF` file and writes it to a new `NetCDF` file in the case of 3D data on a non-regular grid over an ocean subdivision (`tile`).
