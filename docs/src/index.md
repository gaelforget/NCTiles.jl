# NCTiles.jl Documentation

**NCTiles.jl** reads and writes tiled `Netcdf` files that represent e.g. subdivisions of Earth's surface. Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/gaelforget/MeshArrays.jl) is an important prospect. `NCTiles.jl` derives from the `nctiles` design previously implemented in `Matlab` as part of [gcmfaces](https://github.com/gaelforget/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015>)). 

Provided examples (in `Examples/`) demonstrate package functionalities:
- `ex_1.jl` reads a `binary` file containing one interpolated 2D field on a regular grid. It then writes that array to a `NetCDF`/`NCTiles` file.
- `ex_2.jl` reads data from a `NetCDF` file containing one `tile` of model output. It then writes it to a new `NetCDF`/`NCTiles` file. This uses 3D data on a non-regular grid for one ocean subdivision (`tile`).
- `ex_3.jl` is an example of interpolated model output processing in `CBIOMES` where several variables are included in the same `NetCDF`/`NCTiles` file.
- `ex_4.jl` generates a tiled netcdf output (i.e., a `nctiles` output) for a global 2D field on the non-regular `LLC90` grid (see `MeshArrays.jl`). Since the tile width is set to 90, this creates 13 files.

## Contents

```@contents
```

## Main Features

To be continued...

## Index

```@index
```

## API / Functions

```@autodocs
Modules = [NCTiles]
Order   = [:function, :type]
Private = false
```
