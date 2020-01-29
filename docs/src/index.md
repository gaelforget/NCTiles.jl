# `NCTiles.jl` Package Documentation

**NCTiles.jl** reads and writes [NetCDF files](https://en.wikipedia.org/wiki/NetCDF) that represent e.g. subdivisions of Earth's surface (`tiles`). Inter-operability with popular climate model grids and [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and generation of [CF-compliant](http://cfconventions.org) files are key goals of this package. 

## Contents

```@contents
```

_`NCTiles.jl` derives from the earlier `nctiles` implementation in [gcmfaces](https://github.com/MITgcm/gcmfaces) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015))._

## Main Features

`NCTiles.jl` first goes through lazy operations, on data stucture, as it obtains information about variables etc. The second phase calls `createfile` and then the `add*` functions to instantiate files. _Note:_ some of the included functions are interfaces to `MITgcm` output.

Data structures:

- `NCData` contains a string or an array of strings (NetCDF file names) + metadata to read files. 
- `NCvar` contains information needed to write a NetCDF file which can include a list of filenames (see `Bindata`) if the data is not loaded into memory.
- `BinData` is a container for one field.
- `TileData ` contains a `MeshArray` or `BinData` struct in `vals`,
    information about the tile layout in `tileinfo`, and metadata needed for
    reading/writing tile data.

As an example:
    
```
struct TileData{T}
    vals::T
    tileinfo::Dict
    tilesize::Tuple
    precision::Type
    numtiles::Int
end
```    

## Use examples

`DataStructures/06_nctiles.ipynb` in this [GlobalOceanNotebooks repo](https://github.com/gaelforget/GlobalOceanNotebooks/) provides a series of examples that overlap somewhat with the ones found in `Examples/ex*.jl`:

- `ex_1.jl` reads a `binary` file containing one interpolated 2D field on a regular grid. It then writes that array to a `NetCDF`/`NCTiles` file.
- `ex_2.jl` reads data from a `NetCDF` file containing one `tile` of model output. It then writes it to a new `NetCDF`/`NCTiles` file. This uses 3D data on a non-regular grid for one ocean subdivision (`tile`).
- `ex_3.jl` is an example of interpolated model output processing in `CBIOMES` where several variables are included in the same `NetCDF`/`NCTiles` file.
- `ex_4.jl` generates a tiled netcdf output (i.e., a `nctiles` output) for a global 2D field on the non-regular `LLC90` grid (see `MeshArrays.jl`). Since the tile width is set to 90, this creates 13 files.

## Index

```@index
```

## API / Functions

```@autodocs
Modules = [NCTiles]
Order   = [:type,:function]
Private = false
```
