# Use Examples

`DataStructures/03_nctiles.ipynb` in this [GlobalOceanNotebooks repo](https://github.com/gaelforget/GlobalOceanNotebooks/) provides a series of examples. 

Additional examples, based on [nctiles-testcases](https://github.com/gaelforget/nctiles-testcases), are found in the `examples/` folder:

- `Example1.jl` reads two-dimensional fields on a regular grid ("lat-lon") read from binary files, and then writes them to a netcdf file. This example illustrates the use of either `NCDatasets.jl` or `NetCDF.jl` as the backend.
- `Example2.jl` reads two-dimensional fields from the netcdf file generated in `Example1`, and then re-writes them to a new netcdf file.
- `Example3.jl` reads Global Ocean variables which are partitioned into subdomains and writes each one to a collection of `NetCDF` files ( _nctiles_ ).
- `Example4.jl` reads two three-dimensional variables from the netcdf files generated in `Example3`, combines them into a single data structure, and then re-writes them together into a new netcdf file.
- `Example5.jl` writes a `ClimArray` struct from the `ClimateBase.jl` package to a netcdf file using `NCTiles.jl`, and vice versa.
- `Example6.jl` (_Work In Progress_) demonstrates the specification of a climatological time axis.
