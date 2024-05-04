# Examples

One series of use examples is provided [over here](https://nbviewer.jupyter.org/github/JuliaClimate/GlobalOceanNotebooks/blob/master/DataStructures/03_nctiles.ipynb) in the [JuliaClimate Notebooks](https://juliaclimate.github.io/GlobalOceanNotebooks/).

More detailed examples listed below rely on the [nctiles-testcases](https://github.com/gaelforget/nctiles-testcases) and are found in the `NCTiles.jl/examples/` folder.

- [Example1.jl](../generated/Example1/) reads two-dimensional fields on a regular grid ("lat-lon") read from binary files, and then writes them to a netcdf file. This example illustrates the use of either `NCDatasets.jl` or `NetCDF.jl` as the backend.
- [Example2.jl](../generated/Example2/) reads two-dimensional fields from the netcdf file generated in `Example1`, and then re-writes them to a new netcdf file.
- [Example3.jl](../generated/Example3/) reads Global Ocean variables which are partitioned into subdomains and writes each one to a collection of `NetCDF` files ( _nctiles_ ).
- [Example4.jl](../generated/Example4/) reads two three-dimensional variables from the netcdf files generated in `Example3`, combines them into a single data structure, and then re-writes them together into a new netcdf file.
- [Example5.jl](../generated/Example5/) writes a `ClimArray` struct from the `ClimateBase.jl` package to a netcdf file using `NCTiles.jl`, and vice versa.
- [Example6.jl](../generated/Example6/) (_Work In Progress_) demonstrates the specification of a climatological time axis.
