# # Example 3 : Tiled Domain
#
# This example first reads in Global Ocean variables, generated by the [MITgcm](https://mitgcm.readthedocs.io/en/latest/index.html), which are split into a collection of subdomain arrays ( _tiles_ ) using `MeshArrays.jl`. It then writes them to a collection of `NetCDF` files ( _nctiles_ ) using `NCTiles.jl`.

# ### Packages & Helper Functions
#

using NCTiles; 

p=dirname(pathof(NCTiles))
fil = joinpath(p, "../examples/helper_functions.jl")
include(fil);

# ### File Paths & I/O Back-End
#

inputs=NCTiles.NCTILES_TESTCASES
NCTiles.ensure_testcases_installed()

outputs = joinpath(tempdir(),"NCTILES_TESTCASES_OUTPUT/ex3/")
if ~ispath(outputs); mkpath(outputs); end

nc=NCTiles.NCDatasets # I/O Back-End

# ## Process Global Ocean Variables
#
# Here we process a two-dimensional field (`ETAN`), two three-dimensional tracer fields (`THETA`, `SALT`), and the three components of a vector field (`UVELMASS`, `VVELMASS` and `WVELMASS`). In each case, `flds` denotes the lazy representation of the processing chain (incl. all needed metadata) which the `write` function instantiates (i.e. outputs to files).
#
# _Note: on a `C-grid` these components are staggered in space._

writedir=outputs
readme = readlines(inputs)

# ### 2D example
(flds,savename,readme)=prep_nctiles_native(inputs,"state_2d_set1","ETAN",Float32)
write(flds,savename,README=readme);

# ### 3D examples
(flds,savename,readme)=prep_nctiles_native(inputs,"state_3d_set1","THETA",Float32);
write(flds,savename,README=readme);

(flds,savename,readme)=prep_nctiles_native(inputs,"state_3d_set1","SALT",Float32);
write(flds,savename,README=readme);

# ### 3D vector field examples
(flds,savename,readme)=prep_nctiles_native(inputs,"trsp_3d_set1","UVELMASS",Float32);
write(flds,savename,README=readme);

(flds,savename,readme)=prep_nctiles_native(inputs,"trsp_3d_set1","VVELMASS",Float32);
write(flds,savename,README=readme);

(flds,savename,readme)=prep_nctiles_native(inputs,"trsp_3d_set1","WVELMASS",Float32);
write(flds,savename,README=readme);


