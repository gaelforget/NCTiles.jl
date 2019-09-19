using Test
using NCTiles, NCDatasets, MeshArrays
include("testutils.jl")

println("Testing...")


testvars = maketestdata()

ncfilarray2d,tmp = mktemp()
ncfilarray3d,tmp = mktemp()
ncfilbin2d,tmp = mktemp()
ncfilbin3d,tmp = mktemp()
ncfilnc2d,tmp = mktemp()
ncfilnc3d,tmp = mktemp()
ncfiltile2d,tmp = mktemp()

tempfiles = vcat(testvars["fnames2d"], testvars["fnames3d"],
                    testvars["tile_ex"]["fnamestile2d"],
                    ncfilarray2d, ncfilarray3d,
                    ncfilbin2d, ncfilbin3d,
                    ncfilnc2d, ncfilnc3d)

@testset "NCtiles Tests" begin

    # Test in memory write
    @testset "In Memory Data" begin
        fld2d = NCvar("data2d",testvars["units"],[testvars["dims"][1:2]; testvars["dims"][4]],testvars["data2d"],Dict("long_name" =>testvars["longname2d"]),NCDatasets)
        fld3d = NCvar("data3d",testvars["units"],testvars["dims"],testvars["data3d"],Dict("long_name" =>testvars["longname3d"]),NCDatasets)

        writetestfile(ncfilarray2d,fld2d,NCDatasets)
        writetestfile(ncfilarray3d,fld3d,NCDatasets)

        @test testfile(ncfilarray2d,fld2d)
        @test testfile(ncfilarray3d,fld3d)
    end

    # Test BinData write
    @testset "Binary File Data" begin
        bdatafld2d = BinData(testvars["fnames2d"],Float32,(testvars["dims"][1].dims,testvars["dims"][2].dims))
        bdatafld3d = BinData(testvars["fnames3d"],Float32,(testvars["dims"][1].dims,testvars["dims"][2].dims,testvars["dims"][3].dims))

        bfld2d = NCvar("data2d",testvars["units"],[testvars["dims"][1:2]; testvars["dims"][4]],bdatafld2d,Dict("long_name" =>testvars["longname2d"]),NCDatasets)
        bfld3d = NCvar("data3d",testvars["units"],testvars["dims"],bdatafld3d,Dict("long_name" =>testvars["longname3d"]),NCDatasets)

        writetestfile(ncfilbin2d,bfld2d,NCDatasets)
        writetestfile(ncfilbin3d,bfld3d,NCDatasets)

        @test testfile(ncfilbin2d,bfld2d)
        @test testfile(ncfilbin3d,bfld3d)
    end

    #Test NCData write
    @testset "NetCDF File Data" begin
        ncdatafld2d = NCData(ncfilarray2d,"data2d",NCDatasets,Float32)
        ncdatafld3d = NCData(ncfilarray3d,"data3d",NCDatasets,Float32)

        ncfld2d = NCvar("data2d",testvars["units"],[testvars["dims"][1:2]; testvars["dims"][4]],ncdatafld2d,Dict("long_name" =>testvars["longname2d"]),NCDatasets)
        ncfld3d = NCvar("data3d",testvars["units"],testvars["dims"],ncdatafld3d,Dict("long_name" =>testvars["longname3d"]),NCDatasets)

        writetestfile(ncfilnc2d,ncfld2d,NCDatasets)
        writetestfile(ncfilnc3d,ncfld3d,NCDatasets)

        @test testfile(ncfilnc2d,ncfld2d)
        @test testfile(ncfilnc3d,ncfld3d)
    end

    @testset "Tile Data" begin
        tilesize = testvars["tile_ex"]["tilesize"]
        dims = testvars["tile_ex"]["dims"]
        grid = testvars["tile_ex"]["grid"]
        gridvars = GridLoad(grid)
        land = gridvars["hFacC"]
        for f in land.fIndex
            for d in 1:size(land,2)
                land[f,d][land[f,d].==0] .= NaN
                land[f,d][land[f,d].>0] .= 1
            end
        end
        tiledatafld2d = BinData(testvars["tile_ex"]["fnamestile2d"],Float32,Tuple(grid.ioSize))
        tilfld2d = TileData(tiledatafld2d,testvars["tile_ex"]["tilesize"],grid)
        tillat = TileData(gridvars["YC"],tilfld2d.tileinfo,tilfld2d.tilesize,tilfld2d.precision,tilfld2d.numtiles)
        tillon = TileData(gridvars["XC"],tilfld2d.tileinfo,tilfld2d.tilesize,tilfld2d.precision,tilfld2d.numtiles)
        tilarea = TileData(gridvars["RAC"],tilfld2d.tileinfo,tilfld2d.tilesize,tilfld2d.precision,tilfld2d.numtiles)
        tilland = TileData(land,tilfld2d.tileinfo,tilfld2d.tilesize,tilfld2d.precision,tilfld2d.numtiles)
        thic = gridvars["RC"][:,1]
        flds = Dict(["data2d" => NCvar("data2d",testvars["units"],[dims[1:2]; dims[4]],tilfld2d,Dict(),NCDatasets),
                    "lon" => NCvar("lon","degrees_east",dims[1:2],tillon,Dict("long_name" => "longitude"),NCDatasets),
                    "lat" => NCvar("lat","degrees_north",dims[1:2],tillat,Dict("long_name" => "latitude"),NCDatasets),
                    "area" => NCvar("area","m^2",dims[1:2],tilarea,Dict(["long_name" => "grid cell area", "standard_name" => "cell_area"]),NCDatasets),
                    "land" => NCvar("land","1",dims[1:3],tilland,Dict(["long_name" => "land mask", "standard_name" => "land_binary_mask"]),NCDatasets),
                    "thic" => NCvar("thic","m",dims[3],thic,Dict("standard_name" => "cell_thickness"),NCDatasets)
        ])
        writetestfile(ncfiltile2d,flds,NCDatasets)
        savenames = joinpath.(ncfiltile2d*".".*lpad.(string.(1:tilfld2d.numtiles),4,"0").*".nc")

        @test all([testfile(fname,flds[fld]) for fname in savenames for fld in keys(flds)])
    end

end

# Delete temporary files
for fil in tempfiles
    rm(fil)
end
