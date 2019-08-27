using Test
using NCTiles, NCDatasets
include("testutils.jl")

println("Testing...")


testvars = maketestdata()

ncfilarray2d,tmp = mktemp()
ncfilarray3d,tmp = mktemp()
ncfilbin2d,tmp = mktemp()
ncfilbin3d,tmp = mktemp()
ncfilnc2d,tmp = mktemp()
ncfilnc3d,tmp = mktemp()

tempfiles = vcat(testvars["fnames2d"], testvars["fnames3d"], 
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

end

# Delete temporary files
for fil in tempfiles
    rm(fil)
end