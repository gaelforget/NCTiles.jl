using Documenter
using Pkg
using NCTiles

p=dirname(pathof(NCTiles))
artifact_toml = joinpath(p, "../Artifacts.toml")
Pkg.ensure_artifact_installed("NCTILES_TESTCASES","Artifacts.toml")

function gunzip_testcases_if_needed()
    fils=["state_3d_set1.0000000732.data.gz","state_3d_set1.0000001428.data.gz","state_3d_set1.0000002172.data.gz",
    "trsp_3d_set1.0000000732.data.gz","trsp_3d_set1.0000001428.data.gz","trsp_3d_set1.0000002172.data.gz"]
    pth=joinpath(NCTiles.NCTILES_TESTCASES,"diags/")
    for i in 1:length(fils)
        if isfile(joinpath(pth,fils[i]))
            run(`gunzip $(joinpath(pth,fils[i]))`)
        end
    end
end

makedocs(
    sitename = "NCTiles",
    format = Documenter.HTML(),
   pages=[
        "Home" => "index.md",
        "maindocs.md",
        "examples.md",
        "API.md"],
    modules = [NCTiles]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gaelforget/NCTiles.jl.git",
)
