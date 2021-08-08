using Documenter
using NCTiles

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
