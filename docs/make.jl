using Documenter
using NCTiles

makedocs(
    sitename = "NCTiles",
    format = Documenter.HTML(),
    modules = [NCTiles]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
