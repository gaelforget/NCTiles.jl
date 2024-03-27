using Documenter, Literate, NCTiles, Pkg

pth=@__DIR__
lst=("Example1.jl","Example2.jl","Example3.jl","Example4.jl","Example5.jl","Example6.jl")
lstExecute=("Example1.jl","Example2.jl","Example3.jl","Example4.jl","Example5.jl","Example6.jl")
for i in lst
    EXAMPLE = joinpath(pth, "..", "examples", i)
    OUTPUT = joinpath(pth, "src","generated")
    Pkg.activate(joinpath(pth,"..","docs"))
    Literate.markdown(EXAMPLE, OUTPUT, documenter = true)
    cd(pth)
    Pkg.activate(joinpath(pth,"..","docs"))
    tmp=xor(occursin.(i,lstExecute)...)
    Literate.notebook(EXAMPLE, OUTPUT, execute = tmp)
    cd(pth)
    #Literate.notebook(EXAMPLE, OUTPUT, flavor = :pluto)
end

makedocs(
    sitename = "NCTiles",
    format = Documenter.HTML(),
   pages=[
        "Home" => "index.md",
        "Examples" => Any[
            "Guide " => "examples.md",
            "Listing" => [map(s -> "generated/$(s[1:end-2])md",lst)...],
            ],
        "maindocs.md",
        "API.md"],
    warnonly = [:cross_references,:missing_docs],
    modules = [NCTiles]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gaelforget/NCTiles.jl.git",
)
