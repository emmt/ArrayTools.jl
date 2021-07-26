using Documenter

push!(LOAD_PATH, "../src/")
using ArrayTools, ArrayTools.PseudoArrays

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "ArrayTools package for Julia",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = ["index.md", "install.md", "storage.md", "indexing.md",
             "rubberindex.md", "broadcasting.md", "reference.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/ArrayTools.jl.git",
    )
end
