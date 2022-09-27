using Documenter, JADE

makedocs(
    sitename = "JADE: a Julia DOASA Environment",
    modules = [JADE],
    clean = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["JADE" => "index.md", "API Reference" => "api.md"],
)

deploydocs(repo = "github.com/EPOC-NZ/JADE.git", devurl = "docs")
