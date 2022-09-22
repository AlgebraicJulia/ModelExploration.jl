using Documenter
using ModelExploration

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/kris-browns/ModelExploration.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/kris-browns/ModelExploration.jl/blob/main/docs"
end

makedocs(
    sitename = "ModelExploration",
    format = Documenter.HTML(),
    modules = [ModelExploration]
)


@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/kris-browns/ModelExploration.jl.git",
  branch = "gh-pages",
  devbranch = "main"
)