module ModelExploration
export Explore, Dim, Free

using Catlab.CategoricalAlgebra




abstract type Dim end

struct Free <: Dim end

struct Explore
  csettype::Type
  slice::StructACSet
  init::StructACSet
  dims::Vector{Dim}
end

end # module
