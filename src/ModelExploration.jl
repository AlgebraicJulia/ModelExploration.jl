module ModelExploration
export Explore, Dim, Free

using Catlab.CategoricalAlgebra



"""
A dimension implicitly specifies a metric space of structures, (or, think of it
as a preorder if we don't care about distances) ordered by "complexity" in some
context-dependent sense. We generally want to explore forward in the dimension
only to the extent necessary to model a phenomenon.

The data of a Dim allows us to lazily explore this preorder (the set of
structures implied may be infinite, such as an "N-city model" in epidemiology).
"""
abstract type Dim end

struct Finite <: Dim
  vals::Vector{StructACSet}
end

struct Linear <: Dim
  fun::Function # ℕ → CSet
end

struct Free <: Dim end # use ModelEnumeration.jl

struct Rewrite <: Dim
  rules::Vector{Pair{ACSetTransformation}}
end

struct Explore
  csettype::Type
  slice::StructACSet
  init::StructACSet
  dims::Vector{Dim}
end

struct ModelSpace
  dims::Vector{Dim}
end

struct ModelPoint
  coords::Vector{UInt128} # model hashes in each dimension
end

struct ExploreState
  seen::Set{UInt128}
  space::ModelSpace
end

# Magnitude and direction of greatest improvement?
struct Loss
end

eval(mp::ModelPoint)::Loss = Loss()

stopping_criteria(es::ExploreState)::Bool = false

next_mp(es::ExploreState)::Vector{ModelPoint} = []



end # module
