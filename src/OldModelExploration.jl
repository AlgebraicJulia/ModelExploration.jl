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

"""

"""
struct Finite <: Dim
  vals::Vector{StructACSetTransformation}
  ordering::Vector{Vector{Int}}
  function Finite(vs::Vector{StructACSetTransformation})
    new(vs, [Int[] for _ in vs])
  end
end

"""
A sequence of models indexed by ℕ
"""
struct Linear <: Dim
  fun::Function # ℕ → CSetTransformation
end

"""
All possible models up to ACSet isomorphism
"""
struct Free <: Dim end # use ModelEnumeration.jl

struct Rewrite <: Dim
  rules::Vector{Pair{ACSetTransformation}}
end

"""
Reductions of model space from all possible ACSets of a schema
"""
abstract type Constraint end


struct Chase <: Constraint
  rules::Vector{EqConstraint}
end

"""
Last resort for a constraint, because model is hidden implicit in code
"""
struct Filter <: Constraint
  f::Function # StructACsetTransformation -> Bool
end

"""
We can only take pullbacks between slices, but we can default to a trivial
slice if none is explicitly provided.
"""
struct Explore
  csettype::Type
  slice::StructACSet
  dims::Vector{Dim}
  constr::Vector{Constraint}
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
