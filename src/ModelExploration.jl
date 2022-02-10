module ModelExploration

export Product, Literal, Gen, SliceLit, select, unfold, get_infected_states, sum_infected_states, generate_data, make_loss, train, eval_petri_fn

include("Petri.jl")
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra: pullback
using Catlab.Graphs
using Catlab.Theories: attrtype

#############
abstract type Gen end

struct Product <: Gen
  g1::Gen
  g2::Gen
  slice::StructACSet
end

struct Literal <: Gen
  lit::Vector{StructACSet}
end

struct SliceLit <: Gen
  base::StructACSet
  models::Vector{Pair{StructACSet, NamedTuple}}
end

# Not thinking too far ahead on this for now lol
function unfold(g::Literal)::Vector{StructACSet}
  return g.lit
end
getS(::StructACSet{S}) where {S} = S
function unfold(g::SliceLit; loose::Bool=true)::Vector{ACSetTransformation}
  transformations = ACSetTransformation[]
  S = getS(g.base)
  b = map(g.base; Dict([at => x->nothing for at in attrtype(S)])...)
  for (model, constraint) in g.models
    append!(transformations,
            homomorphisms(model, b; loose=loose, initial=constraint))
  end
  return transformations
end

"""Assumes slice generators. TODO: search homomorphisms otherwise"""
function unfold(g::Product)::Vector{StructACSet}
  models = StructACSet[]
  for x in unfold(g.g1)
    for y in unfold(g.g2)
      push!(models, ob(pullback(x, y)))
    end
  end
  return models
end

function select(g::Gen, lossFn::Function)::StructACSet
  
  best = typemax(Float64)
  best_m = nothing
  for (i,m) in enumerate(unfold(g)[1:9])
      println("i $i")
      lss = lossFn(m)
      println("\tloss $lss")
      if lss < best 
        best_m = m 
      end 
  end
  return best_m 
end

end # module