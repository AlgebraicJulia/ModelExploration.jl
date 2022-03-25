module ModelExploration

export Product, Literal, ModelSpace, PullbackSpace, ModelHom, SliceHom, select,
       unfold, to_model_hom, SliceCat

using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra: pullback
using Catlab.Graphs, Catlab.Present
using Catlab.Theories: attrtype

#############
const SliceCat{Ob,Hom} = TypeCat{Slice{Ob,Hom}, SliceMorphism{Ob,Hom}}

abstract type ModelSpace{T,Base} end
abstract type FunctorGen end
abstract type NatGen end


"""
Base case: a explicit diagram considered as a model space

E.g.

1. Shape: •, and diagram maps this object to the petri net 'infectious_type'
2. Shape: • •, and diagram maps these to SIR and SIS.
3. Shape: • •, and diagram maps these to SIR->I.T. and SIS->I.T morphisms
"""
struct Literal{T,Base} <: ModelSpace{T,Base}
  lit::Diagram{T,Base}
end

"""
A recipe for constructing a modelspace (diagram) homomorphism.

E.g.
dom = Literal example #2
codom = Literal example #1
functor F = uniquely determined
nat= homomorphisms from SIR and SIS to Infectious_type (explained below)

Dom => F;Codom :: Functors from shape(dom) -> Petri
data: map from ob(shape(dom)) -> Hom(Petri). αᵢ:: F(i) -> G(i)
only morphisms are identities, so naturality is automatically satisfied.
Therefore, αᵢ merely requires a homomorphism from SIR and SIS to Infectious_type

Therefore, this is equivalent to Literal example #3.
"""
struct ModelHom{T,Base} <: ModelSpace{T,Base}
  dom::ModelSpace{T,Base}
  codom::ModelSpace{T,Base}
  fun::FunctorGen
  nat::NatGen
end

"""
Model space generated by pullback of two diagrams
(sliced over the same base diagram)
"""
struct PullbackSpace{T,Base} <: ModelSpace{T,Base}
  g1::ModelSpace{T,Base}
  g2::ModelSpace{T,Base}
end

struct LitFG <: FunctorGen
  lit::Functor
end
struct LitNG <: NatGen
  lit::Dict
end
struct UniqueFG <: FunctorGen end
struct UniqueNG <: NatGen end
struct AllFG <: FunctorGen end
struct AllNG <: NatGen end
# idea: constrain the search for functors/nat transformations declarative
# with positive/negative application conditions.

"""
Assume that the model space given as an argument is diagram in a slice category.

Generates a diagram hom between diagrams in the underlying category of the slice
category.

For example, Diagram{T,Petri/X} is converted into DiagramHom{T,Petri}, where
the codomain has shape • and sends that to the Petri net X.
"""
struct SliceHom{T,BaseOb,BaseHom} <: ModelSpace{T,SliceCat{BaseOb,BaseHom}}
  slice_diagram::ModelSpace{T,SliceCat{BaseOb,BaseHom}}
end

@present ThOne(FreeSchema) begin
  X::Ob
end
One = FinCat(ThOne)

unfold(g::Literal)::Diagram = g.lit
unfold(g::SliceHom)::Diagram = unfold(g.slice_diagram)

"""
Generate the data of a diagram homomorphism from a slice hom
"""
function to_model_hom(g::SliceHom{T,BaseOb,BaseHom}) where {T,BaseOb,BaseHom}
  Base = TypeCat{BaseOb,BaseHom}
  slicehom = diagram(unfold(g)) # underlying functor of diagram
  d_ob = Dict(map(collect(ob_map(slicehom))) do (k, v)
    k => dom(v.slice)
  end)
  d_hom = Dict(map(collect(hom_map(slicehom))) do (k, v)
    k => v.f
  end)
  X = codom(last(first(collect(ob_map(slicehom)))).slice)
  d = Diagram{T}(FinDomFunctor(d_ob,d_hom,dom(slicehom), Base()))
  cod =Diagram{T}(FinDomFunctor(Dict(:X=>X),nothing,One, Base()))
  shape_m_ob = Dict([k=>:X for k in keys(d_ob)])
  shape_m_hom = Dict([h=>id(One[:X]) for h in hom_generators(dom(slicehom))])
  shape_m = FinDomFunctor(shape_m_ob, shape_m_hom, dom(slicehom), One)
  diagram_m = Dict(map(collect(ob_map(slicehom))) do (k, v)
    k => v.slice
  end)
  ModelHom{T,Base}(Literal{T,Base}(d), Literal{T,Base}(cod),
                   LitFG(shape_m), LitNG(diagram_m),)
end


unfold(l::LitFG) = l.lit
unfold(l::LitNG) = l.lit

unfold(g::ModelHom{T,Base}) where {T,Base} = DiagramHom{T,Base}(
    unfold(g.fun), unfold(g.nat), unfold(d.dom), unfold(g.codom))

#=function unfold(g::SliceLit; loose::Bool=true)::Vector{ACSetTransformation}
  transformations = ACSetTransformation[]
  S = getS(g.base)
  b = map(g.base; Dict([at => x->nothing for at in attrtype(S)])...)
  for (model, constraint) in g.models
    append!(transformations,
            homomorphisms(model, b; loose=loose, initial=constraint))
  end
  return transformations
end=#

"""Assumes slice Generators. TODO: search homomorphisms otherwise"""
function unfold(g::Pullback)::Vector{StructACSet}
  models = StructACSet[]
  for x in unfold(g.g1)
    for y in unfold(g.g2)
      push!(models, ob(pullback(x, y)))
    end
  end
  return models
end

function select(g::ModelSpace, lossFn::Function)::StructACSet

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