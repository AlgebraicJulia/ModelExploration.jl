module Curry

export uncurry, ACSetCat

using DataStructures
using Catlab.CategoricalAlgebra
using Catlab.Present
import Catlab.CategoricalAlgebra: FinDomFunctor, do_ob_map, do_hom_map
import Catlab.CategoricalAlgebra.FinCats: FinTransformationMap
using AutoHashEquals
using ..CatLimits
import ..CatLimits: product, ob_map, hom_map
import Catlab.Theories: dom, codom, curry

const FinSetCat = TypeCat{SetOb,FinDomFunction{Int}}
const ACSetDomCat = FinCats.FinCatPresentation{
  Symbol, Union{FreeSchema.Ob,FreeSchema.AttrType},
  Union{FreeSchema.Hom,FreeSchema.Attr,FreeSchema.AttrType}}

@auto_hash_equals struct ACSetFunctorEq{ACS<:ACSet} <: Functor{ACSetDomCat,FinSetCat}
  acset::ACS
  eqs::Vector{Pair}
end

FinDomFunctor(X::ACSet; eqs=Pair[]) = ACSetFunctorEq(X, eqs)

function dom(F::ACSetFunctorEq)
    p = Presentation(F.acset)
    for (l,r) in F.eqs
      add_equation!(p, l, r)
    end
    FinCat(p)
end
codom(F::ACSetFunctorEq) = TypeCat{SetOb,FinDomFunction{Int}}()

Categories.do_ob_map(F::ACSetFunctorEq, x) = SetOb(F.acset, Symbol(x))
Categories.do_hom_map(F::ACSetFunctorEq, f) = SetFunction(F.acset, Symbol(f))

function (C::Type{ACS})(F::FinTransformationMap) where ACS <: ACSet
    Cd, CCd = C(dom(F)), C(codom(F))
    return CSetTransformation(Cd, CCd; components(F)...)
end
  
  
FinTransformationMap(f::ACSetTransformation; eqs=Pair[]) =
  FinTransformationMap(components(f),
                      FinDomFunctor(dom(f); eqs=eqs),
                      FinDomFunctor(codom(f); eqs=eqs)
)

                      
is_hom_equal(f::ACSetTransformation, g::ACSetTransformation) =
  force(f) == force(g)

# Tensor-hom adjunction (currying of diagrams in C-Set)
#######################################################
const ACSetCat{S} = TypeCat{S, ACSetTransformation}

"""    curry(d::FinFunctor{D, ACSetCat{S}}) where {D,S}
Currying on objects of a functor category
"""
function curry(d::FinFunctor{D, ACSetCat{S}}) where {D,S}
  shapelim = product([dom(d), FinCat(Presentation(S))])
  shape_ind, part_ind = legs(shapelim)
  asl = apex(shapelim)
  omap = Dict(map(ob_generators(asl)) do o
    x = ob_map(shape_ind, o)
    y = ob_map(part_ind, o)
    o => FinSet(ob_map(d, x), Symbol(y))
  end)

  hmap = Dict(map(hom_generators(asl)) do o
    x = hom_map(shape_ind, o)
    y = hom_map(part_ind, o)
    if x isa FreeSchema.Hom{:id}
      o => FinFunction(ob_map(d, only(x.args)), Symbol(y))
    elseif y isa FreeSchema.Hom{:id}
      o => hom_map(d, x)[Symbol(only(y.args))]
    else
      error("x $x y $y")
    end
  end)

  FinDomFunctor(omap,hmap,asl,FinSetCat())
end

"""
Uses an example FinDomFunctor (in the original uncurried format).
"""
function uncurry(d::FinDomFunctor{D1, FinSetCat},
                 old_d::FinDomFunctor{D2, ACSetCat{S}}) where {D1,D2,S}
  # Recover schema for d as a product, not just the apex
  shapelim = product([dom(old_d), FinCat(Presentation(S))])
  asl = apex(shapelim)
  shape_ind, part_ind = legs(shapelim)

  cset_type = typeof(first(old_d.ob_map)[2])
  omap = Dict(map(ob_generators(dom(old_d))) do o
    x = Base.invokelatest(cset_type)
    for o_ in ob_generators(asl)
      if ob_map(shape_ind, o_) == o
        add_parts!(x, Symbol(ob_map(part_ind, o_)), length(ob_map(d, o_)))
      end
    end
    for h in hom_generators(asl)
      h_ = hom_map(shape_ind, h)
      if h_ == id(o)
        set_subpart!(x, Symbol(hom_map(part_ind, h)), collect(hom_map(d, h)))
      end
    end
    o => x
  end)
  hmap = Dict(map(hom_generators(dom(old_d))) do h
    comps = Dict()
    for h_ in hom_generators(asl)
      if hom_map(shape_ind, h_) == h
        comps[Symbol(only(hom_map(part_ind, h_).args))] = hom_map(d, h_)
      end
    end
    dom_, codom_ = [omap[get(h)] for get in [dom, codom]]
    h => ACSetTransformation(dom_,codom_; comps...)
  end)
  FinDomFunctor(omap,hmap,dom(old_d),ACSetCat{S}())
end

"""    curry(d::FinFunctor{D, ACSetCat{S}}) where {D,S}
Currying on morphisms of a functor category with an ACSetCat as codom
"""
function curry(ϕ::FinTransformationMap{D, ACSetCat{S}}) where {D,S}
  cur_d, cur_cd = curry.([dom(ϕ), codom(ϕ)])
  shapelim = product([dom(dom(ϕ)), FinCat(Presentation(S))])
  shape_ind, part_ind = legs(shapelim)
  comps = Dict(map(ob_generators(apex(shapelim))) do o
    oshape, opart = Symbol(shape_ind(o)), Symbol(part_ind(o))
    Symbol(o) => components(ϕ)[oshape][opart]
  end)
  FinTransformationMap(comps,cur_d,cur_cd)
end

"""    uncurry(d::FinTransformationMap, old_d::FinTransformationMap{D, ACSetCat{S}}) where {D, S}
Inverse to currying on morphisms of a functor category with an ACSetCat as codom
"""
function uncurry(d::FinTransformationMap,
                 old_d::FinTransformationMap{D, ACSetCat{S}}) where {D, S}
  # Recover schema for d as a product, not just the apex
  shapelim = product([dom(dom(old_d)), FinCat(Presentation(S))])
  shape_ind, part_ind = legs(shapelim)

  αcomps = Dict(o => DefaultDict{Symbol,Vector{Int}}(()->Int[])
                for o in keys(components(old_d)))

  for o in (ob_generators(apex(shapelim)))
    dic = αcomps[Symbol(ob_map(shape_ind, o))]
    dic[Symbol(ob_map(part_ind, o))] = collect(components(d)[Symbol(o)])
  end

  uc_d, uc_cd = [uncurry(get(d), get(old_d)) for get in [dom, codom]]

  α = Dict(map(collect(αcomps)) do (o, comps)
    o => ACSetTransformation(ob_map(uc_d,   o),
                             ob_map(uc_cd, o); comps...)
  end)


  FinTransformationMap(α, uc_d, uc_cd)
end

end