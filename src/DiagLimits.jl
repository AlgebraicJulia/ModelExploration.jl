module DiagLimits

export leftkan, lk_universal

using DataStructures
using Catlab.Present
using Catlab.Theories
using Catlab.CSetDataStructures: AnonACSetType
using Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra: is_natural, product, coproduct, equalizer, coequalizer, universal, pullback, pushout
using Catlab.CategoricalAlgebra.FinCats: FinCatPresentation, FinTransformationMap
using ..Chase
using ..Curry
using ..CatLimits
using ..Misc
using AutoHashEquals

is_natural(D::DiagramHom; verbose::Bool=false) =
    is_functorial(shape_map(D)) && (verbose ? is_natural_verbose(diagram_map(D)) : is_natural(diagram_map(D)))

# Left Kan extensions of Diagrams in Set or C-set along FinFunctors
################################################################

"""    leftkan(F::FinFunctor, I::StructACSet, cd::Symbol; n=15, verbose=false)
    B
F ↗ η⇑ ↘ LanF(I)
 A  ⟶  C
    I
Computes the left kan extension of I (a functor to Set) by F (a shape functor).
"""
function leftkan(F::FinFunctor, I::StructACSet, cd::Symbol; n=15, verbose=false)
  # Assemble chase input using the collage of F as a schema
  col = apex(collage(F))
  name, cname = Symbol("collage_F_$cd"), Symbol("cset_collage_F_$cd")
  col_eds = pres_to_eds(presentation(col), name)
  col_I = Base.invokelatest(crel_type(presentation(col), name))
  col_I_cset = Base.invokelatest(AnonACSetType(presentation(col)))
  copy_parts!(col_I, to_c_rel(I))

  # Run the chase
  chase_res_, check = chase_crel(col_I, col_eds, n; I_is_crel=true,
                       Σ_is_crel=true, cset_example=col_I_cset, verbose=verbose)
  check || error("Chase failed to terminate")
  chase_res = codom(chase_res_)

  # Project the codom portion of the collage and grab the α components
  res = Base.invokelatest(AnonACSetType(presentation(codom(F))))
  copy_parts!(res, chase_res)
  α = Dict(o=>FinFunction(chase_res, Symbol("α_$o"))
           for o in ob_generators(dom(F)))

  # Return result as a DiagramHom{id}
  ddom = Diagram(FinDomFunctor(I; eqs=equations(dom(F))))
  dcodom = Diagram(FinDomFunctor(res; eqs=equations(codom(F))))
  DiagramHom{id}(F, α, ddom, dcodom)
end

"""
Reduce left kan of a functor to C-Set to a computation for a functor into FinSet
Convert result back to a functor into C-Set.
"""
function leftkan(F::FinFunctor{D,CD}, I::FinDomFunctor{D, ACSetCat{S}},
                 cd::Symbol; n=15, verbose=false) where {D, CD, S}
    Ic = curry(I)
    ctype = AnonACSetType(presentation(dom(Ic)))
    Ic_ = ctype(Ic)

    # Map from D x S -> CD x S based on map F: D->CD
    domprod = product(FinCatPresentation[dom(F), FinCat(Presentation(S))])
    domleg1, domleg2 = legs(domprod)
    codomprod = product(FinCatPresentation[codom(F), FinCat(Presentation(S))])
    Fc = universal(codomprod, Multispan([compose(domleg1, F), domleg2]))

    # Compute result for Functors to Set
    curr_res = leftkan(Fc, Ic_, Symbol("random2_$cd"); n=n, verbose=verbose)

    # Create a meaningless FinFunctor from CD to the C-Set category for uncurry
    cset_rep = last(first(ob_map(I))); cset_hom = id(cset_rep)
    fakeob = Dict([o => cset_rep for o in ob_generators(codom(F))])
    fakehom = Dict([o => cset_hom for o in hom_generators(codom(F))])
    cd_finfun = FinDomFunctor(fakeob, fakehom, codom(F),codom(I))

    # Uncurry result
    ucodom = uncurry(diagram(codom(curr_res)), cd_finfun)
    αcomps = Dict(o => DefaultDict{Symbol,Vector{Int}}(()->Int[])
                  for o in ob_generators(dom(F)))
    for o in ob_generators(apex(domprod))
        αcomps[ob_map(domleg1, o)][Symbol(ob_map(domleg2, o))] = collect(diagram_map(curr_res)[o])
    end
    FU = compose(F, ucodom)
    α = Dict(map(collect(αcomps)) do (o, comps)
      o => ACSetTransformation(ob_map(I,o), ob_map(FU, o); comps...)
    end)
    return DiagramHom{id}(F, α, I, ucodom)
end

"""
A left kan extension LanF(X) is an initial object in the category of extensions
of X along F (X & F viewed as morphisms in the category of diagrams w/ covariant
transformations).
     B              B
F ↗ η⇑ ↘ L      F ↗ α⇑ ↘ M
 A  ⟶  C        A  ⟶  C
    X
    F⋅L  !
  η ↗   ↘
  X  ⟶  F⋅M      where ∀ a ∈ A: η⋅! = α
     α
We don't currently store data in the chase for computing the universal property.
This data is the provenance of each freely added element (for which morphism f
in B did the f_total trigger cause this element to be created?). Instead, we use
an iterative algorithm to recover this information after the fact.
"""
function lk_universal(η::DiagramHom{id,D,CD}, α::DiagramHom{id,D,CD}
                  ) where {D,CD}
  ηα = [η, α]
  L, M = diagram.(codom.(ηα))
  B, C = [only(Set(get.([L, M]))) for get in [dom, codom]]
  X = only(Set(diagram.(dom.(ηα))))
  A = dom(X)
  Fs = force.([x.shape_map for x in ηα])
  F = first(Fs)
  all([F_==F for F_ in Fs]) || error("Fs must match")
  codom(X) == C || error("bad codom")

  # Constraints on components of the natural transformation we wish to compute
  init = Dict([bo => DefaultDict{Symbol,Dict{Int,Int}}(
    ()->Dict{Int,Int}()) for bo in Symbol.(ob_generators(B))])
  # Each element a induces constraints
  for a in ob_generators(A)
    # We don't know ahead of time how many paths F(a)->b we will need to explore
    # before we no longer get any new information, so we maintain a stack.
    stack = [(Symbol(ob_map(F,a)), α.diagram_map[a], η.diagram_map[a])]
    while !isempty(stack)
      # Find constraints on a homomorphism codom(η[a]) ---> codom(α[a])
      # or, given : F(a)->b, we have constraints for η[a];F(f) and α[a];F(f)
      b,f,g = pop!(stack)
      changed = false
      for (ob, mapping) in pairs(components(f))
        for i in parts(codom(g), ob)
          for v in Set(mapping(preimage(g[ob], i)))
            if haskey(init[b][ob],i)
              init[b][ob][i] == v || error("Inconsistent")
            else
              init[b][ob][i] = v
              changed |= true
            end
          end
        end
      end
      if changed
        for h in hom_generators(B)
          if Symbol(dom(h)) == b
            push!(stack, (Symbol(codom(h)), f⋅hom_map(M,h),g⋅hom_map(L,h)))
          end
        end
      end
    end
  end

  σf = Dict(map(Symbol.(ob_generators(B))) do bo
    i = NamedTuple(init[bo])
    hs = homomorphisms(ob_map(L,bo), ob_map(M,bo); initial=i)
    Symbol(bo) => only(hs)
  end)

  res = FinTransformation(L, M; σf...)
  # is_natural(res; verbose=true) || error("res $res")
  return res
end

# (co)Limits in Category of Diagrams
####################################

"""
Product of diagrams in the same category. (coproduct for :op morphisms)
"""
function diagram_hom_product(Xs::AbstractVector{<: Diagram{T}}; kw...) where T

  # Collect / validate type information about the input
  cod = codom(diagram(first(Xs))) # the category the diagrams are in
  is_id = T == id
  is_id || T == op || error("Diagrams of type $T not supported")
  # Equality of typecat too strict, but in theory this should be checked
  # cods = Set([codom(diagram(x)) for x in Xs])
  # length(cods) == 1 || error("Can't take product of diagrams in different cats")

  # Collect data about the product of the shape categories
  P = product([dom(diagram(x)) for x in Xs])
  obs, homs = map([ob_map=>ob_generators, hom_map=>hom_generators]) do (F, gen)
   [tuple([F(l, g) for l in legs(P)]...) for g in gen(apex(P))]
  end;

  # Take (co)products of objects in the underlying category
  omap, nat, base_prod = Dict(), Dict(), Dict()
  for os in obs
    lim = is_id ? product : coproduct
    base_prod[os] = p = lim([ob_map(diagram(X), o) for (o, X) in (zip(os, Xs))])
    s = Symbol(os)
    omap[s] = apex(p)
    nat[s] = legs(p)
  end

  hmap = Dict(map(homs) do hs
    maps = [hom_map(diagram(X), o) for (o, X) in (zip(hs, Xs))]
    Symbol(hs) => (is_id ? otimes : oplus)(maps)
  end)

  # Assemble product diagram
  apx = Diagram{T}(FinDomFunctor(omap, hmap, apex(P), cod))
  ls = map(enumerate(zip(legs(P), Xs))) do (i, (l, X))
    η = Dict([k=>v[i] for (k,v) in nat])
    src, tgt = is_id ? (apx, X) : (X, apx)
    DiagramHom{T}(l, η, src, tgt)
  end;
  return P, base_prod, Vector{DiagramHom{T}}(ls)
end

"""
Computes the equalizer (for id morphisms) / coequalizer (for op morphisms)
"""
function diagram_hom_equalizer(fs::AbstractVector{<: DiagramHom{T}}) where T
  # Collect / validate type information about the input
  X, Y = diagram.(only.(Set.([dom.(fs), codom.(fs)])));
  is_id = T == id
  is_id || T == op || error("diagram hom of type $T not supported")
  eq, eq_incl = is_id ? (equalizer, incl) : (coequalizer, proj)
  cod = codom(diagram(dom(first(fs))))
  # Equality of typecat too strict, but in theory this should be checked
  # cods = Set([codom(diagram(x)) for x in Xs])
  # cod = only(Set(vcat([codom.(diagram.([dom(x),codom(x)])) for x in fs]...)))

  # Shape level equalizer
  shape_eq = equalizer(shape_map.(fs))
  Eshape = apex(shape_eq)

  # Underlying (co)equalizer on natural transformations
  eqs = Dict(map(ob_generators(Eshape)) do o
    o=>eq([diagram_map(f)[o] for f in fs])
  end)
  η = Dict([o=>eq_incl(e) for (o,e) in collect(eqs)])

  om_ = Dict([o=>apex(e) for (o,e) in collect(eqs)])
  hm_ = Dict(map(hom_generators(Eshape)) do h
    h => compose(η[dom(h)], hom_map(is_id ? X : Y, h)) # TO CONFIRM: don't we need to factorize with η[codom(h)]?
  end)

  # Assemble diagram morphism
  E = Diagram{T}(FinDomFunctor(om_,hm_,Eshape, cod))
  src, tgt = is_id ? (diagram(E), X) : (Y, diagram(E))
  l = DiagramHom{T}(incl(shape_eq), η, src, tgt)
  return shape_eq, eqs, l
end


"""
Coproduct of diagrams in the same category. (product for op morphisms)
"""
function diagram_hom_coproduct(Xs::AbstractVector{<: Diagram{T}}; kw...) where {T}
  # Collect / validate type information about the input
  is_id = T == id
  is_id || T == op || error("Diagrams of type $T not supported")
  cod = codom(diagram(first(Xs))) # the category the diagrams are in
  # Equality of typecat too strict, but in theory this should be checked
  # cod = only(Set([codom(diagram(x)) for x in Xs]))

  # Shape level coproduct
  hcprod = coproduct(FinCatPresentation[dom(diagram(x)) for x in Xs])

  # Inclusion of data in underlying category
  obs, homs = Dict(), Dict()
  for (l, X) in zip(legs(hcprod), Xs)
    for (k,v) in ob_map(diagram(X))
      obs[ob_map(l, k)] = v
    end
    for (k,v) in hom_map(diagram(X))
      homs[hom_map(l, k)] = v
    end
  end

  # Assemble diagram morphism
  apx = Diagram{T}(FinDomFunctor(obs,homs, apex(hcprod), cod))
  ls = map(zip(legs(hcprod), Xs)) do (l, X)
    eta = Dict(k => id(v) for (k, v) in ob_map(diagram(X)))
    src, tgt = diagram.(is_id ? [X, apx] : [apx, X])
    DiagramHom{T}(l, eta, src, tgt)
  end;
  return hcprod, Dict(), ls
end

"""
Coequalizer of diagrams in the same category. (equalizer for op morphisms)
"""
function diagram_hom_coequalizer(fs::AbstractVector{<: DiagramHom{T}}; kw...) where {T}
  X, Y = diagram.(only.(Set.([dom.(fs), codom.(fs)])));
  h = Symbol(string(hash(fs))[1:4])
  # Shape level coequalizer
  shape_ceq = coequalizer(shape_map.(fs))
  H = proj(shape_ceq)
  # Left kan extensions
  κ = leftkan(H, Y, Symbol("H$h"); verbose=false)
  FH = shape_map(first(fs))⋅H # F⋅H = G⋅H = ...
  λ = leftkan(FH, X, Symbol("HF$h"); verbose=false)
  αs = [lk_universal(λ, ϕ⋅κ) for ϕ in fs]
  # Take coequalizer of universals, but need to convert them to CSet morphisms
  # TODO: work out doing coequalizers pointwise without the currying/uncurrying
  # This is important to return the coequalizer as a dictionary,
  # which is important for implementing the universal property of coequalizers
  curried = curry.(αs)
  csettypes = [AnonACSetType(presentation(dom(dom(x)))#=, Symbol("x$h")=#)
               for x in curried]
  csethoms = [αtype(αc) for (αtype, αc) in zip(csettypes, curried)]
  ceq = coequalizer(csethoms)
  γ = uncurry(FinTransformationMap(proj(ceq)), first(αs))
  KX = Diagram{T}(codom(γ))
  κγ = Dict([o=> compose(f, components(γ)[Symbol(ob_map(H, o))])
             for (o,f) in components(diagram_map(κ))])
  l = DiagramHom{T}(H, κγ, Diagram{T}(Y), KX)
  return shape_ceq, Dict(nothing=>ceq), l
end

@auto_hash_equals struct DiagLimit{
  Ob,Diagram,LimType<:Union{Limit,Colimit},Cone<:Multispan{Ob}
    } <: AbstractLimit{Ob,Diagram}
  diagram::Diagram
  shapelim::LimType
  baselim::Dict
  cone::Cone
end


@auto_hash_equals struct DiagColimit{
  Ob,Diagram,LimType<:Union{Limit,Colimit},Cocone<:Multicospan{Ob}
    } <: AbstractColimit{Ob,Diagram}
  diagram::Diagram
  shapelim::LimType
  baselim::Dict
  cocone::Cocone
end

product(Xs::AbstractVector{<: Diagram{id}}; kw...) =
  let r = diagram_hom_product(Xs)
  DiagLimit(DiscreteDiagram(Xs), r[1], r[2], Multispan(r[3])) end

coproduct(Xs::AbstractVector{<: Diagram{id}}; kw...) =
  let (sl, dl, csp) = diagram_hom_coproduct(Xs)
  DiagColimit(DiscreteDiagram(Xs), sl, dl, Multicospan(csp)) end

product(Xs::AbstractVector{<: Diagram{op}}; kw...) =
  let r = diagram_hom_coproduct(Xs)
  DiagLimit(DiscreteDiagram(Xs), r[1], r[2], Multispan(r[3])) end

coproduct(Xs::AbstractVector{<: Diagram{op}}; kw...) =
  let r = diagram_hom_product(Xs)
  DiagColimit(DiscreteDiagram(Xs), r[1], r[2], Multicospan(r[3])) end

equalizer(fs::AbstractVector{<: DiagramHom{id}}) =
  let r = diagram_hom_equalizer(fs)
  DiagLimit(ParallelMorphisms(fs), r[1], r[2], Multispan([r[3]])) end

coequalizer(fs::AbstractVector{<: DiagramHom{op}}) =
  let r = diagram_hom_equalizer(fs)
  DiagColimit(ParallelMorphisms(fs), r[1], r[2], Multicospan([r[3]])) end

coequalizer(fs::AbstractVector{<: DiagramHom{id}}) =
  let r = diagram_hom_coequalizer(fs)
  DiagColimit(ParallelMorphisms(fs), r[1], r[2], Multicospan([r[3]])) end

function pullback(fs::Multicospan{<:Diagram{T}, <: DiagramHom{T}}) where {T}
  p = product(Diagram{T}[dom(f) for f in fs])
  equalizer(DiagramHom{T}[compose(l, f) for (l,f) in zip(legs(p), fs)])
end

function pushout(fs::Multispan{<:Diagram{T}, <: DiagramHom{T}}) where {T}
  cp = coproduct(Diagram{T}[codom(f) for f in fs])
  coequalizer(DiagramHom{T}[compose(f,l) for (l,f) in zip(legs(cp), fs)])
end

function universal(p::DiagLimit{<:Diagram{id},<:DiscreteDiagram}, sp::Multispan)
  a_p, a_sp = apex.([p, sp])
  s_map = universal(p.shapelim, Multispan(shape_map.(legs(sp))))
  d_map = Dict(map(ob_generators(dom(diagram(a_sp)))) do o
    d_tgts = [diagram_map(l)[o] for l in sp]
    tgts = tuple([ob_map(shape_map(l),o) for l in sp]...)
    o => universal(p.baselim[tgts], Multispan(d_tgts))
  end)
  DiagramHom{id}(s_map, d_map, a_sp, a_p)
end


function universal(p::DiagLimit{<:Diagram{op},<:DiscreteDiagram}, sp::Multispan)
  a_p, a_sp = apex.([p, sp])
  s_map = universal(p.shapelim, Multicospan(shape_map.(legs(sp))))
  d_map = Dict()
  for (spl, pl) in zip(legs(sp),legs(p))
    for (k,v) in components(diagram_map(spl))
      d_map[ob_map(shape_map(pl),k)] = v
    end
  end
  DiagramHom{op}(s_map, d_map, a_sp, a_p)
end

function universal(cp::DiagColimit{<:Diagram{id},<:DiscreteDiagram}, csp::Multicospan)
  a_cp, a_csp = apex.([cp, csp])
  s_map = universal(cp.shapelim, Multicospan(shape_map.(legs(csp))))
  d_map = Dict()
  for (cspl, cpl) in zip(legs(csp),legs(cp))
    for o in ob_generators(dom(diagram(dom(cpl))))
      d_map[Symbol(ob_map(shape_map(cpl), o))] = diagram_map(cspl)[o]
    end
  end
  DiagramHom{id}(s_map, d_map, a_cp, a_csp)
end

function universal(p::DiagColimit{<:Diagram{op},<:DiscreteDiagram}, sp::Multicospan)
  a_p, a_sp = apex.([p, sp])
  s_map = universal(p.shapelim, Multispan(shape_map.(legs(sp))))
  d_map = Dict(map(ob_generators(dom(diagram(a_sp)))) do o
    d_tgts = [diagram_map(l)[o] for l in sp]
    tgts = tuple([ob_map(shape_map(l),o) for l in sp]...)
    o => universal(p.baselim[tgts], Multicospan(d_tgts))
  end)
  DiagramHom{op}(s_map, d_map, a_p, a_sp)
end

function uni_eq(eq, fs)
  f = only(fs)
  a_eq = apex(eq)
  s_map = universal(eq.shapelim, shape_map(f))
  d_map = Dict(map(collect(ob_map(shape_map(f)))) do (k, v)
    k => factorize(eq.baselim[v], diagram_map(f)[k])
  end)
  s_map, d_map, a_eq
end

function universal(eq::DiagLimit{<:Diagram{id},<:ParallelMorphisms}, fs::SMultispan{1})
  s_map, d_map, a_eq = uni_eq(eq, fs)
  DiagramHom{id}(s_map, d_map, dom(only(fs)), a_eq)
end

function universal(eq::DiagColimit{<:Diagram{op},<:ParallelMorphisms}, fs::SMulticospan{1})
  s_map, d_map, a_eq = uni_eq(eq, fs)
  DiagramHom{op}(s_map, d_map, a_eq, codom(only(fs)))
end

end