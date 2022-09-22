module TestCurry

using Test
using Catlab.CategoricalAlgebra
using Catlab.CSetDataStructures: AnonACSetType
using Catlab.Graphs
using Catlab.Present
using Catlab.Theories
using ModelExploration.Curry

# FinFunctors
#############
const Grph = ACSetCat{Graph}
g1 = Graph(1)
ar = @acset Graph begin V=2; E=2; src=[1,2]; tgt=[2,2] end
t1 = apex(terminal(Graph))
t1_ar = homomorphism(t1, ar)
_, g1_arr2 = homomorphisms(g1, ar)

@present CSpanPres_(FreeSchema) begin
  (C1, C2, C3)::Ob; c1::Hom(C1, C2); c2::Hom(C3,C2)
end
CSpan = FinCat(CSpanPres_)

# Example FinFunctor into Grph
CG_t1ar = FinDomFunctor(Dict(:C1=>t1,:C2=>ar,:C3=>g1),
                        Dict(:c1=>t1_ar,:c2=>g1_arr2),
                        CSpan, Grph());

# FinFunctor to Grph -> FinFunctor to Set
cspan_graph_ex = curry(CG_t1ar);
# (reversible)
@test uncurry(cspan_graph_ex, CG_t1ar) == CG_t1ar

# FinFunctor to Set -> C-Set
cg = AnonACSetType(presentation(dom(cspan_graph_ex)))
cg_cset = cg(cspan_graph_ex)

# C-Set -> FinFunctor to Set -> FinFunctor to Grph
@test uncurry(FinDomFunctor(cg_cset), CG_t1ar) == CG_t1ar


end
