module TestPetri
using ModelExploration

import AlgebraicPetri
using Test
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.Graphs.SimpleGraphs: TheoryGraph

const Petri = AlgebraicPetri.PetriNet

"""
Loose ACSet Transformations

Make a wrapper for Slice and implement limits for that.
"""

# Petri net example
###################
function dezero(s::StructACSet)::StructACSet
  for h in keys(s.homs)
    set_subpart!(s, h, replace(x -> x == 0 ? 1 : x, s[h]))
  end
  return s
end

function graph_to_petri(G::Graph)::Slice
  idE = id(TheoryGraph[:E])
  obs = Dict(:S => :V, :T=>:E, :I => :E, :O=>:E)
  morphs = Dict(:is => :src, :it => idE, :os => :src, :ot => idE)
  return Slice(migrate(PetriNet, G, obs, morphs), A, [:S=>[1=>1]])
end

# Slices over A are PetriNets that have only unary and binary transitions
UnaryT = dezero(@acset Petri begin S = 1; T = 1; I = 1; O = 1 end)
BinaryT = dezero(@acset Petri begin S = 1; T = 1; I = 2; O = 2 end)
S1 = @acset Petri begin S=1 end
A = apex(pushout(homomorphism(S1, UnaryT), homomorphism(S1, BinaryT)))

"""
Create linear Petri Net
Could we do this with functor from Graph?
"""
multi_city(n::Int)::Slice = graph_to_petri(complete_graph(Graph, n))


MultiCity = Linear(multi_city)

young_old, rich_poor = multi_city(2), multi_city(2)
young_mid_old = multi_city(3)
young_old_rich_poor = multi_city(4)

AgeRefine = Finite([young_old, young_mid_old, rich_poor, young_old_rich_poor],
                    [[2,4],    Int[],         [4],       Int[]])

exp = Explore(PetriNet, [Slice(A)], SIR, [Free(), MultiCity, AgeRefine])

data = DataSet(path)

exp.run(data)


end # module