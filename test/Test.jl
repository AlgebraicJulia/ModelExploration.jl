using ModelExploration
using AlgebraicPetri
using Test
using Catlab.CategoricalAlgebra


# Petri net example
###################

# Slices over A are PetriNets that have only unary and binary transitions
A = @acset PetriNet begin
    S = 4
    T = 2
    I = 3
    O = 3
    it = [1,2,2]
    is = [1,2,3]
    ot = [1,2,2]
    os = [2,2,3]
end

"""
Create linear Petri Net
Could we do this with functor from Graph?
"""
function multi_city(n::Int)::StructACSet
    res = PetriNet()
    add_parts!(res, :S, n)
    for sym in [:T, :I, :O]
        add_parts!(res, sym, n-1)
    end
    set_subparts!(res, :it, 1:n-1)
    set_subparts!(res, :ot, 2:n)
    set_subparts!(res, :is, 1:n-1)
    set_subparts!(res, :os, 2:n)
    return res
end


MultiCity = Linear(multi_city)

AgeRefine = Finite([young_old, young_mid_old, rich_poor, youn_old_rich_poor],
                   [1=>2,1=>4,3=>4])

exp = Explore(PetriNet, A, SIR, [Free(), MultiCity, AgeRefine])

data = DataSet(path)

exp.run(data)


