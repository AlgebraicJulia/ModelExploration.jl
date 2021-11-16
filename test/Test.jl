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


exp = Explore(PetriNet, A, SIR, [Free(), MultiCity, AgeRefine])

data = DataSet(path)

exp.run(data)


