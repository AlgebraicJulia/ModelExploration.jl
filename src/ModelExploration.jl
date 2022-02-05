module ModelExploration

export Product, Literal, Gen, select

using Catlab.CategoricalAlgebra
using Catlab.Graphs

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

function select(g::Gen, lossFn::Function)::StructACSet
    return Graph()
end

end # module