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

# Not thinking too far ahead on this for now lol
function unfold(g::Literal)::Vector{StructACSet}
    return g.lit
end

function unfold(g::Product)::Vector{StructACSet}
    models = StructACSet[]
    for x in g.g1.unfold
        for y in g.g2.unfold
            for slice_x in homomorphisms(x, g.slice)
                for slice_y in homomorphisms(y, g.slice)
                    push!(models, ob(pullback(slice_x, slice_y)))
                end
            end
        end
    end
    return models
end

function select(g::Gen, lossFn::Function)::StructACSet
    return Graph()
end

end # module