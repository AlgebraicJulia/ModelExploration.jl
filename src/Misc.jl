module Misc

using Catlab.CategoricalAlgebra, Catlab.Present
import Catlab.CategoricalAlgebra: limit, universal
import Catlab.CategoricalAlgebra.Categories: is_hom_equal
import Base: hash


is_hom_equal(f::ACSetTransformation, g::ACSetTransformation) =
  force(f) == force(g)

function limit(ps::ParallelMorphisms{<:TypeSet})
    err = "equalizers of TypeSets that are not identity not supported"
    eltype(codom(ps)) == Nothing || error(err)
    d = dom(ps)
    Limit(ps, Multispan(d, [IdentityFunction(d, d)]))
  end

universal(ps::Equalizer{<:TypeSet}, span::Multispan) = only(span) ⋅ incl(ps)

Base.hash(pres::Presentation{T,N}, h::UInt) where {T,N} =
  hash(T, hash(N, hash(pres.syntax, hash(pres.generators,
       hash(pres.equations, h)))))

end