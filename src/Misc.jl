module Misc

using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.Sets: IdentityFunction
import Catlab.CategoricalAlgebra: limit, universal
using Catlab.Present
using Catlab.Theories
import Base: hash

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