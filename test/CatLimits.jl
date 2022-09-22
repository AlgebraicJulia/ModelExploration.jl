module CatLimits

using Test
using ModelExploration.CatLimits
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.Present

# Define examples for limits and colimits
#----------------------------------------
@present I_(FreeSchema) begin (I1, I2)::Ob; i::Hom(I1, I2) end
@present M_(FreeSchema) begin (M1, M2)::Ob; m::Hom(M1, M2) end
@present N_(FreeSchema) begin (N1, N2)::Ob; n::Hom(N1, N2) end
@present J_(FreeSchema) begin
  (J1, J2, J3)::Ob; j1::Hom(J1, J2); j2::Hom(J3,J2)
end
@present K_(FreeSchema) begin
  K::Ob; k::Hom(K,K)
  compose(k,k) == k
end

I, J, K, M, N = FinCat.([I_,J_,K_,M_,N_]);
(i,),(j1,j2),(k_,) = hom_generators.([I,J,K])

F_IJ = FinDomFunctor(Dict(:I1=>:J1, :I2=>:J2), Dict(:i=>:j1), I, J)
G_IJ = FinDomFunctor(Dict(:I1=>:J3, :I2=>:J2), Dict(:i=>:j2), I, J)
F_IK = FinDomFunctor(Dict(:I1=>:K, :I2=>:K), Dict(:i=>:k), I, K)
J_I = FinDomFunctor(Dict(:J1=>:I1, :J2=>:I2, :J3 => :I1),
                    Dict(:j1=>:i, :j2=>:i), J, I)

# Limits
#-------
p = product([I, J, K]);
e = equalizer([F_IJ, G_IJ])
e2 = equalizer([F_IJ, F_IJ])

map([p, e, e2]) do lim
  @test all(is_functorial.(legs(lim)))
end

u = universal(p, Multispan([id(I), F_IJ, F_IK]))
@test is_functorial(u)
ij1k,i2j1k,i2j2k = [h for h in hom_generators(apex(p)) if h.args[1] in Symbol.(
  ["(i, id(J1), id(K))","(id(I2), j1, id(K))","(id(I2), id(J2), k)"])]
@test hom_map(u, :i) == compose(ij1k, i2j1k, i2j2k)

u = universal(e2, J_I)
@test is_functorial(u)

# encode a graph homomorphism as a FinCat. It is •→• ⟶ Grph
# by curry adjunction, the shape is the product of •→• and •⇉•
P = apex(product([I, FinCat(SchGraph)]))
@test length(equations(P)) == 2 # the two naturality squares
P = apex(product([I,M,N])) # multidimensional
@test length(equations(P)) == 6
P = apex(product([I,K])) # base equation in K gets multiplied by 2
@test length(equations(P)) == 3 # also have one naturality square


# Colimits
#---------
cp = coproduct([I, I, J]);
@test length(ob_generators(apex(cp))) == 7
cp2 = coproduct([K, K])
@test length(equations(apex(cp2))) == 2
ce1 = coequalizer([G_IJ, G_IJ]);
ce2 = coequalizer([F_IJ, G_IJ]);
map([cp, ce1, ce2]) do colim
  @test all(is_functorial.(legs(colim)))
end

u = universal(cp, Multicospan([F_IJ, G_IJ, id(J)]))
@test is_functorial(u)
@test ob_map(u, Symbol("I1#1")) == ob_generators(J)[1]
@test ob_map(u, Symbol("I1#2")) == ob_generators(J)[3]

u = universal(ce1, J_I)
@test is_functorial(u)

end