using Revise
using Test
using ModelExploration
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories
include("../src/Petri.jl")
using AlgebraicPetri: LabelledPetriNetUntyped

infectious_type_ = LabelledPetriNet([:Pop],
  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :t_infection=>(:Pop=>:Pop),
  :t_strata=>(:Pop=>:Pop))
#Graph(infectious_type)
infectious_type = map(infectious_type_, Name=name->nothing)

IO_help(i::Int) = let d = Dict([j=>j%2+1 for j in 1:2*i]); (I=d,O=d) end
make_slice(p::LabelledPetriNet, n::NamedTuple) =
    Slice{LabelledPetriNet,ACSetTransformation}(
        homomorphism(p, infectious_type; initial=n,
                       type_components=(Name=x->nothing,)))

SIR = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
  :id => (:S => :S),
  :id => (:I => :I),
  :id => (:R => :R)
)
SIR_type = make_slice(SIR, merge((T=[1,2,3,3,3],), IO_help(1)))
#
SIR_travel_ban = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
  :id => (:S => :S),
  :id => (:R => :R)
)
SIR_tb_type = make_slice(SIR_travel_ban, merge((T=[1,2,3,3],), IO_help(1)))
#
SIS = LabelledPetriNet([:S, :I],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:S),
  :id => (:S => :S),
  :id => (:I => :I),
)
SIS_type = make_slice(SIS, merge((T=[1,2,3,3],),IO_help(1)))
#
SVIIvR = LabelledPetriNet([:S, :I, :R, :Iv, :V],
  :inf => ((:S, :I)=>(:I, :I)),
  :inf => ((:V,:I)=>(:I,:I)),
  :inf => (:S, :Iv) => (:Iv,:Iv),
  :inf => ((:V, :Iv) => (:Iv,:Iv)),
  :rec => (:I=>:R),
  :rec => (:Iv=>:R),
  :vax => (:S=>:V),
  :id => (:S => :S),
  :id => (:I => :I),
  :id => (:V => :V),
  :id => (:Iv => :Iv),
  :id => (:R => :R)
)

SVIIvR_type = make_slice(SVIIvR, merge((T=[1,1,1,1,2,2,2,3,3,3,3,3],),IO_help(4)))
#
quarantine = LabelledPetriNet([:Q, :not_Q],
    :interaction => ((:not_Q, :not_Q) => (:not_Q, :not_Q)),
    :enter_quarantine => (:not_Q => :Q),
    :exit_quarantine => (:Q => :not_Q),
    :id => (:Q => :Q),
    :id => (:not_Q => :not_Q),
)
quarantine_type = make_slice(quarantine, (T=[1, 3, 3, 2, 2], I=Dict(1=>1,2=>2), O=Dict(1=>1,2=>2)))
#
age_stratification = LabelledPetriNet([:Child, :Adult],
    :interaction => ((:Child, :Child) => (:Child, :Child)),
    :interaction => ((:Adult, :Adult) => (:Adult, :Adult)),
    :interaction => ((:Child, :Adult) => (:Child, :Adult)),
    :id => (:Child => :Child),
    :id => (:Adult => :Adult),
)
age_s_type = make_slice(age_stratification, merge((T=[1,1,1,2,2],),IO_help(3)))

#
flux_metapopulation = LabelledPetriNet([:Patch1, :Patch2],
    :interaction => ((:Patch1, :Patch1) => (:Patch1, :Patch1)),
    :interaction => ((:Patch2, :Patch2) => (:Patch2, :Patch2)),
	:travel => (:Patch1 => :Patch2),
    :travel => (:Patch2 => :Patch1),
    :id => (:Patch1 => :Patch1),
    :id => (:Patch2 => :Patch2),
)
flux_m_type = make_slice(flux_metapopulation, merge((T=[1,1,3,3,2,2],), IO_help(2)))
simple_trip =  LabelledPetriNet([:P11, :P21, :P12, :P22],
    :interaction11_11 => ((:P11, :P11) => (:P11, :P11)),
    :interaction12_12 => ((:P12, :P12) => (:P12, :P12)),
    :interaction21_21 => ((:P21, :P21) => (:P21, :P21)),
    :interaction22_22 => ((:P22, :P22) => (:P22, :P22)),
    :interaction11_21 => ((:P11, :P21) => (:P11, :P21)),
    :interaction22_12 => ((:P22, :P12) => (:P22, :P12)),
    :interaction12_22 => ((:P12, :P22) => (:P12, :P22)),
    :interaction21_11 => ((:P21, :P11) => (:P21, :P11)),
	:travel11_12 => (:P11 => :P12),
	:travel12_11 => (:P12 => :P11),
	:travel21_22 => (:P21 => :P22),
	:travel22_21 => (:P22 => :P21),
    :id => (:P11 => :P11),
    :id => (:P21 => :P21),
    :id => (:P12 => :P12),
    :id => (:P22 => :P22),
)
simple_t_type = make_slice(simple_trip,
    merge((T=[1,1,1,1,1,1,1,1,3,3,3,3,2,2,2,2],), IO_help(8)))

@present ThFour(FreeCategory) begin
    (X1,X2,X3,X4)::Ob
end
@present ThOne(FreeCategory) begin
  X::Ob
end

const ACSetCat{S} = TypeCat{S, ACSetTransformation}
const ACSetCatSlice{S} = TypeCat{S, ACSetTransformation}
const Petri = ACSetCat{LabelledPetriNet}
const PetriHom = SliceCat{LabelledPetriNet,ACSetTransformation}

One, Four = FinCat.([ThOne,ThFour])

to_slicehom(x) = SliceHom(Literal(Diagram(
  FinDomFunctor(x, nothing, Four, PetriHom()))))

diag_disease = to_slicehom(Dict(:X1=>SIR_type,:X2=>SIS_type,
                               :X3=>SVIIvR_type,:X4=>SIR_tb_type));
diag_strata = to_slicehom(Dict(:X1=>quarantine_type, :X2=>age_s_type,
                               :X3=>flux_m_type, :X4=>simple_t_type));

mh1 = to_model_hom(diag_disease)
mh2 = to_model_hom(diag_strata)
pb = PullbackSpace(mh1, mh2);
upb = unfold(pb);
aupb = apex(upb);
@test length(ob_generators(dom(diagram(aupb)))) == 16

x = ob_map(aupb, Symbol("(X1, X1)"))


