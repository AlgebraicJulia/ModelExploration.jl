using Revise
using Test
using ModelExploration
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories
#include("../src/Petri.jl")
include("../src/ModelSelection.jl")
using AlgebraicPetri: LabelledPetriNetUntyped

infectious_type_ = LabelledPetriNet([:Pop],
  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :t_infection=>(:Pop=>:Pop),
  :t_strata=>(:Pop=>:Pop))
#Graph(infectious_type)
infectious_type = map(infectious_type_, Name=name->nothing)

IO_help(i::Int) = let d = Dict([j=>j%2+1 for j in 1:2*i]); (I=d,O=d) end
make_slice(p::LabelledPetriNet, n::NamedTuple) =
    Slice{ACSetTransformation}(
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
  :inf => ((:V,:I)=>(:Iv,:I)),
  :inf => (:S, :Iv) => (:I,:Iv),
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
const PetriHom = SliceCat{ACSetTransformation}

One, Four = FinCat.([ThOne,ThFour])

to_diag(x) = Diagram(FinDomFunctor(x, nothing, Four, PetriHom()))
to_slicehom(x) = SliceDiagHom(Literal(to_diag(x)))

diag_disease = to_slicehom(Dict(:X1=>SIR_type,:X2=>SIS_type,
                               :X3=>SVIIvR_type,:X4=>SIR_tb_type));
diag_strata = to_slicehom(Dict(:X1=>quarantine_type, :X2=>age_s_type,
                               :X3=>flux_m_type, :X4=>simple_t_type));

pb = PullbackSpace(to_model_hom(diag_disease), to_model_hom(diag_strata));
upb = unfold(pb);
aupb = apex(upb);
@test length(ob_generators(dom(diagram(aupb)))) == 16

true_model = ob_map(aupb, Symbol("(X3, X1)"))
true_rxn = MakeReactionSystem(true_model)

p_real = vcat(repeat([1e-4], 14), repeat([0.01], 6))
u0 = [0.0,0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0,0.0] #|> gpu

#=true_model = ob_map(aupb, Symbol("(X2, X2)"))
true_rxn = MakeReactionSystem(true_model)
p_real = [2e-4,2e-4,2e-4,4e-2,6e-2]
u0 = [500.0, 0.0, 499.0, 1.0]=#
tspan = (0.0,100.0)

#sample_data, sample_times, prob_real, sol_real = generate_data(true_model, p_real, u0, tspan, 50)=#

#=true_model = SIR
true_rxn = MakeReactionSystem(SIR)
p_real = [.1/1000, .01]
u0 = [999.0, 1.0, 0.0]
tspan = (0.0,250.0)=#

sample_data, sample_times, prob_real, sol_real = generate_data(true_model, p_real, u0, tspan, 50)

plt = plot(sol_real, lw=2, label=reshape(map(string, true_model[:, :sname]), 1, ns(true_model)))
plot!(sample_times, sample_data[:,1:2:3], seriestype=:scatter, label="")
display(plt)

minimize_Imax(true_model, u0, tspan, sample_times, 
  rates=Dict(1=>1e-3, 2=>1e-3, 3=>1e-3, 4=>1e-3, 15=>1e-3, 16=>1e-3, 18=>1e-3, 19=>1e-3), target=250)

#small training example to precompile everything
#full_train(SIR, [999.0,1.0,0.0], tspan, sample_data, sample_times);

#assign inital concentrations to models TODO: be more smart about this
models = [
  ob_map(aupb, Symbol("(X1, X1)")) => [0.0, 0.0, 0.0, 999.0, 1.0, 0.0],
  ob_map(aupb, Symbol("(X1, X2)")) => [500.0, 0.0, 0.0, 499.0, 1.0, 0.0],
  # ob_map(aupb, Symbol("(X1, X3)")) => [500.0, 1.0, 0.0, 499.0, 0.0, 0.0],
  # ob_map(aupb, Symbol("(X1, X4)")) => [250.0, 1.0, 0.0, 250.0, 0.0, 0.0, 250.0, 0.0, 0.0, 249.0, 0.0, 0.0],
  # ob_map(aupb, Symbol("(X2, X1)")) => [0.0, 0.0, 999.0, 1.0],
  # ob_map(aupb, Symbol("(X2, X2)")) => [500.0, 0.0, 499.0, 1.0],
  # ob_map(aupb, Symbol("(X2, X3)")) => [500.0, 1.0, 499.0, 0.0],
  # ob_map(aupb, Symbol("(X2, X4)")) => [250.0, 1.0, 250.0, 0.0, 250.0, 0.0, 249.0, 0.0],
  ob_map(aupb, Symbol("(X3, X1)")) => [0.0, 0.0, 0.0, 0.0, 0.0, 999.0, 1.0, 0.0, 0.0, 0.0],
  ob_map(aupb, Symbol("(X3, X2)")) => [500.0, 0.0, 0.0, 0.0, 0.0, 499.0, 1.0, 0.0, 0.0, 0.0],
  # ob_map(aupb, Symbol("(X3, X3)")) => [500.0, 1.0, 0.0, 0.0, 0.0, 499.0, 0.0, 0.0, 0.0, 0.0],
  #ob_map(aupb, Symbol("(X3, X4)")) => [250.0, 1.0, 0.0, 0.0, 0.0,
  #                                     250.0, 0.0, 0.0, 0.0, 0.0,
  #                                     250.0, 0.0, 0.0, 0.0, 0.0,
  #                                     249.0, 0.0, 0.0, 0.0, 0.0],
  # ob_map(aupb, Symbol("(X4, X1)")) => [0.0, 0.0, 0.0, 999.0, 1.0, 0.0],
  # ob_map(aupb, Symbol("(X4, X2)")) => [500.0, 0.0, 0.0, 499.0, 1.0, 0.0],
  # ob_map(aupb, Symbol("(X4, X3)")) => [500.0, 1.0, 0.0, 499.0, 0.0, 0.0],
  # ob_map(aupb, Symbol("(X4, X4)")) => [250.0, 1.0, 0.0, 250.0, 0.0, 0.0, 250.0, 0.0, 0.0, 249.0, 0.0, 0.0]
]

function explore(models, tspan, sample_data, sample_times)
  losses = zeros(length(models))
  Threads.@threads for i in 1:length(models)
    model, u0 = models[i]
    sol, loss = full_train(model, u0, tspan, sample_data, sample_times)
    losses[i] = loss
    println(loss)
    #=plt = plot(sol, lw=2, label=reshape(map(string, model[:, :sname]), 1, ns(model)))
    plot!(sample_times, sample_data, seriestype=:scatter, label="")
    display(plt)=#
  end
  return losses
end

#=println(min_loss)
plt = plot(best_sol, lw=2, label=reshape(map(string, best_model[:, :sname]), 1, ns(best_model)))
plot!(sample_times, sample_data, seriestype=:scatter, label="")
display(plt)=#


# Pushout example
#=Infect = LabelledPetriNet([:I],)
Death = LabelledPetriNet([:I,:D],     :death => (:I => :D))
ISIR = only(homomorphisms(Infect, SIR))
IID = only(homomorphisms(Infect, Death))

I = Diagram(FinDomFunctor(Dict(:X=>Infect), nothing, One, Petri()))
D = FinDomFunctor(Dict(:X1=>SIR,:X2=>SIS,:X3=>SVIIvR,:X4=>SIR_travel_ban),
                  nothing, Four,Petri())
ID = LitModelHom(DiagramHom(FinFunctor(Dict(:X=>:X1), nothing, One, Four), Dict(:X=>ISIR),
                  I, Diagram(D)))
Q = Diagram(FinDomFunctor(Dict(:X=>Death), nothing, One, Petri()))
IQ = LitModelHom(DiagramHom(id(One),Dict(:X=>IID),I,Q))

#po = PushoutSpace(ID,IQ) # fails b/c no chase of ACSets only C-sets
                          # and no FinCats with Attr

# a pushout over a product space: (A+B*C) = A+B * A+C=#