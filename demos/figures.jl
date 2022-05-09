using DataFrames
using CSV
using Plots
using ModelExploration
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories
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
#@test length(ob_generators(dom(diagram(aupb)))) == 16

sample_df = CSV.read("sample_data.csv", DataFrame)

#sample_times = sample_df.times

true_model = ob_map(aupb, Symbol("(X3, X1)"))
true_rxn = MakeReactionSystem(true_model)

p_real = vcat(repeat([1e-4], 14), repeat([0.01], 6))
p = [8.973917428097994e-5, 0.0001363771317481271, 0.0006662112292129558, 0.0006662112292129558, 1.0119613720380602e-6, 1.0000000000000002e-6, 0.004625656837867153, 0.0006662112292129555, 0.0006662112292129546, 0.0006840720342217968, 0.000667710343995489, 0.0004122138440619765, 0.0006662112292129555, 0.0006662112292129555, 0.0007520911121811695, 0.0006662112292129558, 0.00042694183846821024, 0.008117405177620019, 0.0006662112292129547, 0.010214821356644707]
u0 = [0.0,0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0,0.0]
tspan = (0.0,250.0)

sample_data, sample_times, prob_real, sol_real = generate_data(true_model, p_real, u0, tspan, 50)

prob = ODEProblem(true_rxn, u0, tspan, p_real)
sol = solve(prob, Tsit5(), tstops=sample_times)

susc_vals = map(sum, collect(zip([sol[i,:] for i in get_susceptible_states(true_model)]...)))
rec_vals = map(sum, collect(zip([sol[i,:] for i in get_recovered_states(true_model)]...)))
inf_vals = map(sum, collect(zip([sol[i,:] for i in get_infected_states(true_model)]...)))

SVIIvR_Q_df = CSV.read("SVIIvR_Q_traj.csv", DataFrame)

I_samples = sample_df.I_samples
S_samples = sample_df.S_samples
R_samples = sample_df.R_samples

plt = scatter(sample_times, sample_data)
#scatter!(sample_times, S_samples)
#scatter!(sample_times, R_samples)

I_vals = SVIIvR_Q_df.I_vals
S_vals = SVIIvR_Q_df.S_vals

#plot!(collect(0.0:250/(length(I_vals)-1):250.0), I_vals)
#plot!(collect(0.0:250.0/(length(S_vals)-1):250.0), S_vals)

plot!(collect(range(0.0,250.0, length=length(susc_vals))), susc_vals)
plot!(collect(range(0.0,250.0, length=length(inf_vals))), inf_vals)
display(plt)

plot(sol)