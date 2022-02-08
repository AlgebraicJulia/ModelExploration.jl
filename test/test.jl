using Revise
using AlgebraicPetri
using OrdinaryDiffEq
using DiffEqFlux, Flux
using Catalyst
using ModelExploration

SIR =  LabelledPetriNet([:S, :I, :R],
                        :inf => ((:S, :I)=>(:I, :I)),
                        :rec => (:I=>:R));

SIR2 =  LabelledPetriNet([:S, :I, :R],
                        :inf => (:S=>:I),
                        :rec => (:I=>:R));

initial_conc=[1e7,1e5,1e-10]
true_rates = [1e-8,1e-4]
data = generate_data(SIR, true_rates,initial_conc)

final_loss, fitted_params = train(SIR, data)

predicted_traj = generate_data(SIR, fitted_params[1:2],fitted_params[3:5])
plot(data; xaxis="days")
plot!(predicted_traj; xaxis="days")

f = eval_petri_fn(SIR, true_rates, initial_conc)
v1= f(SIR)
v2 = f(SIR2)