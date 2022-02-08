using Revise
using AlgebraicPetri
using OrdinaryDiffEq
using DiffEqFlux, Flux
using Catalyst
using ModelExploration

SIR =  LabelledPetriNet([:S, :I, :R],
                        :inf => ((:S, :I)=>(:I, :I)),
                        :rec => (:I=>:R));

initial_conc=[1e7,1e5,1e-10]
true_rates = [1e-8,1e-4]
data = generate_data(SIR, true_rates,initial_conc)


initial_states = initial_conc;
p = [-10.,-10.,-10.];
g = SIR;
g=SIR

prob = ODEProblem(ReactionSystem(g),
                    zeros(ns(g)), # initial state
                    (0.0,100.0),
                    zeros(ns(g)+nt(g)));

cur_p = exp.(p)
prob′ = remake(prob, p=cur_p, u0=initial_states, tspan=(0.0,100.0))
sol = solve(prob′, Tsit5())
sol_steps = [sol(i) for i in 1:100]
prediction = [sum_infected_states(g, s) for s in sol_steps]
sum(abs2, data .- prediction), sol

#x = train(SIR, data)

function test_train(g::LabelledPetriNet, data::Vector{Float64})
  param_guess = [1e-7,1e-5,1e7,1e5,1e-10] # correct [1e-8,1e-4,1e7,1e5,1e-10]
  prob = ODEProblem(ReactionSystem(g),
                    zeros(ns(g)), # initial state
                    (0.0,100.0),
                    zeros(ns(g)+nt(g))) # current variables in model

  l_func = make_loss(g, prob, data, states=Dict([1=>1e7,2=>1e5,3=>1e-10]),
                      rates=Dict([2=>1e-5]))

  p = DiffEqFlux.sciml_train(l_func,log.(param_guess),ADAM(0.1),maxiters = 1000)
  return p, l_func
end

x, l = test_train(SIR, data)

