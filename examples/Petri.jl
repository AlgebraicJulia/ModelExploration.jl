using AlgebraicPetri
using OrdinaryDiffEq
using DiffEqFlux
#using Catalyst
using LinearAlgebra: dot

"""
We will eventually be more sophisticated and use the actual projection
maps of a pullback, but for now we check whether "I" is in the name

Returns vector int
"""
get_infected_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("I", string(s))]

function get_bounds(g::AbstractLabelledPetriNet)
  lb, ub = Float64[],Float64[]
  for (sym,_) in g[:tname]
    sym_= string(sym)
    if occursin("inf", sym_)
      push!(lb, log(1e-8)); push!(ub, log(1e-5))
    elseif occursin("rec", sym_)
      push!(lb, log(1e-4)); push!(ub, log(1e-1))
    else
      push!(lb, log(1e-10)); push!(ub, log(1))
    end
  end
  for (sym,_) in g[:sname]
    sym_ = string(sym)
    if occursin("I", sym_)
      push!(lb, log(1e2)); push!(ub, log(1e6))
    elseif occursin("R", sym_)
      push!(lb, log(1e-10)); push!(ub, log(1))
    else
      push!(lb, log(1e2)); push!(ub, log(1e8))
    end
  end
  return lb => ub
end
"""
state - Vector{Float}
"""
sum_infected_states(g::AbstractLabelledPetriNet, states) =
  sum([s for (i, s) in enumerate(states) if i ∈ get_infected_states(g)])



"""Generate a sequence of new infections"""
function generate_data(rxn::AbstractLabelledPetriNet, rxn_rates, initial_conc)
  o = ODEProblem(ReactionSystem(rxn),  initial_conc,
                 (0.0,100.0), vcat(rxn_rates, initial_conc))
  trajectory = solve(o, Tsit5())
  t_steps = [trajectory(i) for i in 0:100]
  return [sum_infected_states(rxn, t) for t in t_steps]
end


"""
Produces function from Vector{Float} -> Float.
- The input is the logged params (exponentiate to get the real params)

Outer function inputs:
prob - current problem (ODE problem)
times - each timestep (indices of timestep correspond to indices of data array)
data - Vector{Float} representing Infected at each step - target to hit
states/rates - dictionary whatever you want to be fixed
"""
function make_loss(g::AbstractLabelledPetriNet, prob, data; states=Dict(), rates=Dict())
  function loss(p)
    cur_p = exp.(p)
    u0 = exp.(p[(1:ns(g)) .+ nt(g)])
    for (k,v) in rates
      cur_p[k] = v
    end
    for (k,v) in states
        u0[k] = v
    end

    prob′ = remake(prob, p=cur_p, u0=u0, tspan=(0.0,100.0))
    sol = solve(prob′, Tsit5())
    sol_steps = [sol(i) for i in 0:100]
    prediction = [sum_infected_states(g, s) for s in sol_steps]
    log(sum(abs2, data .- prediction)), sol
  end
end

"""
Train a Petri net against data. Return the optimized parameters and minimum loss.
"""
function train(g::AbstractLabelledPetriNet,data::Vector{Float64})::Pair{Float64, Vector{Float64}}
  param_guess = vcat(repeat([1e-10], nt(g)), repeat([1e-3], ns(g)))
  #param_guess = [1e-10,1e-10,1e6,1e6,1e6] # correct [1e-8,1e-4,1e7,1e5,1e-10]
  prob = ODEProblem(ReactionSystem(g),
                    zeros(ns(g)), # initial state
                    (0.0,100.0),
                    zeros(ns(g)+nt(g))) # current variables in model
  lb, ub = get_bounds(g)
  l_func = make_loss(g, prob, data)
  p = DiffEqFlux.sciml_train(l_func,lb,ADAM(0.1), lower_bounds=lb, upper_bounds=ub, maxiters = 10000)
  return l_func(p.u)[1] => exp.(p.u)
end

function eval_petri_fn(true_net::AbstractLabelledPetriNet, true_rates::Vector{Float64}, true_initial_conc::Vector{Float64})::Function
  data = generate_data(true_net, true_rates, true_initial_conc)
  function model_search_loss_fn(g::AbstractLabelledPetriNet)#::Float64
    return first(train(g, data))
    #return train(g, data)
  end
  return model_search_loss_fn
end

#l_func = make_loss(model, prob, week_times, week_avg; calc_inf=calc_inf, rates=Dict(5=>1/14, 6=>1/14), states=Dict(2=>1e-9, 4=>1e-9));

#p = DiffEqFlux.sciml_train(l_func,log.(param_guess),ADAM(0.1),maxiters = 1000)