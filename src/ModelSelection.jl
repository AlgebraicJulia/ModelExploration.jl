using AlgebraicPetri
using Catlab.CategoricalAlgebra
using Catalyst
using OrdinaryDiffEq
using DiffEqFlux #=, Flux=#
using Plots
#import Catalyst: ReactionSystem

SIR = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
  :id => (:S => :S),
  :id => (:I => :I),
  :id => (:R => :R)
)

counter(a) = [count(==(i),a) for i in unique(a)]
function MakeReactionSystem(pn::AbstractPetriNet, name::Symbol)
    @parameters t k[1:(nt(pn) + ns(pn))]
    @variables S[1:ns(pn)](t)

    rxs = map(1:nt(pn)) do t
      inpts = pn[incident(pn, t, :it),:is]
      otpts = pn[incident(pn, t, :ot),:os]
      in_count = collect(counter(inpts))
      ot_count = collect(counter(otpts))
      Reaction(k[t], [S[i] for i in unique(inpts)],
                     [S[o] for o in unique(otpts)],
                     in_count, ot_count)
    end

    ReactionSystem(rxs, t, S, k; name=name)
end

SIR_rxn = MakeReactionSystem(SIR, :SIR)
#AlgebraicPetri.Graph(SIR)
#Catalyst.Graph(SIR_rxn)

p_real = [.1/1000, .01]
tspan = (0.0,250.0)
u0 = [999.0, 1.0, 0.0]
sample_times = range(tspan[1],stop=tspan[2],length=100)

op = ODEProblem(SIR_rxn, u0, tspan, p_real)
sol_real = solve(op, Tsit5(), tstops=sample_times)

sample_vals = [sol_real.u[findfirst(sol_real.t .>= ts)][var] * (1+(0.1rand()-0.05)) 
                for var in 1:3, ts in sample_times]

plot(sol_real, lw=2)
plot!(sample_times, sample_vals', seriestype=:scatter, label="")

function optimise_p(p_init,tend)
    function loss(p)
        sol = solve(remake(op,tspan=(0.0,tend),p=p), Tsit5(), tstops=sample_times)
        vals = hcat(map(ts -> sol.u[findfirst(sol.t .>= ts)], sample_times[1:findlast(sample_times .<= tend)])...)    
        loss = sum(abs2, vals .- sample_vals[:,1:size(vals)[2]])   
        return loss, sol
    end
    return DiffEqFlux.sciml_train(loss, p_init, maxiters=100)
end

p_estimate = [0.0,0.0]
for i in 50.0:50.0:250.0
    global p_estimate = optimise_p(p_estimate, i).minimizer
    sol_estimate = solve(remake(op,tspan=(0.0,250.0),p=p_estimate), Tsit5())
    
    plot(sample_times, sample_vals', seriestype=:scatter, label="")
    plot!(sol_estimate, lw=2)
end

#p_estimate = optimise_p([0.0,0.0], 50.0).minimizer
#=sol_estimate = solve(remake(op,tspan=(0.0,250.0),p=p_estimate), Tsit5())

plot(sample_times, sample_vals', seriestype=:scatter, label="")
plot!(sol_estimate, lw=2)=#
