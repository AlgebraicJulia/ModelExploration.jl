using AlgebraicPetri
using Catlab.CategoricalAlgebra
using Catalyst
using OrdinaryDiffEq
using DiffEqFlux #=, Flux=#
using Plots
#import Catalyst: ReactionSystem

println("Starting")

#=SIR = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
  :id => (:S => :S),
  :id => (:I => :I),
  :id => (:R => :R)
)=#

get_infected_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("I", string(s))]

get_susceptible_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("S", string(s))]

get_recovered_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("R", string(s)) || occursin("V", string(s))]

counter(a) = [count(==(i),a) for i in unique(a)]
function MakeReactionSystem(pn::AbstractPetriNet)
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

    ReactionSystem(rxs, t, S, k; name=:model)
end

#=SIR_rxn = MakeReactionSystem(SIR)
#AlgebraicPetri.Graph(SIR)
#Catalyst.Graph(SIR_rxn)

p_real = [.1/1000, .01]
tspan = (0.0,250.0)
u0 = [999.0, 1.0, 0.0]=#
#=sample_times = range(tspan[1],stop=tspan[2],length=100)

op = ODEProblem(SIR_rxn, u0, tspan, p_real)
sol_real = solve(op, Tsit5(), tstops=sample_times)

sample_vals = [sol_real.u[findfirst(sol_real.t .>= ts)][var] * (1+(0.1rand()-0.05)) 
                for var in 1:3, ts in sample_times]

plt = plot(sol_real, lw=2, label=reshape(map(string, SIR[:, :sname]), 1, ns(SIR)))
plot!(sample_times, sample_vals', seriestype=:scatter, label="")
display(plt)=#

"""
Generate infection data from a petri net model given:
  - rate parameters p
  - initial concentrations u0
  - timespan
  - number of samples
Returns sample infectious data, problem, and solution
"""
function generate_data(model::AbstractLabelledPetriNet, p, u0, tspan, num_samples)
    sample_times = range(tspan[1], stop=tspan[2], length=num_samples)

    prob = ODEProblem(MakeReactionSystem(model), u0, tspan, p)
    sol = solve(prob, Tsit5(), tstops=sample_times)

    inf_sample_vals = [sol.u[findfirst(sol.t .>= ts)][var] * (1+(0.1rand()-0.05))
                for var in get_infected_states(model), ts in sample_times]

    total_inf_samples = map(sum, eachcol(inf_sample_vals))

    susc_sample_vals = [sol.u[findfirst(sol.t .>= ts)][var] * (1+(0.1rand()-0.05))
                for var in get_susceptible_states(model), ts in sample_times]

    total_susc_samples = map(sum, eachcol(susc_sample_vals))

    rec_sample_vals = [sol.u[findfirst(sol.t .>= ts)][var] * (1+(0.1rand()-0.05))
                for var in get_recovered_states(model), ts in sample_times]

    total_rec_samples = map(sum, eachcol(rec_sample_vals))

    return hcat(total_inf_samples, total_rec_samples, total_susc_samples), sample_times, prob, sol
end


function optimise_p(model::AbstractLabelledPetriNet, op, p_init, tend, sample_data, sample_times)
    # Loss on all populations
    #=loss = function (p)
        sol = solve(remake(op,tspan=(0.0,tend),p=p), Tsit5(), tstops=sample_times)
        vals = hcat(map(ts -> sol.u[findfirst(sol.t .>= ts)], sample_times[1:findlast(sample_times .<= tend)])...)    
        loss = sum(abs2, vals .- sample_vals[:,1:size(vals)[2]])   
        return loss, sol
    end=#

    # Loss on just infected population (requires closer initial values)
    loss = function (p)
        sol = solve(remake(op,tspan=(0.0,tend),p=p), Tsit5(), tstops=sample_times)
        vals = hcat(map(ts -> sol.u[findfirst(sol.t .>= ts)], sample_times[1:findlast(sample_times .<= tend)])...)
        inf_vals = map(sum, collect(zip([vals[i,:] for i in get_infected_states(model)]...)))

        #loss = sum(abs2, vals[2,:] .- sample_vals[:,1:size(vals)[2]][2,:])
        inf_samples = sample_data[:,1]
        rec_samples = sample_data[:,2]
        susc_samples = sample_data[:,3]

        inf_loss = sum(abs2, inf_vals .- inf_samples[1:size(inf_vals)[1]])

        susc_vals = map(sum, collect(zip([vals[i,:] for i in get_susceptible_states(model)]...)))
        susc_loss = sum(abs2, susc_vals .- susc_samples[1:size(susc_vals)[1]])

        rec_vals = map(sum, collect(zip([vals[i,:] for i in get_recovered_states(model)]...)))
        rec_loss = sum(abs2, susc_vals .- rec_samples[1:size(rec_vals)[1]])
        return susc_loss + inf_loss, sol
    end

    callback = function (p, l, pred)
        display(l)
        plt = plot(pred, lw=2, linestyle=:dash, color=[:orange :green :blue])
        plot!(sample_times, sample_vals', seriestype=:scatter, label="")
        display(plt)
        return false
    end

    return DiffEqFlux.sciml_train(loss, p_init, #=ADAM(0.05), cb=callback,=# maxiters=500, abstol=1e-4, reltol=1e-4,
            lower_bounds=repeat([1e-6], length(p_init)), upper_bounds=ones(length(p_init)))
end

function full_train(model, u0, tspan, training_data, sample_times, param_guess)
    rxn = MakeReactionSystem(model)
    #p_estimate = zeros(numreactionparams(rxn))
    #p_estimate = repeat([1e-6], numreactionparams(rxn)-ns(model))
    p_estimate = param_guess
    prob = ODEProblem(rxn, u0, tspan, p_estimate)
    loss = 0
    sol_estimate = Nothing
    # TODO: algorithmically determine this range
    for i in 50:50:250
        res_ode = optimise_p(model, prob, p_estimate, i, training_data, sample_times)
        p_estimate = res_ode.minimizer
        sol_estimate = solve(remake(prob,tspan=tspan,p=p_estimate), Tsit5())
        
        susc_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_susceptible_states(model)]...)))
        rec_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_recovered_states(model)]...)))
        inf_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_infected_states(model)]...)))

        #=plt = plot(sample_times, training_data, seriestype=:scatter, label="")
        plot!(susc_vals, lw=2, label="S")
        plot!(rec_vals, lw=2, label="R")
        plot!(inf_vals, lw=2, label="I")
        display(plt)=#
        plt = plot(sample_times, training_data, seriestype=:scatter, label="")
        plot!(sol_estimate, lw=2, label=reshape(map(string, model[:, :sname]), 1, ns(model)))
        display(plt)

        println("Estimated params: $p_estimate")
        loss = res_ode.minimum
        println("Loss: $loss")
    end
    return  sol_estimate, loss
end

function full_train(model, u0, tspan, training_data, sample_times)
    rxn = MakeReactionSystem(model)
    #p_estimate = zeros(numreactionparams(rxn))
    p_estimate = repeat([1e-5], numreactionparams(rxn)-ns(model))
    #p_estimate = param_guess
    prob = ODEProblem(rxn, u0, tspan, p_estimate)
    loss = 0
    sol_estimate = Nothing
    # TODO: algorithmically determine this range
    for i in 50:50:250
        res_ode = optimise_p(model, prob, p_estimate, i, training_data, sample_times)
        p_estimate = res_ode.minimizer
        sol_estimate = solve(remake(prob,tspan=tspan,p=p_estimate), Tsit5())
        
        susc_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_susceptible_states(model)]...)))
        rec_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_recovered_states(model)]...)))
        inf_vals = map(sum, collect(zip([sol_estimate[i,:] for i in get_infected_states(model)]...)))

        #=plt = plot(sample_times, training_data, seriestype=:scatter, label="")
        plot!(susc_vals, lw=2, label="S")
        plot!(rec_vals, lw=2, label="R")
        plot!(inf_vals, lw=2, label="I")
        display(plt)=#
        plt = plot(sample_times, training_data, seriestype=:scatter, label="")
        plot!(sol_estimate, lw=2, label=reshape(map(string, model[:, :sname]), 1, ns(model)))
        display(plt)

        println("Estimated params: $p_estimate")
        loss = res_ode.minimum
        println("Loss: $loss")
    end
    return  sol_estimate, loss
end

function run_param_est()
    println("Generating data...")

    data, sample_times, op, sol_real = generate_data(SIR, p_real, u0, tspan, 100)
    plt = plot(sol_real, lw=2, label=reshape(map(string, SIR[:, :sname]), 1, ns(SIR)))
    plot!(sample_times, data, seriestype=:scatter, label="")
    display(plt)

    println("Running param estimation...")
    p_estimate = [0.0,0.0]
    for i in 50.0:50.0:250.0
        res_ode = optimise_p(SIR, op, p_estimate, i, data, sample_times)
        p_estimate = res_ode.minimizer
        sol_estimate = solve(remake(op,tspan=(0.0,250.0),p=p_estimate), Tsit5())
    
        plt = plot(sample_times, sample_vals', seriestype=:scatter, label="")
        plot!(sol_estimate, lw=2, label=reshape(map(string, SIR[:, :sname]), 1, ns(SIR)))
        display(plt)

        println("Estimated params: $p_estimate")
        loss = res_ode.minimum
        println("Loss: $loss")
    end

    sol_estimate = solve(remake(op,tspan=(0.0,250.0),p=p_estimate), Tsit5())
    return sol_estimate, p_estimate
end

#sol_estimate, p_estimate = run_param_est()

#=plot(sample_times, sample_vals', seriestype=:scatter, label="")
plot!(sol_estimate, lw=2)=#
#println("Final param estimate: $p_estimate")
