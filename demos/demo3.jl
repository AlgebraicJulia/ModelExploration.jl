### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 915f4df0-bc3a-11ec-2560-c9772e143679
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using AlgebraicPetri
	using AlgebraicPetri: Graph
    using ModelExploration
	using Plots
	using Test
	using Logging
	include("../src/ModelSelection.jl")
	Logging.disable_logging(Logging.Warn)
end;

# ╔═╡ 6b70e364-152e-49c4-ae1b-381cbbcb6fc3
md"""Define typing system for models"""

# ╔═╡ 4dcde238-98a4-4982-9409-cc486a70b1e8
begin
	infectious_type_ = LabelledPetriNet([:Pop],
	  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
	  :t_infection=>(:Pop=>:Pop),
	  :t_strata=>(:Pop=>:Pop))
	
	infectious_type = map(infectious_type_, Name=name->nothing)
	
	Graph(infectious_type_)
end

# ╔═╡ e2e3c56d-3350-4ca8-b073-6e6b027ca984
IO_help(i::Int) = let d = Dict([j=>j%2+1 for j in 1:2*i]); (I=d,O=d) end;

# ╔═╡ 179b3378-66ff-43b9-920a-100c48b65a40
make_slice(p::LabelledPetriNet, n::NamedTuple) =
    Slice{ACSetTransformation}(
        homomorphism(p, infectious_type; initial=n,
                       type_components=(Name=x->nothing,)))

# ╔═╡ 712853b5-a679-4c5e-a47b-34af8ee8ea0b
md"""Define petri nets for disease dynamics dimension"""

# ╔═╡ 45511aa5-4a45-4d30-94c7-bc47c7fbfe58
begin
	SIR = LabelledPetriNet([:S, :I, :R],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
	SIR_type = make_slice(SIR, merge((T=[1,2,3,3,3],), IO_help(1)))
	Graph(SIR)
end

# ╔═╡ eedc21a6-5bcd-47e4-8e1b-814c232ab479
begin
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
	Graph(SVIIvR)
end

# ╔═╡ 78d94905-418c-4cfc-91c5-9dac9bcd4049
md"""Define petri nets for stratification dimension"""

# ╔═╡ 6f2b8b20-a866-49cd-be09-f659269cbcc3
begin
	quarantine = LabelledPetriNet([:Q, :not_Q],
	    :interaction => ((:not_Q, :not_Q) => (:not_Q, :not_Q)),
	    :enter_quarantine => (:not_Q => :Q),
	    :exit_quarantine => (:Q => :not_Q),
	    :id => (:Q => :Q),
	    :id => (:not_Q => :not_Q),
	)
	quarantine_type = make_slice(quarantine, (T=[1, 3, 3, 2, 2], I=Dict(1=>1,2=>2), O=Dict(1=>1,2=>2)))
	Graph(quarantine)
end

# ╔═╡ 075fa8b1-b45c-43b4-8c24-151cbd3ac7f9
begin
	age_stratification = LabelledPetriNet([:Child, :Adult],
	    :interaction => ((:Child, :Child) => (:Child, :Child)),
	    :interaction => ((:Adult, :Adult) => (:Adult, :Adult)),
	    :interaction => ((:Child, :Adult) => (:Child, :Adult)),
	    :id => (:Child => :Child),
	    :id => (:Adult => :Adult),
	)
	age_s_type = make_slice(age_stratification, merge((T=[1,1,1,2,2],),IO_help(3)))
	Graph(age_stratification)
end

# ╔═╡ 63596b64-b7aa-42ed-a980-99a684c98d4b
begin
	@present ThTwo(FreeCategory) begin
	    (X1,X2)::Ob
	end
	@present ThOne(FreeCategory) begin
	  X::Ob
	end
	
	const ACSetCat{S} = TypeCat{S, ACSetTransformation}
	const ACSetCatSlice{S} = TypeCat{S, ACSetTransformation}
	const Petri = ACSetCat{LabelledPetriNet}
	const PetriHom = SliceCat{ACSetTransformation}
	
	One, Two = FinCat.([ThOne,ThTwo])
	
	to_diag(x) = Diagram(FinDomFunctor(x, nothing, Two, PetriHom()))
	to_slicehom(x) = SliceDiagHom(Literal(to_diag(x)))
end;

# ╔═╡ cc7a7deb-af14-4d15-a444-7a021d765168
md"""Construct the model space as the product of two dimensions"""

# ╔═╡ cea17847-544d-4b9d-89a2-2f2a5337df8b
begin
	diag_disease = to_slicehom(Dict(:X1=>SIR_type, :X2=>SVIIvR_type));
	diag_strata = to_slicehom(Dict(:X1=>quarantine_type, :X2=>age_s_type));
end;

# ╔═╡ c5578104-3cbc-4f3d-b468-dfb702759f89
begin
	# This object contains the recipe to build the composite model space
	pb = PullbackSpace(to_model_hom(diag_disease), to_model_hom(diag_strata));
	# This function evaluates the recipe
	upb = unfold(pb);
	# This discard the projection maps, just keeps the product models
	aupb = apex(upb);
end;

# ╔═╡ 5a4b68ce-0546-43ac-911a-00e638b0b0df
md"""For example, combining SIR with quarantine gives us two copies of SIR (one quarantined, one not) with enter/exit quarantine transitions for each"""

# ╔═╡ 131476ed-ed71-4d5f-9fc0-708f1490c19d
Graph(ob_map(aupb, Symbol("(X1, X1)")))

# ╔═╡ 3dcc4486-0961-44ee-ad90-0c286972fa12
md"""Simulate a trajectory based on the combination of the vaccinated model + quarantine. Arbitrary but realistic parameters, 250 time steps."""

# ╔═╡ e15eb36c-1f27-4394-98e9-77a96b75282b
begin
	true_model = ob_map(aupb, Symbol("(X2, X1)"))
	true_rxn = MakeReactionSystem(true_model)
	
	p_real = vcat(repeat([1e-4], 14), repeat([0.01], 6))
	p_real = vcat([1.1e-4,1.3e-4,1.2e-4,1.4e-4,1.5e-4,1.2e-4,
				   1.e-4,1.3e-4,1.5e-4,1.2e-4,1.7e-4,1.6e-4,
		           1.2e-4,1.1e-4,], [0.01,0.02,0.01,0.01,0.01,0.02])
	u0 = [0.0,0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0,0.0] 
	tspan = (0.0,250.0)
end;

# ╔═╡ 1f47a2fb-8311-45c0-a178-93ddfe7cbcbf
begin
	sample_data, sample_times, prob_real, sol_real = generate_data(true_model, p_real, u0, tspan, 50)
	
	plt = plot(sol_real, lw=2, label=reshape(map(string, true_model[:, :sname]), 1, ns(true_model)))
	plot!(sample_times, sample_data, seriestype=:scatter, label="")
end

# ╔═╡ b44ff4bc-92c3-4312-aad3-c1d0878cdf46
md"""Equip models with initial guesses ![](https://i.imgur.com/VKwUP9w.png)
"""

# ╔═╡ 8a0d7fb4-49f5-47be-a890-ea7900ddc58d
begin 
	
disease_initial = (S=999., I=1., R=0.,V=0.,Iv=0.) 
strata_initial = (Q=0., not_Q=1.,Child=0.5,Adult=0.5)
	
models = map(ob_generators(dom(diagram(aupb)))) do o
	model = ob_map(aupb, o)
	init_pop =  [disease_initial[n1] * strata_initial[n2] 
				 for (n1, n2) in model[:sname]]
	model => init_pop
end
end;

# ╔═╡ f3b0d954-2eef-4f21-9ba5-fd32a2fbdfb8
function explore(models, tspan, sample_data, sample_times)
  losses = zeros(length(models))
  sols = repeat(Any[nothing], length(models))
  Threads.@threads for i in 1:length(models)
    model, u0 = models[i]
    sol, loss = full_train(model, u0, tspan, sample_data, sample_times)
    losses[i] = loss
	sols[i] = sol
  end
  return losses, sols
end;

# ╔═╡ 5e12700b-d7d9-44c3-9728-71d3c08b18c1
md"""Optimizing all models can take minutes"""

# ╔═╡ 3466eca0-f9d6-4ecf-9eeb-5b9eddc2d061
losses, sols = explore(models, tspan, sample_data, sample_times);

# ╔═╡ 6e81b7d5-f225-42ee-a60d-9076c274a33b
function plot_res(i::Int)
	plot(sample_times, sample_data, seriestype=:scatter, label="")
	        plot!(sols[i], lw=2, label=reshape(map(string, models[i][1][:, :sname]), 1, ns(models[i][1])))
end;

# ╔═╡ da5ab695-e5c3-49b7-b849-afb03a81c598
md"""Model (1,1)"""

# ╔═╡ 4cd3c2f4-5ecb-45f5-b6e0-cbeab132da3a
losses[1]

# ╔═╡ 90ee4cd1-2957-49cf-8f29-fb0c26db3219
plot_res(1)

# ╔═╡ 5d3ae1b2-969a-407c-af4f-452212e4c6f8
md"""Model (1,2)"""

# ╔═╡ 9858e6fc-60a6-417a-8092-280ecd9c139c
losses[3]

# ╔═╡ a9d58251-e0f3-4a77-8b0e-afd4f24c31c4
plot_res(3)

# ╔═╡ b836c19d-b5aa-45d1-8487-586042cbcda2
md"""Model (2,2)"""

# ╔═╡ 65cf7a22-c2ca-4a44-9671-51aab7a0f644
losses[4]

# ╔═╡ 8e5b3b93-d023-4e79-b174-1356b560fd20
plot_res(4)

# ╔═╡ a47b8757-8723-4dfe-a81c-6adf6981f356
md"""Model (2,1): The correct model"""

# ╔═╡ 5999a061-ca85-4002-9b72-5e26bfc97939
losses[2]

# ╔═╡ 3322866c-61bd-40c0-9a59-50610fe37a0d
plot_res(2)

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╟─6b70e364-152e-49c4-ae1b-381cbbcb6fc3
# ╠═4dcde238-98a4-4982-9409-cc486a70b1e8
# ╟─e2e3c56d-3350-4ca8-b073-6e6b027ca984
# ╟─179b3378-66ff-43b9-920a-100c48b65a40
# ╟─712853b5-a679-4c5e-a47b-34af8ee8ea0b
# ╠═45511aa5-4a45-4d30-94c7-bc47c7fbfe58
# ╠═eedc21a6-5bcd-47e4-8e1b-814c232ab479
# ╟─78d94905-418c-4cfc-91c5-9dac9bcd4049
# ╠═6f2b8b20-a866-49cd-be09-f659269cbcc3
# ╠═075fa8b1-b45c-43b4-8c24-151cbd3ac7f9
# ╟─63596b64-b7aa-42ed-a980-99a684c98d4b
# ╟─cc7a7deb-af14-4d15-a444-7a021d765168
# ╠═cea17847-544d-4b9d-89a2-2f2a5337df8b
# ╠═c5578104-3cbc-4f3d-b468-dfb702759f89
# ╟─5a4b68ce-0546-43ac-911a-00e638b0b0df
# ╠═131476ed-ed71-4d5f-9fc0-708f1490c19d
# ╟─3dcc4486-0961-44ee-ad90-0c286972fa12
# ╠═e15eb36c-1f27-4394-98e9-77a96b75282b
# ╠═1f47a2fb-8311-45c0-a178-93ddfe7cbcbf
# ╟─b44ff4bc-92c3-4312-aad3-c1d0878cdf46
# ╠═8a0d7fb4-49f5-47be-a890-ea7900ddc58d
# ╟─f3b0d954-2eef-4f21-9ba5-fd32a2fbdfb8
# ╟─5e12700b-d7d9-44c3-9728-71d3c08b18c1
# ╠═3466eca0-f9d6-4ecf-9eeb-5b9eddc2d061
# ╠═6e81b7d5-f225-42ee-a60d-9076c274a33b
# ╟─da5ab695-e5c3-49b7-b849-afb03a81c598
# ╟─4cd3c2f4-5ecb-45f5-b6e0-cbeab132da3a
# ╠═90ee4cd1-2957-49cf-8f29-fb0c26db3219
# ╟─5d3ae1b2-969a-407c-af4f-452212e4c6f8
# ╟─9858e6fc-60a6-417a-8092-280ecd9c139c
# ╠═a9d58251-e0f3-4a77-8b0e-afd4f24c31c4
# ╟─b836c19d-b5aa-45d1-8487-586042cbcda2
# ╟─65cf7a22-c2ca-4a44-9671-51aab7a0f644
# ╠═8e5b3b93-d023-4e79-b174-1356b560fd20
# ╟─a47b8757-8723-4dfe-a81c-6adf6981f356
# ╟─5999a061-ca85-4002-9b72-5e26bfc97939
# ╠═3322866c-61bd-40c0-9a59-50610fe37a0d
