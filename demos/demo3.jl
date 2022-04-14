### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 915f4df0-bc3a-11ec-2560-c9772e143679
begin 
	using Pkg, Revise
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using AlgebraicPetri
	using AlgebraicPetri: Graph
    using ModelExploration
	using Plots
	using Test
	include("../src/ModelSelection.jl")
end

# ╔═╡ 4dcde238-98a4-4982-9409-cc486a70b1e8
begin
	infectious_type_ = LabelledPetriNet([:Pop],
	  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
	  :t_infection=>(:Pop=>:Pop),
	  :t_strata=>(:Pop=>:Pop))
	
	infectious_type = map(infectious_type_, Name=name->nothing)
	
	Graph(infectious_type)
end

# ╔═╡ e2e3c56d-3350-4ca8-b073-6e6b027ca984
IO_help(i::Int) = let d = Dict([j=>j%2+1 for j in 1:2*i]); (I=d,O=d) end

# ╔═╡ 179b3378-66ff-43b9-920a-100c48b65a40
make_slice(p::LabelledPetriNet, n::NamedTuple) =
    Slice{ACSetTransformation}(
        homomorphism(p, infectious_type; initial=n,
                       type_components=(Name=x->nothing,)))

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
end

# ╔═╡ cea17847-544d-4b9d-89a2-2f2a5337df8b
begin
	diag_disease = to_slicehom(Dict(:X1=>SIR_type,:X2=>SVIIvR_type));
	diag_strata = to_slicehom(Dict(:X1=>quarantine_type, :X2=>age_s_type));
end

# ╔═╡ c5578104-3cbc-4f3d-b468-dfb702759f89
begin
	pb = PullbackSpace(to_model_hom(diag_disease), to_model_hom(diag_strata));
	upb = unfold(pb);
	aupb = apex(upb);
end;

# ╔═╡ e15eb36c-1f27-4394-98e9-77a96b75282b
begin
	true_model = ob_map(aupb, Symbol("(X2, X1)"))
	true_rxn = MakeReactionSystem(true_model)
	
	p_real = vcat(repeat([1e-4], 14), repeat([0.01], 6))
	u0 = [0.0,0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0,0.0] #|> gpu
	tspan = (0.0,250.0)
end

# ╔═╡ 1f47a2fb-8311-45c0-a178-93ddfe7cbcbf
begin
	sample_data, sample_times, prob_real, sol_real = generate_data(true_model, p_real, u0, tspan, 50)
	
	plt = plot(sol_real, lw=2, label=reshape(map(string, true_model[:, :sname]), 1, ns(true_model)))
	plot!(sample_times, sample_data, seriestype=:scatter, label="")
end

# ╔═╡ 4fa7b851-389d-4341-845b-6bf64729edf3
models = [
  ob_map(aupb, Symbol("(X1, X1)")) => [0.0, 0.0, 0.0, 999.0, 1.0, 0.0],
  ob_map(aupb, Symbol("(X1, X2)")) => [500.0, 0.0, 0.0, 499.0, 1.0, 0.0],
  ob_map(aupb, Symbol("(X2, X1)")) => [0.0, 0.0, 0.0, 0.0, 0.0, 999.0, 1.0, 0.0, 0.0, 0.0],
  ob_map(aupb, Symbol("(X2, X2)")) => [500.0, 0.0, 0.0, 0.0, 0.0, 499.0, 1.0, 0.0, 0.0, 0.0]
];

# ╔═╡ f3b0d954-2eef-4f21-9ba5-fd32a2fbdfb8
function explore(models, tspan, sample_data, sample_times)
  losses = zeros(length(models))
  sols = repeat(Any[nothing], length(models))
  println(typeof(sols))
  Threads.@threads for i in 1:length(models)
    model, u0 = models[i]
    sol, loss = full_train(model, u0, tspan, sample_data, sample_times)
    losses[i] = loss
    #println(loss)
	sols[i] = sol
  end
  return losses, sols
end

# ╔═╡ 3466eca0-f9d6-4ecf-9eeb-5b9eddc2d061
losses, sols = explore(models, tspan, sample_data, sample_times);

# ╔═╡ 98255319-3eb6-4467-ba55-4fd14cd54bef
begin
	plot(sample_times, sample_data, seriestype=:scatter, label="")
	        plot!(sols[1], lw=2, label=reshape(map(string, models[1][1][:, :sname]), 1, ns(models[1][1])))
end

# ╔═╡ 4cd3c2f4-5ecb-45f5-b6e0-cbeab132da3a
losses[1]

# ╔═╡ c56b1786-96f1-4969-ae8f-533fc0a628eb
begin
	plot(sample_times, sample_data, seriestype=:scatter, label="")
	        plot!(sols[2], lw=2, label=reshape(map(string, models[2][1][:, :sname]), 1, ns(models[2][1])))
end

# ╔═╡ 9858e6fc-60a6-417a-8092-280ecd9c139c
losses[2]

# ╔═╡ f4d7151c-35ef-4d74-a2c2-4f3e54c4461a
begin
	plot(sample_times, sample_data, seriestype=:scatter, label="")
	        plot!(sols[4], lw=2, label=reshape(map(string, models[4][1][:, :sname]), 1, ns(models[4][1])))
end

# ╔═╡ 65cf7a22-c2ca-4a44-9671-51aab7a0f644
losses[4]

# ╔═╡ e5cbb97b-a320-4835-97f1-7ff781a55998
begin
	plot(sample_times, sample_data, seriestype=:scatter, label="")
	        plot!(sols[3], lw=2, label=reshape(map(string, models[3][1][:, :sname]), 1, ns(models[3][1])))
end

# ╔═╡ 5999a061-ca85-4002-9b72-5e26bfc97939
losses[3]

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═4dcde238-98a4-4982-9409-cc486a70b1e8
# ╠═e2e3c56d-3350-4ca8-b073-6e6b027ca984
# ╠═179b3378-66ff-43b9-920a-100c48b65a40
# ╠═45511aa5-4a45-4d30-94c7-bc47c7fbfe58
# ╠═eedc21a6-5bcd-47e4-8e1b-814c232ab479
# ╠═6f2b8b20-a866-49cd-be09-f659269cbcc3
# ╠═075fa8b1-b45c-43b4-8c24-151cbd3ac7f9
# ╠═63596b64-b7aa-42ed-a980-99a684c98d4b
# ╠═cea17847-544d-4b9d-89a2-2f2a5337df8b
# ╠═c5578104-3cbc-4f3d-b468-dfb702759f89
# ╠═e15eb36c-1f27-4394-98e9-77a96b75282b
# ╠═1f47a2fb-8311-45c0-a178-93ddfe7cbcbf
# ╠═4fa7b851-389d-4341-845b-6bf64729edf3
# ╠═f3b0d954-2eef-4f21-9ba5-fd32a2fbdfb8
# ╠═3466eca0-f9d6-4ecf-9eeb-5b9eddc2d061
# ╠═98255319-3eb6-4467-ba55-4fd14cd54bef
# ╠═4cd3c2f4-5ecb-45f5-b6e0-cbeab132da3a
# ╠═c56b1786-96f1-4969-ae8f-533fc0a628eb
# ╠═9858e6fc-60a6-417a-8092-280ecd9c139c
# ╠═f4d7151c-35ef-4d74-a2c2-4f3e54c4461a
# ╠═65cf7a22-c2ca-4a44-9671-51aab7a0f644
# ╠═e5cbb97b-a320-4835-97f1-7ff781a55998
# ╠═5999a061-ca85-4002-9b72-5e26bfc97939
