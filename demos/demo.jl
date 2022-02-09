### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ bd814b22-c126-40eb-94ff-02eeac4a2965
begin 
	using Pkg
	Pkg.activate(Base.current_project())
    # Pkg.add([
    #     Pkg.PackageSpec(name="Catalyst", version="10.1"),
    # ])
    # using Catalyst

	Pkg.instantiate()
	#Pkg.develop(path="..")
	#Pkg.add(url="https://github.com/kris-brown/Catlab.jl", rev="looseacset_csp")
	#Pkg.add(url="https://github.com/AlgebraicJulia/AlgebraicPetri.jl", rev="Catalyst_v10")
	using Revise
	using Catlab.CategoricalAlgebra
	using AlgebraicPetri
	using AlgebraicPetri: Graph
    using ModelExploration
	using Plots
	using Test
end

# ╔═╡ 28a03065-8d42-402d-b140-13b7f87001b6
md"""Define typing system for models"""

# ╔═╡ 36849bbb-b1e0-4fca-b633-fe04ac3377ac
begin
	infectious_type = LabelledPetriNet([:Pop],
	  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
	  :t_infection=>(:Pop=>:Pop),
	  :t_strata=>(:Pop=>:Pop))	
	Graph(infectious_type)
end

# ╔═╡ 5c946c97-7484-4a58-89fa-eba5f903faac
md"""Define petri nets for epidemiology models dimension"""

# ╔═╡ 43112b41-b539-481e-9a29-fea386439544
begin
	IO_help(i::Int) = let d = Dict([j=>j%2+1 for j in 1:2*i]); (I=d,O=d) end
	
	SIR = LabelledPetriNet([:S, :I, :R],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
	SIR_type = merge((T=[1,2,3,3,3],), IO_help(1))
	
	SIR_travel_ban = LabelledPetriNet([:S, :I, :R],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :id => (:S => :S),
	  :id => (:R => :R)
	)

	SIR_tb_type = merge((T=[1,2,3,3],), IO_help(1))
	
	SIS = LabelledPetriNet([:S, :I],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:S),
	  :id => (:S => :S),
	  :id => (:I => :I),
	)
	
	SIS_type = merge((T=[1,2,3,3],),IO_help(1))

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
	
	SVIIvR_type = merge((T=[1,1,1,1,2,2,2,3,3,3,3,3],),IO_help(4))
	
end;

# ╔═╡ 148087ee-ed18-4e1e-aecb-a3736614c4de
	epi_models = SliceLit(infectious_type, [
		SIR=>SIR_type, SIS=>SIS_type, 
		SIR_travel_ban=>SIR_tb_type, SVIIvR=>SVIIvR_type]);


# ╔═╡ 8eac847e-d590-47d4-8b20-b407b620c559
md"""Define petri nets for stratification dimension"""

# ╔═╡ c34e9ad3-e43a-4d6b-981d-4dbcb82d7c0e
begin
	quarantine = LabelledPetriNet([:Q, :not_Q],
	    :interaction => ((:not_Q, :not_Q) => (:not_Q, :not_Q)),
	    :enter_quarantine => (:not_Q => :Q),
	    :exit_quarantine => (:Q => :not_Q),
	    :id => (:Q => :Q),
	    :id => (:not_Q => :not_Q),
	)
	quarantine_type = (T=[1, 3, 3, 2, 2], I=Dict(1=>1,2=>2), O=Dict(1=>1,2=>2))
	
	age_stratification = LabelledPetriNet([:Child, :Adult],
	    :interaction => ((:Child, :Child) => (:Child, :Child)),
	    :interaction => ((:Adult, :Adult) => (:Adult, :Adult)),
	    :interaction => ((:Child, :Adult) => (:Child, :Adult)),
	    :id => (:Child => :Child),
	    :id => (:Adult => :Adult),
	)
	age_s_type = merge((T=[1,1,1,2,2],),IO_help(3))
	
	flux_metapopulation = LabelledPetriNet([:Patch1, :Patch2],
	    :interaction => ((:Patch1, :Patch1) => (:Patch1, :Patch1)),
	    :interaction => ((:Patch2, :Patch2) => (:Patch2, :Patch2)),	    
		:travel => (:Patch1 => :Patch2),
	    :travel => (:Patch2 => :Patch1),
	    :id => (:Patch1 => :Patch1),
	    :id => (:Patch2 => :Patch2),
	)
	flux_m_type = merge((T=[1,1,3,3,2,2],), IO_help(2))

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
	
	simple_t_type = merge((T=[1,1,1,1,1,1,1,1,3,3,3,3,2,2,2,2],), IO_help(8))
end;

# ╔═╡ 35a823ec-260b-4a15-876b-30d3b6a090a6
	strat_models = SliceLit(infectious_type, [quarantine => quarantine_type, age_stratification => age_s_type, flux_metapopulation => flux_m_type, simple_trip => simple_t_type ]);


# ╔═╡ dc439087-d393-416f-b0b4-57a3c2c8d69c
md"""Construct the model space as Product([dim1, dim2], slice)"""

# ╔═╡ 634c732d-6c68-4a2d-a570-11e149358973
prodSpace= ModelExploration.Product(epi_models, strat_models, infectious_type);

# ╔═╡ b6a6d27d-cae4-4fad-9f59-e1feccad5cea
begin 
	true_model = unfold(prodSpace)[2]
	Graph(true_model)
end

# ╔═╡ f3092427-d6b8-4241-a319-b17e968ca735
true_model[:name]

# ╔═╡ 0af79678-19bd-4ddd-a5e7-08843dcbd289
ns(true_model)

# ╔═╡ 7664c213-a569-4069-9b0e-45db9d35344f
begin 
	true_rates = Float64[1e-3,1e-5,1e-4,1e-7,1e-4]
	true_initial_pop = Float64[1e6,1e6,1e3,1e2,1e-2,1e-2]
	loss_fun = eval_petri_fn(true_model, true_rates,true_initial_pop)
end

# ╔═╡ 7fb8a18a-d0f3-40a2-bc78-e40651701d65
plot(generate_data(true_model, true_rates,true_initial_pop))


# ╔═╡ 21f0b44f-6868-4b2c-b235-35ea91fd3046
md"""Define the loss function
We've already computed the simulated data from the desired Petri net (load from file)
"""

# ╔═╡ 33ae7c95-ded4-45d6-a600-4d8b409b16a0
md"""Run selection, should return the Petri net that generated the data"""

# ╔═╡ 2ec8c472-d3de-4fc9-875b-c5815b7d22f7
begin 
	best_model = Graph(ModelExploration.select(prodSpace, loss_fun))
    @test best_model == true_model 
end

# ╔═╡ Cell order:
# ╠═bd814b22-c126-40eb-94ff-02eeac4a2965
# ╟─28a03065-8d42-402d-b140-13b7f87001b6
# ╠═36849bbb-b1e0-4fca-b633-fe04ac3377ac
# ╟─5c946c97-7484-4a58-89fa-eba5f903faac
# ╠═43112b41-b539-481e-9a29-fea386439544
# ╠═148087ee-ed18-4e1e-aecb-a3736614c4de
# ╟─8eac847e-d590-47d4-8b20-b407b620c559
# ╟─c34e9ad3-e43a-4d6b-981d-4dbcb82d7c0e
# ╠═35a823ec-260b-4a15-876b-30d3b6a090a6
# ╠═dc439087-d393-416f-b0b4-57a3c2c8d69c
# ╠═634c732d-6c68-4a2d-a570-11e149358973
# ╠═b6a6d27d-cae4-4fad-9f59-e1feccad5cea
# ╠═f3092427-d6b8-4241-a319-b17e968ca735
# ╠═0af79678-19bd-4ddd-a5e7-08843dcbd289
# ╠═7664c213-a569-4069-9b0e-45db9d35344f
# ╠═7fb8a18a-d0f3-40a2-bc78-e40651701d65
# ╟─21f0b44f-6868-4b2c-b235-35ea91fd3046
# ╟─33ae7c95-ded4-45d6-a600-4d8b409b16a0
# ╠═2ec8c472-d3de-4fc9-875b-c5815b7d22f7
