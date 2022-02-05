### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ bd814b22-c126-40eb-94ff-02eeac4a2965
begin 
	using Revise, Pkg
	Pkg.develop(path="..")
	#Pkg.develop(path="../../Catlab.jl")
	#Pkg.develop(path="../../AlgebraicPetri.jl")
	using Catlab.CategoricalAlgebra
	using AlgebraicPetri
    using ModelExploration
end

# ╔═╡ 5c946c97-7484-4a58-89fa-eba5f903faac
md"""Define petri nets for epidemiology models dimension"""

# ╔═╡ 43112b41-b539-481e-9a29-fea386439544
begin
	SIR = LabelledPetriNet([:S, :I, :R],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
	
	SIR_travel_ban = LabelledPetriNet([:S, :I, :R],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :id => (:S => :S),
	  :id => (:R => :R)
	)
	
	SIS = LabelledPetriNet([:S, :I],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:S),
	  :id => (:S => :S),
	  :id => (:I => :I),
	)

	# TODO: use semagrams to define SVIIvR model
	
	epi_models = [SIS, SIR, SIR_travel_ban]
	Graph(SIR_travel_ban)
end

# ╔═╡ 8eac847e-d590-47d4-8b20-b407b620c559
md"""Define petri nets for stratification dimension"""

# ╔═╡ c34e9ad3-e43a-4d6b-981d-4dbcb82d7c0e
begin
	quarantine = LabelledPetriNet([:Q, :not_Q],
	    :id => (:Q => :Q),
	    :id => (:not_Q => :not_Q),
	    :enter_quarantine => (:not_Q => :Q),
	    :exit_quarantine => (:Q => :not_Q),
	    :interaction => ((:not_Q, :not_Q) => (:not_Q, :not_Q))
	)
	
	age_stratification = LabelledPetriNet([:Child, :Adult],
	    :id => (:Child => :Child),
	    :id => (:Adult => :Adult),
	    :interaction => ((:Child, :Child) => (:Child, :Child)),
	    :interaction => ((:Adult, :Adult) => (:Adult, :Adult)),
	    :interaction => ((:Child, :Adult) => (:Child, :Adult))
	)
	
	flux_metapopulation = LabelledPetriNet([:Patch1, :Patch2],
	    :travel => (:Patch1 => :Patch2),
	    :travel => (:Patch2 => :Patch1),
	    :id => (:Patch1 => :Patch1),
	    :id => (:Patch2 => :Patch2),
	    :interaction => ((:Patch1, :Patch1) => (:Patch1, :Patch1)),
	    :interaction => ((:Patch2, :Patch2) => (:Patch2, :Patch2))
	)
	
	# TODO: use semagrams to make a simple trip stratification
	strat_models = [quarantine, age_stratification, flux_metapopulation]
	Graph(quarantine)
end

# ╔═╡ 384f2536-e54d-48e8-bd6e-2c58ce81aa26
md"""Declare "literal" generators"""

# ╔═╡ 70d70a1e-63ba-495a-a7e0-33d677ea56cf
begin
	dim1 = Literal(epi_models)
	dim2 = Literal(strat_models)
end

# ╔═╡ dc439087-d393-416f-b0b4-57a3c2c8d69c
md"""Construct the model space as Product([dim1, dim2], slice)"""

# ╔═╡ 36849bbb-b1e0-4fca-b633-fe04ac3377ac
begin
	infectious_type = LabelledPetriNet([:Pop],
	  :interaction=>((:Pop, :Pop)=>(:Pop, :Pop)),
	  :t_infection=>(:Pop=>:Pop),
	  :t_strata=>(:Pop=>:Pop)
	)
	
	prodSpace= ModelExploration.Product(dim1, dim2, infectious_type)
	Graph(infectious_type)
end

# ╔═╡ 21f0b44f-6868-4b2c-b235-35ea91fd3046
md"""Define the loss function
We've already computed the simulated data from the desired Petri net (load from file)
"""

# ╔═╡ 33ae7c95-ded4-45d6-a600-4d8b409b16a0
md"""Run selection, should return the Petri net that generated the data"""

# ╔═╡ 2ec8c472-d3de-4fc9-875b-c5815b7d22f7
ModelExploration.select(prodSpace, ()->1.)

# ╔═╡ Cell order:
# ╠═bd814b22-c126-40eb-94ff-02eeac4a2965
# ╟─5c946c97-7484-4a58-89fa-eba5f903faac
# ╠═43112b41-b539-481e-9a29-fea386439544
# ╟─8eac847e-d590-47d4-8b20-b407b620c559
# ╠═c34e9ad3-e43a-4d6b-981d-4dbcb82d7c0e
# ╟─384f2536-e54d-48e8-bd6e-2c58ce81aa26
# ╠═70d70a1e-63ba-495a-a7e0-33d677ea56cf
# ╟─dc439087-d393-416f-b0b4-57a3c2c8d69c
# ╠═36849bbb-b1e0-4fca-b633-fe04ac3377ac
# ╟─21f0b44f-6868-4b2c-b235-35ea91fd3046
# ╟─33ae7c95-ded4-45d6-a600-4d8b409b16a0
# ╠═2ec8c472-d3de-4fc9-875b-c5815b7d22f7
