### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ bd814b22-c126-40eb-94ff-02eeac4a2965
begin 
	using Revise, Pkg
	Pkg.develop(path="..")
	Pkg.develop(path="../../Catlab.jl")
	Pkg.develop(path="../../AlgebraicPetri.jl")
	using Catlab.CategoricalAlgebra, Catlab.Graphs
    using ModelExploration
end

# ╔═╡ b63b47a0-85fd-11ec-15ca-e3b0dacf37f1
md"""Hello $\frac{1}{1}$"""

# ╔═╡ 5c946c97-7484-4a58-89fa-eba5f903faac
md"""Define petri nets for two dimensions"""

# ╔═╡ 384f2536-e54d-48e8-bd6e-2c58ce81aa26
md"""Declare "literal" generators"""

# ╔═╡ 70d70a1e-63ba-495a-a7e0-33d677ea56cf


# ╔═╡ dc439087-d393-416f-b0b4-57a3c2c8d69c
md"""Construct the model space as Product([dim1, dim2], slice)"""

# ╔═╡ ffeed7a4-1414-4ee5-85d6-ed74c0dab7c5
dim1  = dim2 = Literal(StructACSet[])


# ╔═╡ 36849bbb-b1e0-4fca-b633-fe04ac3377ac
prodSpace= ModelExploration.Product(dim1,dim2, Graph())

# ╔═╡ 21f0b44f-6868-4b2c-b235-35ea91fd3046
md"""Define the loss function
We've already computed the simulated data from the desired Petri net (load from file)
"""

# ╔═╡ 33ae7c95-ded4-45d6-a600-4d8b409b16a0
md"""Run selection, should return the Petri net that generated the data"""

# ╔═╡ 2ec8c472-d3de-4fc9-875b-c5815b7d22f7
ModelExploration.select(prodSpace, ()->1.)

# ╔═╡ Cell order:
# ╟─b63b47a0-85fd-11ec-15ca-e3b0dacf37f1
# ╠═bd814b22-c126-40eb-94ff-02eeac4a2965
# ╟─5c946c97-7484-4a58-89fa-eba5f903faac
# ╠═384f2536-e54d-48e8-bd6e-2c58ce81aa26
# ╠═70d70a1e-63ba-495a-a7e0-33d677ea56cf
# ╟─dc439087-d393-416f-b0b4-57a3c2c8d69c
# ╠═ffeed7a4-1414-4ee5-85d6-ed74c0dab7c5
# ╠═36849bbb-b1e0-4fca-b633-fe04ac3377ac
# ╟─21f0b44f-6868-4b2c-b235-35ea91fd3046
# ╟─33ae7c95-ded4-45d6-a600-4d8b409b16a0
# ╠═2ec8c472-d3de-4fc9-875b-c5815b7d22f7
