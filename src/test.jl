using ModelingToolkit, Catalyst #, LinearAlgebra
#using DiffEqBase, DiffEqJump
#using Plots, SpecialFunctions

## Parameter
N = 10                       # maximum cluster size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / bulk volume
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³

integ(x) = Int(floor(x))
n        = integ(N/2)
nr       = N%2 == 0 ? (n*(n + 1) - n) : (n*(n + 1)) # No. of forward reactions

# possible pairs of reactant multimers
pair = []
for i = 2:N
    push!(pair,[1:integ(i/2)  i .- (1:integ(i/2))])
end
pair = vcat(pair...)
vᵢ = @view pair[:,1]   # Reactant 1 indices
vⱼ = @view pair[:,2]   # Reactant 2 indices
volᵢ = Vₒ*vᵢ           # cm⁻³
volⱼ = Vₒ*vⱼ           # cm⁻³
sum_vᵢvⱼ = @. vᵢ + vⱼ  # Product index

# set i to  1 for additive kernel, 2  for constant
i = 1
if i==1
    B = 1.53e03                # s⁻¹
    kv = @. B*(volᵢ + volⱼ)/V  # dividing by volume as its a bi-molecular reaction chain
elseif i==2
    C = 1.84e-04               # cm³ s⁻¹
    kv = fill(C/V, nr) 
end

# state variables are X, pars stores rate parameters for each rx
@parameters t       
@variables k[1:nr]  X[1:N](t)
pars = Pair.(collect(k), kv)

# time-span
if i == 1
    tspan = (0. ,2000.)   
elseif i == 2
    tspan = (0. ,350.)
end

 # initial condition of monomers
u₀    = zeros(Int64, N)
u₀[1] = uₒ  
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value

# vector to store the Reactions in
rx = []              
for n = 1:nr
    # for clusters of the same size, double the rate
    if (vᵢ[n] == vⱼ[n]) 
        push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
    else
        push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]], 
                           [1, 1], [1]))
    end
end
@named rs = ReactionSystem(rx, t, collect(X), collect(k))
Catalyst.Graph(rs)