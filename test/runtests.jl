using Test

@testset "param_est" begin
  include("param_est.jl")
end

@testset "ttest" begin
  include("ttest.jl")
end

@testset "Petri" begin
  include("Petri.jl")
end

@testset "ModelExploration" begin
  include("ModelExploration.jl")
end

