using Test

#=@testset "ModelExploration" begin
  include("ModelExploration.jl")
end
  
@testset "CatLimits" begin
  include("CatLimits.jl")
end=#


#=@testset "Chase" begin
  include("Chase.jl")
end

@testset "Currying" begin
  include("Curry.jl")
end=#

@testset "DiagLimits.jl" begin
  include("DiagLimits.jl")
end

