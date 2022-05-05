using Catlab.Present
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs
using Catlab.Graphs
using Catlab.Theories
using Catlab.Graphics
import Catlab.Graphics.Graphviz.to_graphviz
#using Catlab.Graphics.Graphviz
import Catlab.WiringDiagrams: oapply
using Flux

base = @acset Graph begin
    V = 3
    E = 5
    src = [1,2,2,1,3]
    tgt = [2,2,3,1,3]
end
g = to_graphviz(base)

wide_boi = @acset Graph begin
    V = 5
    E = 11
    src = [1,1,1,2,3,4, 1, 2, 3, 4, 5]
    tgt = [2,3,4,5,5,5, 1, 2, 3, 4, 5]
end
wg = to_graphviz(wide_boi)
two = @acset Graph begin
    V = 2
    E = 2
    src = [1,2]
    tgt = [1,2]
end
length_boi = @acset Graph begin
    V = 4
    E = 7
    src = [1,2,3, 1, 2 , 3 , 4]
    tgt = [2,3,4, 1 , 2 , 3 , 4]
end
lg = to_graphviz(length_boi)

h1 = homomorphisms(length_boi, base)[2]
h2 = homomorphisms(wide_boi, base)[2]
h3 = homomorphisms(two, base)[5]

p = apex(pullback(h1,h2))
pg = to_graphviz(p)

if false 

#test = apex(coproduct(wide_boi, length_boi))
#testg = to_graphviz(test)=#

A, B, X, Y = Ob(FreeTracedMonoidalCategory, :A, :B, :X, :Y)
f = Hom(:f, otimes(X,A), otimes(X,B))

#base = trace(X, A, B, f)
#bd = Graphviz.to_graphviz(base)

show_diagram(d::WiringDiagram) = Graphviz.to_graphviz(
    add_junctions(d),
    orientation=LeftToRight,
    labels=false, label_attr=:xlabel,
    node_attrs=Graphviz.Attributes(
        :fontname => "Courier",
    ),
    edge_attrs=Graphviz.Attributes(
        :fontname => "Courier",
    )
)

@present Primitives(FreeCartesianCategory) begin
    T::Ob
    SepConv3x3::Hom(T, T)
    DepConv3x3::Hom(T, T)
    Conv1x1::Hom(T, T)
    MaxPool3x3::Hom(T, T)
    AvgPool3x3::Hom(T, T)
end

chain = @program Primitives (x::T) begin
    return Conv1x1(AvgPool3x3(SepConv3x3(DepConv3x3(x))))
end
chain_exp = to_hom_expr(FreeTracedMonoidalCategory, chain)

#h1 = homomorphisms(chain_exp, base)

struct NNDom
    nchannels::Vector{Int} # records how many channels for each layer in a parallel net
end

struct NeuralNet
    dom::NNDom
    codom::NNDom
    net
end

struct Layer
    l
end

#=@instance CartesianCategory{NNDom, NeuralNet} begin
    id(A::NNDom) = NeuralNet(A, A, x->x)
    dom(f::NeuralNet) = f.dom
    codom(f::NeuralNet) = f.codom

    compose(f::NeuralNet, g::NeuralNet) = begin
        @assert f.codom == g.dom
        try
            NeuralNet(Chain(f.net..., g.net...), f.dom, g.codom)
        catch try
            NeuralNet(Chain(f.net..., g.net), f.dom, g.codom)
        catch try
            NeuralNet(Chain(f.net, g.net...), f.dom, g.codom)
        catch try
            NeuralNet(Chain(f.net, g.net), f.dom, g.codom)
        catch
        end
        end
        end
        end
    end

    otimes(A::NNDom, B::NNDom) = NNDom(vcat(A.nchannels, B.nchannels))
    otimes(f::NeuralNet, g::NeuralNet) = begin
        
    end
end=#

#=net = @program Primitives (x::T) begin
    top_net = MaxPool3x3(SepConv3x3(DepConv3x3(x)))
    bot_net = MaxPool3x3(Conv1x1(x))
    y = [top_net, bot_net]
    return AvgPool3x3(y)
end
#net = pair(net...)

function nearest_common_neighbor(d::WiringDiagram, x::Int, y::Int)
    

function oapply_nn(d::WiringDiagram, gens::Dict)
    function parallel_helper(ids)


    function rec_helper(id)
        if length(outneighbors(d, id)) > 1
            Flux.Parallel(hcat, )

function oapply_nn(d::WiringDiagram, gens::Dict)
    function rec_helper(id, res)
        if length(outneighbors(d, id)) > 1
            for n in outneighbors(d, id)
                v = Vector()
                push!(res, v)
                rec_helper(n, v)
            end
        else
            for n in outneighbors(d, id)
                @assert haskey(gens, boxes(d)[id].value)
                push!(res, gens[boxes(d)[id].value])
                rec_helper(n, res)
            end
        end
    end
    res = Any[]
    rec_helper(input_id(d), res)
    return res
end

gens = Dict(
    :DepConv3x3 => 1,
    :SepConv3x3 => 2,
    :MaxPool3x3 => 3,
    :Conv1x1 => 4,
    :AvgPool3x3 => 5
)=#
end