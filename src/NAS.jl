using Catlab.Present
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Theories
using Catlab.Graphics
using Catlab.Graphics.Graphviz
import Catlab.WiringDiagrams: oapply
using Flux

show_diagram(d::WiringDiagram) = to_graphviz(
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

net = @program Primitives (x::T) begin
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
)
