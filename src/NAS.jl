using Catlab.Present
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Theories
import Catlab.WiringDiagrams: oapply
using Flux

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
    return [top_net, bot_net]
end
#net = pair(net...)

#function oapply(d::WiringDiagram, gens::AbstractDict{S,L}) where {S, L <: Layer}


