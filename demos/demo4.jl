using Revise
using Test
using ModelExploration
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories
using Plots
#include("../src/Petri.jl")
include("../src/ModelSelection.jl")
using AlgebraicPetri: LabelledPetriNetUntyped


SIR = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
)
SIRS = LabelledPetriNet([:S, :I, :R],
  :inf => ((:S, :I)=>(:I, :I)),
  :rec => (:I=>:R),
  :sus => (:R=>:S),
)

SIR2 = LabelledPetriNet([:S1, :I1, :R1, :S2, :I2, :R2],
    :inf1 => ((:S1, :I1)=>(:I1, :I1)),
    :inf2 => ((:S2, :I2)=>(:I2, :I2)),
    :rec1 => (:I1=>:R1),
    :rec2 => (:I2=>:R2),
    :S12 =>  (:S1=>:S2),
    :S21 =>  (:S2=>:S1),
    :I12 =>  (:I1=>:I2),
    :I21 =>  (:I2=>:I1),
    :R12 =>  (:R1=>:R2),
    :R21 =>  (:R2=>:R1),
)

SIRS2 = LabelledPetriNet([:S1, :I1, :R1, :S2, :I2, :R2],
    :inf1 => ((:S1, :I1)=>(:I1, :I1)),
    :inf2 => ((:S2, :I2)=>(:I2, :I2)),
    :rec1 => (:I1=>:R1),
    :rec2 => (:I2=>:R2),
    :sus1 => (:R1=>:S1),
    :sus2 => (:R2=>:S2),
    :S12 =>  (:S1=>:S2),
    :S21 =>  (:S2=>:S1),
    :I12 =>  (:I1=>:I2),
    :I21 =>  (:I2=>:I1),
    :R12 =>  (:R1=>:R2),
    :R21 =>  (:R2=>:R1),
)


SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :death => (:I=>:D)
)


SIRSD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :sus => (:R=>:S),
    :death => (:I=>:D)
)


SIRD2 = LabelledPetriNet([:S1, :I1, :R1, :D1, :S2, :I2, :R2, :D2],
    :inf1 => ((:S1, :I1)=>(:I1, :I1)),
    :inf2 => ((:S2, :I2)=>(:I2, :I2)),
    :rec1 => (:I1=>:R1),
    :rec2 => (:I2=>:R2),
    :death1 => (:I1=>:D1),
    :death2 => (:I2=>:D2),
    :S12 =>  (:S1=>:S2),
    :S21 =>  (:S2=>:S1),
    :I12 =>  (:I1=>:I2),
    :I21 =>  (:I2=>:I1),
    :R12 =>  (:R1=>:R2),
    :R21 =>  (:R2=>:R1),
)

SIRSD2 = LabelledPetriNet([:S1, :I1, :R1, :D1, :S2, :I2, :R2, :D2],
    :inf1 => ((:S1, :I1)=>(:I1, :I1)),
    :inf2 => ((:S2, :I2)=>(:I2, :I2)),
    :rec1 => (:I1=>:R1),
    :rec2 => (:I2=>:R2),
    :sus1 => (:R1=>:S1),
    :sus2 => (:R2=>:S2),
    :death1 => (:I1=>:D1),
    :death2 => (:I2=>:D2),
    :S12 =>  (:S1=>:S2),
    :S21 =>  (:S2=>:S1),
    :I12 =>  (:I1=>:I2),
    :I21 =>  (:I2=>:I1),
    :R12 =>  (:R1=>:R2),
    :R21 =>  (:R2=>:R1),
)

models = [SIR, SIRS, SIRD, SIRSD, SIR2, SIRS2, SIRD2, SIRSD2];
tspan = (0.0, 50.0);
#p = repeat([1e-4], 12)
p = [.001, .002, .01, .02, .005, .01, .005, 0.01, 0.0025, .001, 0.001, 0.01];
u0 = [400.0, 10.0, 0.0, 0.0, 200.0, 0.0, 0.0, 0.0];

# Generate data
true_model = SIRD2;
sample_data, sample_times, prob, sol = generate_data(true_model, p, u0, tspan, 50);

label=reshape(map(string, true_model[:, :sname]), 1, ns(true_model))
plot(sol; label=label)
