# # ![Logo](docs/src/assets/logo.png)  ModelExploration.jl [![Documentation](https://github.com/kris-brown/ModelExploration.jl/workflows/Documentation/badge.svg)](https://kris-brown.github.io/ModelExploration.jl/dev/) ![Tests](https://github.com/kris-brown/ModelExploration.jl/workflows/Tests/badge.svg)

This package leverages [AlgebraicJulia](https://www.algebraicjulia.org/) to provide an interface for scientists to explore spaces of models.

Many search / optimization problems can be thought of as “exploring a space of models.” E.g.:

- We have a notion of what a chemical reaction network is and want to find one that fits some experimental data
- We have a notion of an electrical circuit and want to explore candidates that have the same voltage response to some black box we are probing.
- Neural networks depend heavily on their architecture - can we explore the space of possible architectures in a meaningful way to help us find better ones?

We address this process generically by formulating a language for model exploration. This language accomodates arbitrary kinds of “primitive” model spaces and instead is focused on how we can compose smaller spaces into larger ones, as described in our [ACT 2022 paper](https://arxiv.org/abs/2206.08755).


## NOTE
This library is currently under active development, and so is not yet at a
point where a constant API/behavior can be assumed. That being said, if this
project looks interesting/relevant please contact us and
[let us know](https://www.algebraicjulia.org/#contributing)!