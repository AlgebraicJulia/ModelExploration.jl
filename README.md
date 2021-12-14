# ModelExploration.jl
Leveraging AlgebraicJulia to provide an interface for scientists to explore spaces of models.

## Motivation

Here, we abstract the notion of a model and a search space. For us, a model is an instance of a ACSet. A model space is some subset of the possible instances for a fixed ACSet schema. We think of exploring this space indexed along user-specified *dimensions*.
- The most basic instance of this is a grid search, where each dimension is a finite number of points and we search among all combinations of all features.
- If the dimensions are *linearly ordered*, we can order the search space: e.g., for a three-dimensional search space with each dimension being ℕ, we would search: 000, 100, 010, 001, 110, 101, 011, 200, etc. Thus, we have a means of searching *infinite* search spaces (provided we have some *stopping criterion*).
- We can relax the constraint of the dimensions being *linearly* ordered. If each dimension is a preorder, we still induce an ordered search space that we can iteratively explore. For example, if one dimension is a square (one corner as bottom and the opposite corner as top) and we take the product with the dimension `1 ⟶ 2`, then the result is a cube (again, one corner as bottom and the opposite as top).
- We can think of the edges of the preorder as weighted (the above scenario being the case of all edges weighted equally). This allows us to think of each dimension not as a preorder but as a `Lawvere metric space`. We then get a product metric space, which gives us a more refined notion of the 'nearest' point yet to be explored.

At each point in the search space, indexed by a point for each dimension, we need a composite model that combines the models of each of the indices. This can be generically performed for any ACSet via pullback.

The required domain-specific pieces of data to run this procedure are:
- Defining the search space:
    - the ACSet schema for the class of models being explored
    - (optionally) any constraints on this class of models (e.g. objects of a slice category)
- Exploring the search space
    - an evaluation function for models in the search space (returning some metric of how far off from ideal the model is + (optionally) in what way it is off)
    - stopping criteria: when is a model 'good enough'?

## Three Examples

### Optimal boolean circuits
Given truth table find shortest formula. We have strict stopping criteria.

Possible schema:
AND -> -> -> 𝔹ool ⇇ NOT
        Input↗  ↖Output

CONSTRAINT IDEAS:
- # input = `n`, # output = `1`  (CAN BE ENFORCED WITH CHASE)
- output vertex/wire should not be the input to any function
- what do to about graphs that don't feed the output anything?

- initial model just input and output points, no other wires/boxes?
- TODO brainstorm dimensions

CONSTRAINT TYPES:
- Chase/slice === good (transparent)
- BOOLEAN FUNCTION FILTER / REPAIR FUNCTION (not transparent, maybe needed)

data: opaque function from 𝔹oolⁿ ⟶ 𝔹ool

E-graph approach to optimize circuits? Need to state logical laws in diagramatic form.

### Epidemiology models
Given experimental data, we want to find the model that best explains it.
Our model space is Petri Nets that have up to two inputs and two outputs.




### Neural architecture search
Given a dataset of `n` features, we want to learn a function `ℝⁿ ⟶ ℝ` that fits the data without overfitting it. Fixing the method of training the model, the free parameter that decides what our function will be is the architecture of the network. The model space has networks with `n` inputs and one output.


### (Optional) Resistor networks
Oscilloscope data + graphs of resistors. Kind of a combination of the first two.
