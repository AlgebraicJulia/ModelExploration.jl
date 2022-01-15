# ModelExploration.jl
Leveraging AlgebraicJulia to provide an interface for scientists to explore spaces of models.

## Basic idea

Here, we abstract the notion of a *model* and a *search space*. For us, a model is an instance of an `ACSet`.

- The first key notion is that of a `Generator`: something that can generate a sequence of models.
   - <img src="/img/gen.png" width=50% height=50%>
   - There are many different ways we can imagine specifying this. E.g.:
       - One could literally give a finite sequence of precomputed models
           - <img src="/img/preorder.png" width=50% height=50%>
       - One could write a function of type `ℕ ⟶ ACSet` that induces an infinite sequence
           - <img src="/img/Ncity.png" width=50% height=50%>
       - One could enumerate all possible models
       - One could specify some rewrite rules and initial model
    - Call these *primitive* `Generators`
- A *composite* `Generator` can be either "addition-like" (`AddLayer`) or "multiplication-like" (`MulLayer`).
    - Composite `Generators` depend on other `Generators`.
    - Thus, the overall model space is specified by a set of `Generators` (whose dependencies must form a *rooted* DAG: the root being the top-level model).
- An `AddLayer` is specified by an undirected wiring diagram (UWD).
   - <img src="/img/uwd.png" width=50% height=50%>
   - UWDs have `Boxes`, `Junctions`, `Ports`, and `Wires`.[^1]
   - Each `Box` represents a placeholder for submodels to be filled in (by a `Generator`).
   - Each `Junction` is given by an `ACSet`.
       - This characterizes an overlap along which models are to be glued.
       - This is optional to specify. It defaults to the empty `ACSet`.
   - `Ports` live on `Boxes`. `Wires` connect `Junctions` and `Ports`
       - A `Wire` says that, given a model generated by the `Port`'s `Box`, it is to be glued to all other models connected to that `Junction`.
    - The model is glued together via a construction called a *pushout*
       - <img src="/img/pushout.png" width=50% height=50%>
       - It requires a *homomorphism* from the `ACSet` of the `Junction` into the model that is generated by the `Box`.
            - If there are *multiple* homomorphisms possible, we pick one at random.
            - We'll allow constraints on each `Port` to help guide this selection process
                - Example: we're making models with data flowing on wires and through functions.
                - We want to glue two submodels along a wire, *head-to-tail*.
                - We put a constraint on one port ("this wire is connected *only* to a function input") and the other port ("this wire is connected *only* to a function output").
            - We really want there to be *some* homomorphism, so we implicitly create a constraint on the `Generator` referred to by the `Box`.
               - Any model it produces must satisfy the `Box`'s interface.
               - For example: if models `[A,B,C,D,E,...]` are a stream of models that would generated by a particular `Generator`, but `B` and `D` don't satisfy a given interface, then the `AddLayer` will only see the sequence `[A,C,E,...]`
- A `MulLayer` is specified in terms of a set of dimensions, which are actually `Generators`.
   - This is in analogy to grid-search. Along one dimension you have models `[A B C]`, and along another dimension you have models `[1;2]`, which induces a grid of product models `[A1 B1 C1; A2 B2 C2]`.
        - This resulting grid can be traversed from the base point as a *preorder* (i.e. radiating outward from the origin: `(A1)->(B1,A2)->(B2,C1)->(C2)`).
   - However, we can generalize this:
        - Rather than have each dimension be a linear sequence of models, it can be preorder itself
            - The resulting product space still has the structure of a preorder.
        - We also could (but don't yet) give distances to the edges of the preorder and get a *metric space*.
            - Now, the analogy of 'radiating outward from the origin' becomes more literal.
   - Given a set of dimensions (i.e. `Generators`) and a choice for a model along each dimension, we construct a product model via *pullback*.
        - We need to interpret the models along each dimension as slices over a common base in order to take a pullback. This slice is optionally part of the data of the `MulLayer` (default: terminal object).
        - Suppose we slice over a particular A. If there are multiple homomorphisms X->A, we pick one at random.
        - Implicitly there is constraint on the generator output that there exist at least one homomorphism.
- Every `Generator` can be equipped with a `Loss` function, which evaluates the generated models against some criterion and possibly directs search in productive directions.
   - For example, imagine a task that we know has the structure of having three subtasks whose answers are combined to get the final result.
   - If we have criteria for both the final result *and* the subtasks, we can generate a sequence of models using four different `Loss` functions, one at the top level `AddLayer`, and one for each of the three `Generators` that fit into that top layer's three `Boxes`.
   - A `Generator` can have a stopping criterion that halts the sequence of models based on the loss function.
       - Because models are getting more and more complex, we want to stop as soon as we get the functionality we need.
- We can further constrain the outputs of `Generators` either by mere `Filters` or by `Chase` constraints.

As an algebraic data type, then, the required data is:
```
SearchSpace := {schema :: ACSetSchema,  gens :: [Generator]}

Generator := {name   :: Symbol,
              gen    :: Gen,
              constr :: [Constraint],
              loss   :: Maybe LossFn}

Gen := PrimitiveGen | CompositeGen

PrimitiveGen := ExplicitPreorder | RewriteRules
                | FreeGeneration | FunctionGenerated
                | etc.

CompositeGen := AddLayer | MulLayer

AddLayer := {pattern :: UWD}

UWD := {boxes :: [Box], ports :: [Port],
        junctions :: [Junction], wires :: [Wire]}
  Box      := {gen_id :: Symbol}
  Port     := {box_id :: Int, constraint :: [InterfaceConstraint]}
  Junction := {overlap :: Maybe ACSet}
  Wire     := {junc_id :: Int, port_id :: Int}

MulLayer   := {dim_gen_ids :: [Symbol], slice :: Maybe ACSet}

LossFn     := {fn :: RealValuedFunction, stop :: Maybe StopCriteria}
Constraint := Filter | Chase
```

This might seem overwhelming, but the full complexity is not needed in every case. An epidemiology model might be well generated from just one `MulLayer` containing three primitive `Dimensions`, and a circuit might be generated well from one `AddLayer`. It's likely not needed to have constraints on layers or to constrain the interfaces used in `AddLayers`.


[^1]: Normally UWDs also have 'outer ports', which allow you to know ahead of time what `Boxes` the result can be plugged into. However, for us, checking whether a particular submodel fits into a `Box` with a particular interface is something done at 'run-time' rather than 'compile-time', as it were. More details will come below.

## Three Examples

### Optimal boolean circuits
Given truth table find shortest formula. We have strict stopping criteria.

Possible schema:
AND -> -> -> 𝔹ool ⇇ NOT
        Input↗  ↖Output

CONSTRAINT:
- \# input = `n`, \# output = `1`  (CAN BE ENFORCED WITH CHASE)
- output vertex/wire should not be the input to any function
data: opaque function from 𝔹oolⁿ ⟶ 𝔹ool

### Epidemiology models
Given experimental data, we want to find the model that best explains it.
Our model space is Petri Nets that have up to two inputs and two outputs.


### Neural architecture search
Given a dataset of `n` features, we want to learn a function `ℝⁿ ⟶ ℝ` that fits the data without overfitting it. Fixing the method of training the model, the free parameter that decides what our function will be is the architecture of the network. The model space has networks with `n` inputs and one output.


### (Optional) Resistor networks
Oscilloscope data + graphs of resistors. Kind of a combination of the first two.


