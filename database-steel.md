# Databases and optimization: a match made in ~~heaven~~ steel
production
Sean L. Wu
2023-12-03

There’s a great paper [“Database structures for mathematical programming
models”](https://www.sciencedirect.com/science/article/abs/pii/S0167923697000079)
that describes how databases and optimziation models can interact to
make complex modeling tasks easier, using as a jumping off point a
project regarding steel production. Let’s figure it out.

## LP problem

The basic setup considering LPs of the form:

$$
\begin{align*}
\text{max      } & \sum_{j=1}^{n} c_{j}x_{j} \\
\text{subject to      } & l_{i}^{\text{row}} \leq \sum_{j=1}^{n} a_{ij}x_{j} \leq u_{i}^{\text{row}}, i=1,\dots,m \\
& l_{j}^{\text{col}} \leq x_{j} \leq u_{j}^{\text{col}}, j=1,\dots,n
\end{align*}
$$

It has a “relational” interpretation by defining the following sets:

- $\mathcal{I}=\{1,\dots,m\}$: set of constraints
- $\mathcal{J}=\{1,\dots,n\}$: set of decision variables

For each $i\in\mathcal{I}$ we have the data:

- $l_{i}^{\text{row}}$: lower limit on constraint $i$
- $u_{i}^{\text{row}}$: upper limit on constraint $i$

For each $j\in\mathcal{J}$ we have the data:

- $x_{j}$: actual value of the decision variable
- $c_{j}$: profit (if positive) or cost (if negative) per unit of the
  variable
- $l_{j}^{\text{col}}$: lower limit on variable $j$
- $u_{j}^{\text{col}}$: upper limit on variable $j$

We also define a set of constraint-variable pairs:

- $\mathcal{C}\subseteq\mathcal{I}\times\mathcal{J}$ such that
  $(i,j)\in\mathcal{C}$ means variable $j$ is used in constraint $i$
- $a_{ij}$ is the coefficient of variable $j$ in constraint $i$

Using the relational point of view, we can rewrite the LP as:

$$
\begin{align*}
\text{max      } & \sum_{j\in\mathcal{J}} c_{j}x_{j} \\
\text{subject to      } & l_{i}^{\text{row}} \leq \sum_{(i,j)\in\mathcal{C}} a_{ij}x_{j} \leq u_{i}^{\text{row}}, \forall i\in\mathcal{I} \\
& l_{j}^{\text{col}} \leq x_{j} \leq u_{j}^{\text{col}}, \forall j\in\mathcal{J}
\end{align*}
$$

If the middle sum is slightly ambigious, it means there is an outer loop
over $i\in\mathcal{I}$ and the inner summation is over $j\in\mathcal{J}$
associated to that $i$.

## Basic process model

The basic idea of the models is: raw materials enter, various
transformations to intermediate materials occur, and finished materials
leave. Profit is total revenue from finished products, minus costs of
raw materials and transformations.

Model has some stuff:

- Materials: all the physical stuff, including raw, intermediate, and
  finished materials. Each have data on:
  - cost per unit for buying, min/max quantities that can be bought
  - revenue per unit for selling, and min/quantities that can be sold
  - raw materials can only be bought, finished can only be sold,
    intermediates often neither. We can also have list of allowed
    conversions to other materials, with given yield and cost per unit
    converted.
- Facilities: places at which transformations occur. They may have data
  on:
  - min/max on overall capacity
  - min/max on total use of each input
  - min/max on total use of each output
- Activities: housed in facilities. Activities use and produce materials
  in certain proportions. Each activity at a facilities has the
  following data:
  - Amount of each input for unit of activity
  - Amount of each output for unit of activity
  - Cost per unit of activity
  - min/max on units of activity
  - units of activity that can be accomodated by one unit of the owning
    facility’s overall capacity

## First “algebraic” formulation

Two fundamental sets:

- $\mathcal{M}$: Materials
- $\mathcal{F}$: Facilities