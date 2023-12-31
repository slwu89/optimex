# IJKLM model
Sean L. Wu
2023-11-13

## The IJKLM model

An (in)famous [Julia Discourse
thread](https://discourse.julialang.org/t/performance-julia-jump-vs-python-pyomo/92044)
once asked how to make model building in JuMP faster compared to other
platforms (Pyomo, and I guess GAMS). This turned into a [blog post by
GAMS](https://www.gams.com/blog/2023/07/performance-in-optimization-models-a-comparative-analysis-of-gams-pyomo-gurobipy-and-jump/),
a company which produces a modeling language that plugs into various
commercial and open source solvers. The original formulation of the JuMP
model showed very bad performance, but that was due to the way the
author (employee of GAMS) was using Julia, rather than anything JuMP
specific. The JuMP dev team responded with their own [blog post on the
JuMP website](https://jump.dev/2023/07/20/gams-blog/), which was
announced on a [Discourse
thread](https://discourse.julialang.org/t/jump-developers-response-to-benchmarks-by-gams/101920).
They used a DataFrames.jl based solution that was very fast.

The backstory is filled with intrigue, no doubt, but I’m interested to
see if acsets can provide an alterative way to generate this model. The
original data is reproducible at
[justine18/performance_experiment](https://github.com/justine18/performance_experiment),
but I decided to just redo it in Julia as the writing/reading to/from
JSON files is a pain in the butt.

By the way, the model is given as:

$$\text{min} \ z = 1$$

$$\sum_{(j,k):(i,j,k) \in \mathcal{IJK}} \ \sum_{l:(j,k,l) \in \mathcal{JKL}} \ \sum_{m:(k,l,m) \in \mathcal{KLM}} x_{i,j,k,l,m} \ge 0 \hspace{1cm} \forall \ i \in \mathcal{I}$$

$$x_{i,j,k,l,m} \ge 0 \hspace{1cm} \forall \ (i,j,k) \in \mathcal{IJK}, l:(j,k,l) \in \mathcal{JKL}, m:(k,l,m) \in \mathcal{KLM} $$

The blog post calls subsets of Cartesian products “maps”, which seems to
be confused as a “map” is generally understood to be a function in math.
General subsets of products are known as “relations”.

## Data generation

First we load some packages we’ll need. `DataFrames` for data frames,
`Distributions` for sampling binomial random variates, `JuMP` to set up
the model, `HiGHS` for a solver. `Catlab` provides a new data structure,
the C-Set which will be compared to the data frames method.

``` julia
using DataFrames
using Distributions
using JuMP, HiGHS
using Catlab, DataMigrations
using BenchmarkTools, MarkdownTables
```

We first generate synthetic “data”. This should follow the data
generation as I understand it from the original repo. The probability of
all zeros with the given model sizes is incomprehensibly small but I
added a check for it anyway. Who knows.

``` julia
SampleBinomialVec = function(A,B,C,p=0.05)
    vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    while sum(vec) == 0
        vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    end
    return vec
end

n=100 # something large
m=20 # 20

# Sets IJKLM 
I = ["i$x" for x in 1:n]
J = ["j$x" for x in 1:m]
K = ["k$x" for x in 1:m]
L = ["l$x" for x in 1:m]
M = ["m$x" for x in 1:m]

# make IJK
IJK = DataFrame(Iterators.product(I,J,K))
rename!(IJK, [:i,:j,:k])
IJK.value = SampleBinomialVec(I,J,K)

# make JKL
JKL = DataFrame(Iterators.product(J,K,L))
rename!(JKL, [:j,:k,:l])
JKL.value = SampleBinomialVec(J,K,L)

# make KLM
KLM = DataFrame(Iterators.product(K,L,M))
rename!(KLM, [:k,:l,:m])
KLM.value = SampleBinomialVec(K,L,M)

# make the products just general sparse relations
IJK_sparse = [(x.i, x.j, x.k) for x in eachrow(IJK) if x.value == true]
JKL_sparse = [(x.j, x.k, x.l) for x in eachrow(JKL) if x.value == true]
KLM_sparse = [(x.k, x.l, x.m) for x in eachrow(KLM) if x.value == true]
```

## The “intuitive” formulation

As we know this is the slow one.

``` julia
@benchmark let 
    x_list = [
        (i, j, k, l, m)
        for (i, j, k) in IJK_sparse
        for (jj, kk, l) in JKL_sparse if jj == j && kk == k
        for (kkk, ll, m) in KLM_sparse if kkk == k && ll == l
    ]
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[x_list] >= 0)
    @constraint(
        model,
        [i in I], 
        sum(x[k] for k in x_list if k[1] == i) >= 0
    )
    optimize!(model)
end
```

    BenchmarkTools.Trial: 167 samples with 1 evaluation.
     Range (min … max):  27.004 ms … 35.219 ms  ┊ GC (min … max):  7.59% … 14.08%
     Time  (median):     29.646 ms              ┊ GC (median):    13.17%
     Time  (mean ± σ):   29.978 ms ±  1.554 ms  ┊ GC (mean ± σ):  12.29% ±  2.25%

                     █▃   ▃                                        
      ▅▁▄▆▃▃▅▁▁▄▃▄▃▅███▄▆▆█▆▇▃▆▃█▄▇▆▆▅▄▆▇▆▃▄▄▃▅▃▆▅▃▃▁▃▄▁▃▃▃▁▁▃▁▁▃ ▃
      27 ms           Histogram: frequency by time          34 ms <

     Memory estimate: 79.13 MiB, allocs estimate: 1679535.

## The DataFrames version

The fast one at the JuMP blog link.

``` julia
IJK_sparse_df = filter(x -> x.value == 1, IJK)
select!(IJK_sparse_df, Not(:value))

JKL_sparse_df = filter(x -> x.value == 1, JKL)
select!(JKL_sparse_df, Not(:value))

KLM_sparse_df = filter(x -> x.value == 1, KLM)
select!(KLM_sparse_df, Not(:value))

ijklm_df = DataFrames.innerjoin(
    DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
    KLM_sparse_df;
    on = [:k, :l],
)
```

Let’s benchmark it.

``` julia
@benchmark let
    ijklm = DataFrames.innerjoin(
        DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
        KLM_sparse_df;
        on = [:k, :l],
    )
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
    for df in DataFrames.groupby(ijklm, :i)
        @constraint(model, sum(df.x) >= 0)
    end
    optimize!(model)
end
```

    BenchmarkTools.Trial: 1105 samples with 1 evaluation.
     Range (min … max):  4.009 ms … 10.604 ms  ┊ GC (min … max): 0.00% … 39.24%
     Time  (median):     4.197 ms              ┊ GC (median):    0.00%
     Time  (mean ± σ):   4.521 ms ±  1.120 ms  ┊ GC (mean ± σ):  4.92% ±  9.81%

      ▅█▆▅▄▃▂▁                                                    
      ████████▆▅▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▆▅▆▆██▆▇▆▅ █
      4.01 ms      Histogram: log(frequency) by time     9.11 ms <

     Memory estimate: 2.73 MiB, allocs estimate: 23833.

## The acsets version

Acsets (Attributed C-Sets) are a nifty data structure coming from
applied category theory, but its not too far off to think of them just
as in-memory relational databases. They are provided in the
[ACSets.jl](https://github.com/AlgebraicJulia/ACSets.jl) library, and
imported and extended with machinery from applied category theory in
[Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl). One of the
advantages of using acsets and categorical machinery in general is that
they have a natural graphical presentation, which in many cases is very
readable and illuminating. To learn more about this data structure and
many more topics, you can consult the package documentation and also you
can begin with some [introductory articles on the AlgebraicJulia
blog](https://blog.algebraicjulia.org/post/2020/09/cset-graphs-1/).

We use `@present` to make a schema for the acset which will store the
data. Note that each “set” has turned into an object in the schema, and
that the relations are also objects. There are projection maps from the
relations into the sets which are involved in each relation. Each
relation also has an arrow into `IntAttr` which will store the results
of the binomial random draw that is used to “sparsify” the data.

``` julia
@present IJKLMSch(FreeSchema) begin
    (I,J,K,L,M,IJK,JKL,KLM)::Ob
    IJK_I::Hom(IJK,I)
    IJK_J::Hom(IJK,J)
    IJK_K::Hom(IJK,K)
    JKL_J::Hom(JKL,J)
    JKL_K::Hom(JKL,K)
    JKL_L::Hom(JKL,L)
    KLM_K::Hom(KLM,K)
    KLM_L::Hom(KLM,L)
    KLM_M::Hom(KLM,M)
    IntAttr::AttrType
    value_ijk::Attr(IJK,IntAttr)
    value_jkl::Attr(JKL,IntAttr)
    value_klm::Attr(KLM,IntAttr)
end

Catlab.to_graphviz(IJKLMSch, graph_attrs=Dict(:dpi=>"72",:ratio=>"expand",:size=>"8"))
```

![](gams_files/figure-commonmark/cell-7-output-1.svg)

Now we programatically generate the data type (and functions to work
with it) for our schema, and fill it with data. The code is verbose, but
we’re storing all the sets and relations in a single data structure.

``` julia
@acset_type IJKLMData(IJKLMSch, index=[:IJK_I,:IJK_J,:IJK_K,:JKL_J,:JKL_K,:JKL_L,:KLM_K,:KLM_L,:KLM_M])

# the basic sets
I = collect(1:n)
J = collect(1:m)
K = collect(1:m)
L = collect(1:m)
M = collect(1:m)

# make the data
ijklm_dat = IJKLMData{Int}()

add_parts!(ijklm_dat, :I, length(I))
add_parts!(ijklm_dat, :J, length(J))
add_parts!(ijklm_dat, :K, length(K))
add_parts!(ijklm_dat, :L, length(L))
add_parts!(ijklm_dat, :M, length(M))

# add the relations...first as product...then sparsify

# IJK
ijk_prod = Iterators.product(I,J,K)

add_parts!(
    ijklm_dat, 
    :IJK,
    length(ijk_prod),
    IJK_I=vec([e[1] for e in ijk_prod]),
    IJK_J=vec([e[2] for e in ijk_prod]),
    IJK_K=vec([e[3] for e in ijk_prod]),
    value_ijk=IJK.value
)

rem_parts!(
    ijklm_dat, 
    :IJK, 
    findall(ijklm_dat[:,:value_ijk] .== 0)
)

# JKL
jkl_prod = Iterators.product(J,K,L)

add_parts!(
    ijklm_dat, 
    :JKL,
    length(jkl_prod),
    JKL_J=vec([e[1] for e in jkl_prod]),
    JKL_K=vec([e[2] for e in jkl_prod]),
    JKL_L=vec([e[3] for e in jkl_prod]),
    value_jkl=JKL.value
)

rem_parts!(
    ijklm_dat, 
    :JKL, 
    findall(ijklm_dat[:,:value_jkl] .== 0)
)

# KLM
klm_prod = Iterators.product(K,L,M)

add_parts!(
    ijklm_dat, 
    :KLM,
    length(klm_prod),
    KLM_K=vec([e[1] for e in klm_prod]),
    KLM_L=vec([e[2] for e in klm_prod]),
    KLM_M=vec([e[3] for e in klm_prod]),
    value_klm=KLM.value
)

rem_parts!(
    ijklm_dat, 
    :KLM, 
    findall(ijklm_dat[:,:value_klm] .== 0)
)
```

### conjunctive queries on acsets

Now, the critical thing that the JuMP devs did to speed thing up was to
replace the for loops with 2 inner joins, to get the “paths” through the
relations. How to do this with acsets? Well one thing we can do is
execute a conjunctive query on the acset to get the same thing. This is
described in a [post at the AlgebraicJulia
blog](https://blog.algebraicjulia.org/post/2020/12/cset-conjunctive-queries/).

``` julia
connected_paths_query = @relation (i=i,j=j,k=k,l=l,m=m) begin
    IJK(IJK_I=i, IJK_J=j, IJK_K=k)
    JKL(JKL_J=j, JKL_K=k, JKL_L=l)
    KLM(KLM_K=k, KLM_L=l, KLM_M=m)
end

Catlab.to_graphviz(connected_paths_query, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:dpi=>"72",:size=>"3.5",:ratio=>"expand"))
```

![](gams_files/figure-commonmark/cell-9-output-1.svg)

While the blog post should be consulted for a complete explanation, the
conjunctive query is expressed using an undirected wiring diagram (UWD)
which is visualized using Catlab. Nodes (labeled ovals) in the UWD
correspond to tables (primary keys) in the acset. Junctions (labeled
dots) correspond to variables. Ports, which are unlabed in this
graphical depiction, are where wires connect junctions to nodes. These
correspond to columns of the table they are connected to. Outer ports,
which are wires that run “off the page”, are the columns of the table
that will be returned as the result of the query. The rows that are
returned from the query come from filtering the Cartesian product of the
tables (nodes) such that variables in columns match according to ports
that share a junction. In this case, it is equivalent to the two inner
joins done above using DataFrames.

The JuMP blog post notes that while the data frames version doesn’t
resemble the nested summation it is arguably just as readable,
especially if the columns were related to the process that was being
modeled. I suggest that the acsets version is also just as readable, if
not more, as the data schema and query diagram directly represent the
data that parameterizes the optimization model.

The query is evaluated on the acset by computing its limit, which is a
kind of generalization of meet. We can confirm that both the acsets and
data frame methods return the same number of rows, and look at the
output.

``` julia
ijklm_query = query(ijklm_dat, connected_paths_query)
size(ijklm_query) == size(ijklm_df)
```

    true

``` julia
ijklm_query[1:5,:] |> markdown_table
```

| i   | j   | k   | l   | m   |
|-----|-----|-----|-----|-----|
| 73  | 10  | 16  | 16  | 11  |
| 97  | 15  | 15  | 14  | 13  |
| 97  | 15  | 15  | 12  | 1   |
| 69  | 19  | 3   | 10  | 15  |
| 69  | 19  | 3   | 20  | 4   |

Now that we know they are equal, we can go ahead and see how fast the
acsets version is. The fact that the acsets based query is right on the
tails of the `DataFrames` version is a performance win for the acsets
library, as it is usually a generic conjunctive query engine across very
general data structures (i.e., acsets are in general much more complex
than a single dataframe, due to presence of multiple tables connected
via foreign keys).

``` julia
@benchmark let
    ijklm = query(ijklm_dat, connected_paths_query)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
    for df in DataFrames.groupby(ijklm, :i)
        @constraint(model, sum(df.x) >= 0)
    end
    optimize!(model)
end
```

    BenchmarkTools.Trial: 1024 samples with 1 evaluation.
     Range (min … max):  4.359 ms … 10.520 ms  ┊ GC (min … max): 0.00% … 38.41%
     Time  (median):     4.482 ms              ┊ GC (median):    0.00%
     Time  (mean ± σ):   4.881 ms ±  1.223 ms  ┊ GC (mean ± σ):  5.71% ± 10.57%

      ▆█▄▂▂▂                                                      
      ███████▆▄▆▄▅▅▄▄▄▁▄▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▄▇█▆▇█▇▇▆▄ ▇
      4.36 ms      Histogram: log(frequency) by time     9.34 ms <

     Memory estimate: 3.47 MiB, allocs estimate: 30406.

### data migrations of acsets

While the query execution in the previous section is already quite
useful for practical application, it does have one downside, and that is
the type of the return object is a `DataFrame`. While this is
appropriate for many cases, one of the benefits of acsets is that, from
one point of view, they are in-memory relational databases, and are
therefore capable of representing data that is not possible to do with a
single table. Therefore, it would be nice if one could execute a *data
migration* from one type of acset, where “type” means the schema, to
another, that is able to carry along further data we need. For more
details on contravariant data migration, please see Evan Patterson’s
Topos Institute colloquium presentation [“Categories of diagrams in data
migration and computational
physics”](https://www.youtube.com/live/Ra-PLnog_M0?si=1_ex4wLud2hSR7be).

In this context, if the “set” objects (`I`, `J`, etc) were further
connected to other tables, maybe, say, a list of suppliers, or a list of
materials, or even a process graph of downstream work, it would be
inconvenient at least, if we lost that relational information during a
query. In that case, we’d really want to return *another acset* on a
different schema that is precisely the right shape for what we want to
do.

In this simple case, we have a quite simple schema which looks like a
starfish. However, if some of the “set” objects had further morphisms to
other objects (i.e., other tables), the advantages of using data
migration become apparent.

``` julia
@present IJKLMRelSch(FreeSchema) begin
    (IJKLM,I,J,K,L,M)::Ob
    i::Hom(IJKLM,I)
    j::Hom(IJKLM,J)
    k::Hom(IJKLM,K)
    l::Hom(IJKLM,L)
    m::Hom(IJKLM,M)
end

@acset_type IJKLMRelType(IJKLMRelSch)

Catlab.to_graphviz(IJKLMRelSch, graph_attrs=Dict(:dpi=>"72",:ratio=>"expand",:size=>"3.5"))
```

![](gams_files/figure-commonmark/cell-13-output-1.svg)

Now we formulate the data migration, using tools from the
[AlgebraicJulia/DataMigrations.jl](https://github.com/AlgebraicJulia/DataMigrations.jl)
package. While we will not be able to rigorously explain data migration
here, if one has $C$-Set (instance of data on schema $C$) and wants to
migrate it to a $D$-Set, a data migration functor $F$ needs to be
specified.

Here, $C$ is our schema `IJKLMSch` and $D$ is `IJKLMRelSch`. The functor
$F$ is a mapping from $D$ to the category of diagrams on $C$; formally
we denote it $F:D\rightarrow \text{Diag}^{\text{op}}(C)$. Each object in
$D$ gets assigned a diagram into $C$, and morphisms in $D$ get assigned
to contravariant morphisms of diagrams.

``` julia
M = @migration IJKLMRelSch IJKLMSch begin
    IJKLM => @join begin
        ijk::IJK
        jkl::JKL
        klm::KLM
        i::I
        j::J
        k::K
        l::L
        m::M
        IJK_I(ijk) == i
        IJK_J(ijk) == j
        JKL_J(jkl) == j
        IJK_K(ijk) == k
        JKL_K(jkl) == k
        KLM_K(klm) == k
        JKL_L(jkl) == l
        KLM_L(klm) == l
        KLM_M(klm) == m
    end
    I => I
    J => J
    K => K
    L => L
    M => M
    i => i
    j => j
    k => k
    l => l
    m => m
end
```

    DataMigration{Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, TypeCat{SimpleDiagram{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, D} where D<:(Functor{<:Category{Ob, Hom, Catlab.CategoricalAlgebra.FinCats.FinCatSize} where {Ob, Hom}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}), Catlab.CategoricalAlgebra.Diagrams.SimpleDiagramHom{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinTransformationMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:id}}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:generator}}}}}, Dict{Symbol, SimpleDiagram{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, D} where D<:(Functor{<:Category{Ob, Hom, Catlab.CategoricalAlgebra.FinCats.FinCatSize} where {Ob, Hom}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}})}, Dict{Symbol, Catlab.CategoricalAlgebra.Diagrams.SimpleDiagramHom{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinTransformationMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:id}}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:generator}}}}}}, Dict{Any, Union{}}}(FinDomFunctor(Dict{Symbol, SimpleDiagram{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, D} where D<:(Functor{<:Category{Ob, Hom, Catlab.CategoricalAlgebra.FinCats.FinCatSize} where {Ob, Hom}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}})}(:IJKLM => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:l => L, :m => M, :klm => KLM, :k => K, :ijk => IJK, :j => J, :jkl => JKL, :i => I), Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:generator}}(Symbol("##unnamedhom#9") => KLM_M, Symbol("##unnamedhom#5") => JKL_K, Symbol("##unnamedhom#1") => IJK_I, Symbol("##unnamedhom#4") => IJK_K, Symbol("##unnamedhom#3") => JKL_J, Symbol("##unnamedhom#2") => IJK_J, Symbol("##unnamedhom#6") => KLM_K, Symbol("##unnamedhom#8") => KLM_L, Symbol("##unnamedhom#7") => JKL_L), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[ijk, jkl, klm, i, j, k, l, m], Hom = Hom{:generator}[##unnamedhom#1, ##unnamedhom#2, ##unnamedhom#3, ##unnamedhom#4, ##unnamedhom#5, ##unnamedhom#6, ##unnamedhom#7, ##unnamedhom#8, ##unnamedhom#9], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:j=>(:Ob=>5), Symbol("##unnamedhom#4")=>(:Hom=>4), Symbol("##unnamedhom#8")=>(:Hom=>8), Symbol("##unnamedhom#9")=>(:Hom=>9), :l=>(:Ob=>7), :jkl=>(:Ob=>2), Symbol("##unnamedhom#7")=>(:Hom=>7), Symbol("##unnamedhom#5")=>(:Hom=>5), :ijk=>(:Ob=>1), :klm=>(:Ob=>3)…), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[])))), :I => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:I => I), Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}(), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I], Hom = Hom{:generator}[], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:I=>(:Ob=>1)), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[])))), :M => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:M => M), Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}(), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[M], Hom = Hom{:generator}[], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:M=>(:Ob=>1)), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[])))), :J => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:J => J), Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}(), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[J], Hom = Hom{:generator}[], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:J=>(:Ob=>1)), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[])))), :K => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:K => K), Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}(), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[K], Hom = Hom{:generator}[], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:K=>(:Ob=>1)), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[])))), :L => Diagram{op}(FinFunctor(Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}(:L => L), Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}(), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[L], Hom = Hom{:generator}[], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:L=>(:Ob=>1)), Pair[])), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[I, J, K, L, M, IJK, JKL, KLM], Hom = Hom{:generator}[IJK_I, IJK_J, IJK_K, JKL_J, JKL_K, JKL_L, KLM_K, KLM_L, KLM_M], AttrType = AttrType{:generator}[IntAttr], Attr = Attr{:generator}[value_ijk, value_jkl, value_klm]), Dict(:KLM=>(:Ob=>8), :IJK_J=>(:Hom=>2), :IntAttr=>(:AttrType=>1), :JKL=>(:Ob=>7), :JKL_J=>(:Hom=>4), :IJK=>(:Ob=>6), :K=>(:Ob=>3), :KLM_L=>(:Hom=>8), :JKL_K=>(:Hom=>5), :M=>(:Ob=>5)…), Pair[]))))), Dict{Symbol, Catlab.CategoricalAlgebra.Diagrams.SimpleDiagramHom{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinTransformationMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:id}}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:generator}}}}}(:l => DiagramHom{op}([(l, id(L))], [], FinFunctor(Dict{Symbol, Ob{:generator}}(:l=>L, :m=>M, :klm=>KLM, :k=>K, :ijk=>IJK, :j=>J, :jkl=>JKL, :i=>I), Dict{Symbol, Hom{:generator}}(Symbol("##unnamedhom#9")=>KLM_M, Symbol("##unnamedhom#5")=>JKL_K, Symbol("##unnamedhom#1")=>IJK_I, Symbol("##unnamedhom#4")=>IJK_K, Symbol("##unnamedhom#3")=>JKL_J, Symbol("##unnamedhom#2")=>IJK_J, Symbol("##unnamedhom#6")=>KLM_K, Symbol("##unnamedhom#8")=>KLM_L, Symbol("##unnamedhom#7")=>JKL_L), …), FinFunctor(Dict{Symbol, Ob{:generator}}(:L=>L), Dict{Symbol, Union{Attr, AttrType, Hom}}(), …)), :m => DiagramHom{op}([(m, id(M))], [], FinFunctor(Dict{Symbol, Ob{:generator}}(:l=>L, :m=>M, :klm=>KLM, :k=>K, :ijk=>IJK, :j=>J, :jkl=>JKL, :i=>I), Dict{Symbol, Hom{:generator}}(Symbol("##unnamedhom#9")=>KLM_M, Symbol("##unnamedhom#5")=>JKL_K, Symbol("##unnamedhom#1")=>IJK_I, Symbol("##unnamedhom#4")=>IJK_K, Symbol("##unnamedhom#3")=>JKL_J, Symbol("##unnamedhom#2")=>IJK_J, Symbol("##unnamedhom#6")=>KLM_K, Symbol("##unnamedhom#8")=>KLM_L, Symbol("##unnamedhom#7")=>JKL_L), …), FinFunctor(Dict{Symbol, Ob{:generator}}(:M=>M), Dict{Symbol, Union{Attr, AttrType, Hom}}(), …)), :k => DiagramHom{op}([(k, id(K))], [], FinFunctor(Dict{Symbol, Ob{:generator}}(:l=>L, :m=>M, :klm=>KLM, :k=>K, :ijk=>IJK, :j=>J, :jkl=>JKL, :i=>I), Dict{Symbol, Hom{:generator}}(Symbol("##unnamedhom#9")=>KLM_M, Symbol("##unnamedhom#5")=>JKL_K, Symbol("##unnamedhom#1")=>IJK_I, Symbol("##unnamedhom#4")=>IJK_K, Symbol("##unnamedhom#3")=>JKL_J, Symbol("##unnamedhom#2")=>IJK_J, Symbol("##unnamedhom#6")=>KLM_K, Symbol("##unnamedhom#8")=>KLM_L, Symbol("##unnamedhom#7")=>JKL_L), …), FinFunctor(Dict{Symbol, Ob{:generator}}(:K=>K), Dict{Symbol, Union{Attr, AttrType, Hom}}(), …)), :j => DiagramHom{op}([(j, id(J))], [], FinFunctor(Dict{Symbol, Ob{:generator}}(:l=>L, :m=>M, :klm=>KLM, :k=>K, :ijk=>IJK, :j=>J, :jkl=>JKL, :i=>I), Dict{Symbol, Hom{:generator}}(Symbol("##unnamedhom#9")=>KLM_M, Symbol("##unnamedhom#5")=>JKL_K, Symbol("##unnamedhom#1")=>IJK_I, Symbol("##unnamedhom#4")=>IJK_K, Symbol("##unnamedhom#3")=>JKL_J, Symbol("##unnamedhom#2")=>IJK_J, Symbol("##unnamedhom#6")=>KLM_K, Symbol("##unnamedhom#8")=>KLM_L, Symbol("##unnamedhom#7")=>JKL_L), …), FinFunctor(Dict{Symbol, Ob{:generator}}(:J=>J), Dict{Symbol, Union{Attr, AttrType, Hom}}(), …)), :i => DiagramHom{op}([(i, id(I))], [], FinFunctor(Dict{Symbol, Ob{:generator}}(:l=>L, :m=>M, :klm=>KLM, :k=>K, :ijk=>IJK, :j=>J, :jkl=>JKL, :i=>I), Dict{Symbol, Hom{:generator}}(Symbol("##unnamedhom#9")=>KLM_M, Symbol("##unnamedhom#5")=>JKL_K, Symbol("##unnamedhom#1")=>IJK_I, Symbol("##unnamedhom#4")=>IJK_K, Symbol("##unnamedhom#3")=>JKL_J, Symbol("##unnamedhom#2")=>IJK_J, Symbol("##unnamedhom#6")=>KLM_K, Symbol("##unnamedhom#8")=>KLM_L, Symbol("##unnamedhom#7")=>JKL_L), …), FinFunctor(Dict{Symbol, Ob{:generator}}(:I=>I), Dict{Symbol, Union{Attr, AttrType, Hom}}(), …))), FinCat(Presentation{T, Symbol}(Catlab.Theories.FreeSchema, (Ob = Ob{:generator}[IJKLM, I, J, K, L, M], Hom = Hom{:generator}[i, j, k, l, m], AttrType = AttrType{:generator}[], Attr = Attr{:generator}[]), Dict(:j=>(:Hom=>2), :l=>(:Hom=>4), :IJKLM=>(:Ob=>1), :K=>(:Ob=>4), :M=>(:Ob=>6), :k=>(:Hom=>3), :m=>(:Hom=>5), :I=>(:Ob=>2), :J=>(:Ob=>3), :L=>(:Ob=>5)…), Pair[])), TypeCat(SimpleDiagram{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, D} where D<:(Functor{<:Category{Ob, Hom, Catlab.CategoricalAlgebra.FinCats.FinCatSize} where {Ob, Hom}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}), Catlab.CategoricalAlgebra.Diagrams.SimpleDiagramHom{GATlab.Stdlib.StdModels.op, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinTransformationMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Any, Any}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:id}}}, Catlab.CategoricalAlgebra.FinCats.FinDomFunctorMap{Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Catlab.CategoricalAlgebra.FinCats.FinCatPresentation{Catlab.Theories.ThSchema.Meta.T, Union{Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Ob}, Union{Catlab.Theories.FreeSchema.Attr, Catlab.Theories.FreeSchema.AttrType, Catlab.Theories.FreeSchema.Hom}}, Dict{Symbol, Catlab.Theories.FreeSchema.Ob{:generator}}, Dict{Symbol, Catlab.Theories.FreeSchema.Hom{:generator}}}})), Dict{Any, Union{}}())

A diagram is itself a functor $D:J\rightarrow C$, where $J$ is (usually)
a small category, and $D$ will point at some instance of the diagram in
$C$. We can plot what the largest diagram looks like, that which object
`IJKLM` in $D$ is mapped to. Note the similarity to the conjunctive
query visualized as a UWD previously. In particular, note that
“relation” elements must agree upon the relevant “set” elements via
their morphisms. The object in $C$ that each object in $J$ corresponds
to is given by the text after the colon in the relevant node.

``` julia
F = functor(M)
to_graphviz(F.ob_map[:IJKLM],node_labels=true)
```

![](gams_files/figure-commonmark/cell-15-output-1.svg)

Because of the simplicity of the schema `IJKLMRelSch`, the contravariant
morphisms of diagrams simply pick out the object in $D$ associated with
the source of the morphism. Likewise, the natural transformation part of
morphisms of diagrams simply selects for each object its identity
morphism.

We run the data migration to move data from the schema `IJKLMSch` to
`IJKLMRelSch` using the function `migrate`, and check that the result
has the same number of records as other methods.

``` julia
ijklm_migrate_acset = migrate(IJKLMRelType, ijklm_dat, M)
nparts(ijklm_migrate_acset, :IJKLM) == size(ijklm_query,1)
```

    true

Let’s look at the first few rows.

``` julia
pretty_tables(ijklm_migrate_acset, tables=[:IJKLM], max_num_of_rows=5)
```

    ┌───────┬────┬────┬────┬────┬────┐
    │ IJKLM │  i │  j │  k │  l │  m │
    ├───────┼────┼────┼────┼────┼────┤
    │     1 │ 73 │ 10 │ 16 │ 16 │ 11 │
    │     2 │ 97 │ 15 │ 15 │ 14 │ 13 │
    │     3 │ 97 │ 15 │ 15 │ 12 │  1 │
    │     4 │ 69 │ 19 │  3 │ 10 │ 15 │
    │     5 │ 69 │ 19 │  3 │ 20 │  4 │
    └───────┴────┴────┴────┴────┴────┘

Once again, let’s benchmark:

``` julia
@benchmark let
    ijklm = migrate(IJKLMRelType, ijklm_dat, M)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[parts(ijklm,:IJKLM)] >= 0)
    for i in parts(ijklm,:I)
        @constraint(model, sum(x[incident(ijklm,i,:i)]) >= 0)
    end
    optimize!(model)
end
```

    BenchmarkTools.Trial: 514 samples with 1 evaluation.
     Range (min … max):  8.592 ms … 21.361 ms  ┊ GC (min … max): 0.00% … 22.67%
     Time  (median):     8.991 ms              ┊ GC (median):    0.00%
     Time  (mean ± σ):   9.723 ms ±  1.594 ms  ┊ GC (mean ± σ):  6.39% ±  9.62%

      ██▅▅▅▆▃▁▂▁                           ▂▄▃ ▁                  
      ███████████▇▁▆▄▇█▁▄▅▅▄▁▁▁▇▁▁▁▄▁▄▁▁▄▆▇████████▅▄▁▁▆▁▁▄▄▁▄▁▄ █
      8.59 ms      Histogram: log(frequency) by time     13.8 ms <

     Memory estimate: 9.34 MiB, allocs estimate: 66729.
