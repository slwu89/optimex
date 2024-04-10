using Catlab, JuMP

@present SchGraphJuMP <: SchGraph begin
    Var::AttrType
    var::Attr(V,Var)
end

@abstract_acset_type AbstractGraphJuMP <: AbstractGraph

@acset_type GraphJuMP(SchGraphJuMP, index=[:src,:tgt]) <: AbstractGraphJuMP

graphopt = @acset GraphJuMP{JuMP.Containers.DenseAxisArray} begin
    V = 3
    E = 2
    src = [1,2]
    tgt = [2,3]
end

jumpmod = JuMP.Model()

for v in vertices(graphopt)
    var = @variable(
        jumpmod,
        [t ∈ 5:10],
    )
    graphopt[v,:var] = var
end

graphopt[2, (:tgt,:var)]
graphopt[2, [:tgt,:var]]
graphopt[3, :var]

subpart(graphopt, 2, (:tgt,:var))
subpart(graphopt, 3, (:var,))
subpart(graphopt, 3, :var)

subpart(graphopt, 2, [:tgt,:var])
subpart(graphopt, 3, [:var])
subpart(graphopt, 3, :var)

# unroll the foldl
ACSets.ACSetInterface.collect_or_id(subpart(graphopt, 2, :tgt))

subpart(graphopt, 3, :var)
ACSets.ACSetInterface.collect_or_id(subpart(graphopt, 3, :var))

# fig 3.12 as AoN
df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:E,:F,:G,:end],
    Predecessor = [
        [], [:start], [:start], [:A], [:B], [:B], [:C,:D], [:D,:E], [:F,:G]
    ],
    Duration = zeros(9)
)
g = make_ProjGraph(df)
to_graphviz(g, node_labels=:label)