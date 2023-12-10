# from question: https://discourse.julialang.org/t/speeding-up-jump-model-creation-with-sets-that-depend-on-other-indexes/107333

using Catlab, JuMP
using BenchmarkTools

@present NetSch(FreeSchema) begin
    (V,E)::Ob
    src::Hom(E,V)
    tgt::Hom(E,V)
    TimeType::AttrType
    te::Attr(E,TimeType)
    tv::Attr(V,TimeType)
end

to_graphviz(NetSch, graph_attrs=Dict(:dpi=>"72",:size=>"3.5",:ratio=>"expand"))

@acset_type NetType(NetSch, index=[:src,:tgt])

num_vertices = 100

T1 = [1:3, 4:6, 7:9, 10:12]
T2 = [1:4, 5:8, 9:12]
T3 = [1:6, 7:12]
T_options = [T1, T2, T3]

# make the acset
NetData = NetType{typeof(T1)}()

# add the verticies
V = 1:num_vertices
add_parts!(
    NetData, :V, length(V),
    tv = [T_options[v%3+1] for v in V]
)

# add the edges
E_list = [(i,j) for i in 1:num_vertices for j in i+1:num_vertices]
add_parts!(
    NetData, :E, length(E_list),
    src = first.(E_list), tgt = last.(E_list),
    te = [T_options[(a+b)%3+1] for (a,b) in E_list]
)

# ------------------------------------------------------------
# comparisions

V = parts(NetData,:V)
E = parts(NetData,:E)

# build model from acset
test_acs() = begin
    model = JuMP.Model()
    @variable(
        model, 
        flow[e=E, NetData[e,:te]] ≥ 0
    )
    @objective(
        model, 
        Min,
        sum(flow[e,t] for e=E, t=NetData[e,:te])
    )
    @expression(
        model,
        incoming[v=V,Tᵥ=NetData[v,:tv]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[e,Tₑ]
            for e in incident(NetData,v,:tgt), Tₑ in NetData[e,:te]
        )
    )
    @expression(
        model,
        outgoing[v=V,Tᵥ=NetData[v,:tv]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[e,Tₑ]
            for e in incident(NetData,v,:src), Tₑ in NetData[e,:te]
        )
    )
    @constraint(
        model,
        balance[v=V,Tᵥ=NetData[v,:tv]],
        incoming[v,Tᵥ] == outgoing[v,Tᵥ]
    )
    return model
end

# no acs
TV = Dict(
    v => T_options[v%3+1]
    for v in V
)
TE = Dict(
    (a, b) => T_options[(a+b)%3+1]
    for (a, b) in E_list
)

test_noacs() = begin
    model = JuMP.Model()
    @variable(model, flow[e ∈ E_list, TE[e]] ≥ 0)
    @objective(model, Min,
        sum(flow[e, T] for e in E_list, T ∈ TE[e])
    )
    @expression(model,
        incoming[v ∈ V, Tᵥ ∈ TV[v]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[(src, v), Tₑ]
            for src in [src for (src, dst) in E_list if dst == v], Tₑ ∈ TE[(src, v)]
        )
    )
    @expression(model,
        outgoing[v ∈ V, Tᵥ ∈ TV[v]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[(v, dst), Tₑ]
            for dst in [dst for (src, dst) in E_list if src == v], Tₑ ∈ TE[(v, dst)]
        )
    )
    @constraint(model,
        balance[v ∈ V, Tᵥ ∈ TV[v]],
        incoming[v, Tᵥ] == outgoing[v, Tᵥ]
    )
    return model
end

@benchmark test_acs()
@benchmark test_noacs()