# from question: https://discourse.julialang.org/t/speeding-up-jump-model-creation-with-sets-that-depend-on-other-indexes/107333

using Catlab, JuMP
using BenchmarkTools
using MacroTools

@present NetSch(FreeSchema) begin
    (V,E)::Ob
    src::Hom(E,V)
    tgt::Hom(E,V)
    TimeType::AttrType
    te::Attr(E,TimeType)
    tv::Attr(V,TimeType)
end

# to_graphviz(NetSch, graph_attrs=Dict(:dpi=>"72",:size=>"3.5",:ratio=>"expand"))

@acset_type NetType(NetSch, index=[:src,:tgt])

num_vertices = 100

T1 = [1:3, 4:6, 7:9, 10:12]
T2 = [1:4, 5:8, 9:12]
T3 = [1:6, 7:12]
T_options = [T1, T2, T3]

# make the acset
const NetData = NetType{typeof(T1)}()

# add the verticies
const V = 1:num_vertices
add_parts!(
    NetData, :V, length(V),
    tv = [T_options[v%3+1] for v in V]
)

# add the edges
const E_list = [(i,j) for i in 1:num_vertices for j in i+1:num_vertices]
add_parts!(
    NetData, :E, length(E_list),
    src = first.(E_list), tgt = last.(E_list),
    te = [T_options[(a+b)%3+1] for (a,b) in E_list]
)

# dicts for the other method
const TV = Dict(
    v => T_options[v%3+1]
    for v in V
)
const TE = Dict(
    (a, b) => T_options[(a+b)%3+1]
    for (a, b) in E_list
)

# list of part IDs
const E = parts(NetData,:E)


# ------------------------------------------------------------
# comparision 1: make `flow` dv

test1() = begin
    model = JuMP.Model()
    @variable(
        model, 
        flow[e=E, NetData[e,:te]] ≥ 0
    )
end

test2() = begin
    model = JuMP.Model()
    @variable(model, flow[e ∈ E_list, TE[e]] ≥ 0)
end

@benchmark test1()
@benchmark test2()


# ------------------------------------------------------------
# comparision 2: make objective

test1() = begin
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
end

test2() = begin
    model = JuMP.Model()
    @variable(model, flow[e ∈ E_list, TE[e]] ≥ 0)
    @objective(model, Min,
        sum(flow[e, T] for e in E_list, T ∈ TE[e])
    )
end

@benchmark test1()
@benchmark test2()

# ------------------------------------------------------------
# comparision 3: make expression

test1() = begin
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
end

test2() = begin
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
end

@benchmark test1()
@benchmark test2()


# ------------------------------------------------------------
# the expression for `incoming` is the performance bottleneck.

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

MacroTools.prettify(
    @macroexpand @expression(
        model,
        incoming[v=V,Tᵥ=NetData[v,:tv]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[e,Tₑ]
            for e in incident(NetData,v,:tgt), Tₑ in NetData[e,:te]
        )
    )
)

# the results are side by side below



# compare to
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

MacroTools.prettify(
    @macroexpand @expression(model,
        incoming[v ∈ V, Tᵥ ∈ TV[v]],
        sum(
            length(Tₑ ∩ Tᵥ) * flow[(src, v), Tₑ]
            for src in [src for (src, dst) in E_list if dst == v], Tₑ ∈ TE[(src, v)]
        )
    )
)

# acs
dotterel = let model = model
        container(((v, Tᵥ)->begin
                    locust = let
                            manatee = MutableArithmetics.Zero()
                            for e = incident(NetData, v, :tgt)
                                for Tₑ = NetData[e, :te]
                                    squirrel = Tₑ ∩ Tᵥ
                                    sanddollar = length(squirrel)
                                    manatee = operate!!(add_mul, manatee, sanddollar, flow[e, Tₑ])
                                end
                            end
                            manatee
                        end
                    mouse = flatten!(locust)
                    JuMP._replace_zero(model, mouse)
                end), (JuMP.Containers).nested((()->V), ((v,)->NetData[v, :tv])), JuMP.Containers.AutoContainerType, Any[:v, :Tᵥ])
end

# dicts
dotterel = let model = model
        container(((v, Tᵥ)->begin
                    locust = let
                            manatee = MutableArithmetics.Zero()
                            for src = [src for (src, dst) = E_list if (NonlinearOperator(==, :==))(dst, v)]
                                for Tₑ = TE[(src, v)]
                                    squirrel = Tₑ ∩ Tᵥ
                                    sanddollar = length(squirrel)
                                    manatee = operate!!(add_mul, manatee, sanddollar, flow[(src, v), Tₑ])
                                end
                            end
                            manatee
                        end
                    mouse = flatten!(locust)
                    JuMP._replace_zero(model, mouse)
                end), (JuMP.Containers).nested((()->V), ((v,)->TV[v])), JuMP.Containers.AutoContainerType, Any[:v, :Tᵥ])
end

# benchmark incident versus this lookup method
test1(acs, vertices) = begin
    for v in vertices
        incident(acs,v,:tgt)        
    end
end

test2(E_list, vertices) = begin
    for v in vertices
        [src for (src, dst) in E_list if dst == v]
    end
end

test1(NetData, V)
test2(E_list, V)


# ------------------------------------------------------------
# ye olde stuff

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