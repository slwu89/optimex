using Catlab
import Catlab.Graphics.Graphviz.Html
using JuMP
using DataFrames
import HiGHS
import SQLite
import SQLite.DBInterface


# --------------------------------------------------------------------------------
# load and transform data

filename = joinpath(@__DIR__, "commodity_nz.db");
db = SQLite.DB(filename)

function get_table(db, table)
    query = DBInterface.execute(db, "SELECT * FROM $table")
    return DataFrames.DataFrame(query)
end

df_shipping = get_table(db, "shipping")
df_products = get_table(db, "products")
df_supply = get_table(db, "supply")
df_demand = get_table(db, "demand")

df_cost = DataFrames.leftjoin(df_shipping, df_products; on = [:product])
df_cost.flow_cost = df_cost.cost_per_km .* df_cost.distance_km

places = unique([df_cost.origin; df_cost.destination])
products = unique(df_cost.product)

df_edges=unique(df_cost[:, [:origin, :destination]])
transform!(df_edges, :origin => ByRow(o -> begin
    findfirst(o .== places)
end) => :src)
transform!(df_edges, :destination => ByRow(o -> begin
    findfirst(o .== places)
end) => :tgt)

capacity=30

df_productvertex = crossjoin(DataFrame(place=places), DataFrame(product=products))
leftjoin!(df_productvertex, df_supply, on=[:place=>:origin, :product])
leftjoin!(df_productvertex, df_demand, on=[:place=>:destination, :product])
df_productvertex[ismissing.(df_productvertex.capacity), :capacity] .= 0.0
df_productvertex[ismissing.(df_productvertex.cost), :cost] .= 0.0
df_productvertex[ismissing.(df_productvertex.demand), :demand] .= 0.0


# --------------------------------------------------------------------------------
# define acset types and useful methods

@present SchMultiCommodity <: SchGraph begin    
    (ProdVertex,Product,Shipping)::Ob
    pv_p::Hom(ProdVertex,Product)
    pv_v::Hom(ProdVertex,V)
    s_p::Hom(Shipping,Product)
    s_e::Hom(Shipping,E)

    Label::AttrType
    vlabel::Attr(V,Label)
    plabel::Attr(Product,Label)

    NumVar::AttrType
    supplycap::Attr(ProdVertex,NumVar) # supply capacity
    demand::Attr(ProdVertex,NumVar) # demand
    purcost::Attr(ProdVertex,NumVar) # purchase cost
    shipcost::Attr(Shipping,NumVar) # shipment cost
    flowcap::Attr(E,NumVar) # flow capacity

    DecisionVar::AttrType
    s::Attr(ProdVertex,DecisionVar) # optimal supply
    x::Attr(Shipping,DecisionVar) # optimal shipment
end

to_graphviz(SchMultiCommodity, graph_attrs=Dict(:size=>"6",:ratio=>"fill"))

@abstract_acset_type AbstractMultiCommodity <: HasGraph

@acset_type MultiCommodity(SchMultiCommodity, index=[:src,:tgt,:pv_p,:pv_v,:s_p,:s_e]) <: AbstractMultiCommodity

"""
Make a node label for a model that has been optimized.
"""
function make_node_label_optim(v, acs, label_dict)
    prodv = incident(acs, v, :pv_v)
    label = Vector{String}()
    push!(
        label, """
            <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
            <TR><TD COLSPAN="2">$(acs[v, :vlabel])</TD></TR>
            <TR><TD>Product</TD><TD>Purchased</TD></TR>
        """
    )
    for pv in prodv
        if acs[pv, :s] > 0.0
            emoji = label_dict[acs[pv, (:pv_p, :plabel)]]
            push!(
                label, """
                    <TR><TD>$(emoji)</TD><TD>$(acs[pv, :s])</TD></TR>
                """
            )
        end
    end
    push!(
        label, """
            </TABLE>
        """
    )
end

"""
Make a node label for the model parameters.
"""
function make_node_label_structure(v, acs, label_dict)
    prodv = incident(acs, v, :pv_v)
    label = Vector{String}()
    push!(
        label, """
            <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
            <TR><TD COLSPAN="4">$(acs[v, :vlabel])</TD></TR>
            <TR><TD>Product</TD><TD>Supply Capacity</TD><TD>Demand</TD><TD>Purchase Cost</TD></TR>
        """
    )
    for pv in prodv
        emoji = label_dict[acs[pv, (:pv_p, :plabel)]]
        push!(
            label, """
                <TR><TD>$(emoji)</TD><TD>$(acs[pv, :supplycap])</TD><TD>$(acs[pv, :demand])</TD><TD>$(acs[pv, :purcost])</TD></TR>
            """
        )
    end
    push!(
        label, """
            </TABLE>
        """
    )
end

make_node_label_dict = Dict(:optim=>make_node_label_optim, :struct=>make_node_label_structure)

"""
Make an edge label for an optmized model.
"""
function make_edge_label_optim(e, acs, label_dict)
    ship = incident(acs, e, :s_e)
    label = Vector{String}()
    edge_label = acs[e, (:s_e, :src, :vlabel)] * " → " * acs[e, (:s_e, :tgt, :vlabel)]
    push!(
        label, """
            <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
            <TR><TD COLSPAN="2">$(edge_label)</TD></TR>
            <TR><TD>Product</TD><TD>Shipped</TD></TR>
        """
    )
    for s in ship
        if acs[s, :x] > 0.0
            emoji = label_dict[acs[s, (:s_p, :plabel)]]
            push!(
                label, """
                    <TR><TD>$(emoji)</TD><TD>$(acs[s, :x])</TD></TR>
                """
            )
        end
    end
    push!(
        label, """
            </TABLE>
        """
    )
end

"""
Make an edge label for the parameters of the model.
"""
function make_edge_label_structure(e, acs, label_dict)
    ship = incident(acs, e, :s_e)
    label = Vector{String}()
    edge_label = acs[e, (:s_e, :src, :vlabel)] * " → " * acs[e, (:s_e, :tgt, :vlabel)]
    push!(
        label, """
            <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
            <TR><TD COLSPAN="2">$(edge_label)</TD></TR>
            <TR><TD>Product</TD><TD>Shipping Cost</TD></TR>
        """
    )
    for s in ship
        emoji = label_dict[acs[s, (:s_p, :plabel)]]
        push!(
            label, """
                <TR><TD>$(emoji)</TD><TD>$(acs[s, :shipcost])</TD></TR>
            """
        )
    end
    push!(
        label, """
            </TABLE>
        """
    )
end

make_edge_label_dict = Dict(:optim=>make_edge_label_optim,:struct=>make_edge_label_structure)

function make_property_graph(acs::T; labels, kwargs...) where {T<:AbstractMultiCommodity}
    pg = to_graphviz_property_graph(acs; kwargs...)   
    label_dict = Dict("milk"=>"🥛", "kiwifruit"=>"🥝")
    for v in parts(acs, :V)
        label = make_node_label_dict[labels](v, acs, label_dict)
        set_vprops!(pg, v, label=Html(join(label)), shape="plain")
    end
    for e in parts(acs, :E)
        label = make_edge_label_dict[labels](e, acs, label_dict)
        set_eprops!(pg, e, label=Html(join(label)), shape="plain")
    end
    return pg
end


# --------------------------------------------------------------------------------
# build acset

multinet_ = @acset MultiCommodity{String,Float64,Union{Float64,JuMP.VariableRef}} begin
    V=length(places)
    vlabel=places
    Product=length(products)
    plabel=products
    E=nrow(df_edges)
    src=df_edges.src
    tgt=df_edges.tgt
    flowcap=capacity
end

for r in eachrow(df_productvertex)
    v=only(incident(multinet_, r.place, :vlabel))
    p=only(incident(multinet_, r.product, :plabel))
    add_part!(
        multinet_, :ProdVertex, 
        pv_v=v, pv_p=p,
        supplycap=r.capacity, demand=r.demand, purcost=r.cost
    )
end

for r in eachrow(df_cost)
    src=only(incident(multinet_, r.origin, :vlabel))
    tgt=only(incident(multinet_, r.destination, :vlabel))
    product=only(incident(multinet_, r.product, :plabel))
    add_part!(
        multinet_, :Shipping,
        s_p=product, s_e=only(edges(multinet_, src, tgt)),
        shipcost=r.flow_cost
    )
end

to_graphviz(make_property_graph(multinet_, labels=:struct, graph_attrs=Dict(:rankdir=>"TD")))


# --------------------------------------------------------------------------------
# make JuMP model

"""
Adds `s` and `x`, decision variables for how much supply to procure of a product at each node,
and how much of each product to ship along edges.
"""
function add_variables!(acs, model)
    acs[:, :s] = @variable(model, 0 ≤ s[i=parts(acs, :ProdVertex)] ≤ acs[i, :supplycap])
    acs[:, :x] = @variable(model, 0 ≤ x[i=parts(acs, :Shipping)])
end

"""
Add constraints to respect total flow capacity of all products along each edge.
"""
function add_flow_constraint!(acs, model)
    for i in parts(acs, :E)
        @constraint(
            model,
            sum(acs[incident(acs, i, :s_e), :x]) ≤ acs[i, :flowcap]
        )
    end
end

"""
Add constraint that, for each product, node pair, the supply plus inbound
shipping, minus outbound shipping, must equal the demand.
"""
function add_conservation_constraint!(acs, model)
    for i in parts(acs, :ProdVertex)
        v = acs[i, :pv_v]
        p = acs[i, :pv_p]

        # inbound shipping
        inbound = intersect(incident(acs, v, (:s_e, :tgt)), incident(acs, p, :s_p))

        # outbound shipping
        outbound = intersect(incident(acs, v, (:s_e, :src)), incident(acs, p, :s_p))

        @constraint(
            model,
            acs[i, :s] + 
                sum(acs[inbound, :x], init=zero(JuMP.VariableRef)) - 
                sum(acs[outbound, :x], init=zero(JuMP.VariableRef)) == acs[i, :demand]
        )    
    end
end

"""
Minimize costs of shipping and costs of purchase
"""
function add_objective!(acs, model)
    @objective(
        model, Min,
        sum(acs[:, :shipcost] .* acs[:, :x]) + 
            sum(acs[:, :purcost] .* acs[:, :s])
    )
end


# make the JuMP model
multinet1 = deepcopy(multinet_)
model = JuMP.Model(HiGHS.Optimizer)

add_variables!(multinet1, model)
add_flow_constraint!(multinet1, model)
add_conservation_constraint!(multinet1, model)
add_objective!(multinet1, model)

optimize!(model)
solution_summary(model)

multinet1[:, :x] = value.(multinet1[:, :x])
multinet1[:, :s] = value.(multinet1[:, :s])

to_graphviz(make_property_graph(multinet1, labels=:optim, graph_attrs=Dict(:rankdir=>"TD")))

"""
Format optimized decision variables for comparison. 
"""
function format_output(acs)
    df_solution_flow = DataFrame(tables(acs).Shipping)
    df_solution_flow = df_solution_flow[df_solution_flow.x .> 0.0, :]
    df_solution_flow.origin = acs[df_solution_flow.s_e, (:src, :vlabel)]
    df_solution_flow.destination = acs[df_solution_flow.s_e, (:tgt, :vlabel)]
    df_solution_flow.product = acs[df_solution_flow.s_p, :plabel]
    select!(df_solution_flow, Not([:s_p,:s_e,:shipcost]))
    rename!(df_solution_flow, Dict(:x => :x_flow))
    return df_solution_flow
end

format_output(multinet1)


# --------------------------------------------------------------------------------
# once more, with feeling

multinet2 = deepcopy(multinet_)

model = JuMP.Model(HiGHS.Optimizer)

add_variables!(multinet2, model)
add_flow_constraint!(multinet2, model)

inbound_uwd = @relation (pv=pv_id, v=place, p=prod, x=x) begin
    ProdVertex(_id=pv_id, pv_p=p_id, pv_v=v_id)
    Product(_id=p_id, plabel=prod)
    V(_id=v_id, vlabel=place)
    E(_id=e_id, tgt=v_id)
    Shipping(_id=ship_id, s_e=e_id, s_p=p_id, x=x)
end

query(multinet2, inbound_uwd)
query(multinet2, inbound_uwd, (place="auckland", ))

for i in parts(acs, :ProdVertex)
    v = acs[i, :pv_v]
    p = acs[i, :pv_p]

    # inbound shipping
    inbound = intersect(incident(acs, v, (:s_e, :tgt)), incident(acs, p, :s_p))

    # outbound shipping
    outbound = intersect(incident(acs, v, (:s_e, :src)), incident(acs, p, :s_p))

    @constraint(
        model,
        acs[i, :s] + 
            sum(acs[inbound, :x], init=zero(JuMP.VariableRef)) - 
            sum(acs[outbound, :x], init=zero(JuMP.VariableRef)) == acs[i, :demand]
    )    
end