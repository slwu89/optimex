using Catlab
using JuMP
using DataFrames
import HiGHS
import SQLite
import SQLite.DBInterface

filename = joinpath(@__DIR__, "commodity_nz.db");
db = SQLite.DB(filename)
SQLite.tables(db)

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

places = unique([df_shipping.origin; df_shipping.destination])
product = unique(df_shipping.product)
df_edges=unique(df_shipping[:, [:origin, :destination]])
transform!(df_edges, :origin => ByRow(o -> begin
    findfirst(o .== places)
end) => :src)
transform!(df_edges, :destination => ByRow(o -> begin
    findfirst(o .== places)
end) => :tgt)

capacity=30

multinet = @acset MultiCommodity{String,Float64,JuMP.VariableRef} begin
    V=length(places)
    vlabel=places
    Product=length(product)
    plabel=product
    E=nrow(df_edges)
    src=df_edges.src
    tgt=df_edges.tgt
    flowcap=capacity
end

df_prodvertex = outerjoin(df_supply, df_demand, on=[:origin=>:destination, :product])
replace!(df_prodvertex.capacity, missing=>0.0)
replace!(df_prodvertex.cost, missing=>0.0)
replace!(df_prodvertex.demand, missing=>0.0)

for r in eachrow(df_prodvertex)
    v=only(incident(multinet, r.origin, :vlabel))
    p=only(incident(multinet, r.product, :plabel))
    add_part!(
        multinet, :ProdVertex, 
        pv_v=v, pv_p=p,
        supplycap=r.capacity, demand=r.demand, purcost=r.cost
    )
end

for r in eachrow(df_cost)
    src=only(incident(multinet, r.origin, :vlabel))
    tgt=only(incident(multinet, r.destination, :vlabel))
    product=only(incident(multinet, r.product, :plabel))
    add_part!(
        multinet, :Shipping,
        s_p=product, s_e=only(edges(multinet, src, tgt)),
        shipcost=r.flow_cost
    )
end
