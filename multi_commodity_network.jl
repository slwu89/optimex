using Catlab
using JuMP
import DataFrames
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



@present NetworkSch <: SchGraph begin
    NameType::AttrType
    NumberType::AttrType
    VectorType::AttrType
    VectorJumpType::AttrType # hold Vector{JuMP.VariableRef}
    # data
    name::Attr(V,NameType)
    c_ip::Attr(V,VectorType) # cost to purchase
    u_ip::Attr(V,VectorType) # supply capacity
    d_ip::Attr(V,VectorType) # demand capacity
    c_ijp::Attr(E,VectorType) # shipping cost
    u_ij::Attr(E,NumberType) # flow capacity
    # decision variables
    s_ip::Attr(V,VectorJumpType)
    x_ijp::Attr(E,VectorJumpType)
end

@acset_type NetworkData(NetworkSch, index=[:src,:tgt]) <: AbstractGraph

to_graphviz(NetworkSch, graph_attrs=Dict(:dpi=>"72",:size=>"4",:ratio=>"expand"))

